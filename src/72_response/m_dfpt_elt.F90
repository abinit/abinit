!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dfpt_elt
!! NAME
!!  m_dfpt_elt
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (DRH, DCA, XG, GM, AR, MB)
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

module m_dfpt_elt

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_dtset

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes, only : MPI_type
 use m_time,        only : timab
 use m_special_funcs,  only : abi_derfc
 use m_geometry,    only : metric
 use m_cgtools,     only : dotprod_vn
 use m_pawtab,      only : pawtab_type,pawtab_free,pawtab_nullify
 use m_pawrad,      only : pawrad_type,pawrad_init,pawrad_free
 use m_pawpsp,      only : pawpsp_cg
 use m_paw_numeric, only : paw_spline
 use m_spacepar,    only : redgr
 use m_atm2fft,     only : atm2fft, dfpt_atm2fft
 use m_mkcore,      only : dfpt_mkcore
 use m_dfpt_mkvxcstr, only : dfpt_mkvxcstr
 use m_paral_atom, only : get_my_atmtab, free_my_atmtab
 use m_mpinfo,  only : ptabs_fourdp, proc_distrb_cycle
 use m_fftcore,      only : sphereboundary
 use m_fft,             only : fourdp

 implicit none

 private
!!***

 public :: dfpt_eltfrxc
 public :: dfpt_eltfrloc
 public :: dfpt_eltfrkin
 public :: dfpt_eltfrhar
 public :: elt_ewald
 public :: dfpt_ewald
!!***

contains
!!***

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
!!     see ~abinit/doc/variables/vargs.htm#ngfft
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

subroutine dfpt_eltfrxc(atindx,dtset,eltfrxc,enxc,gsqcut,kxc,mpi_enreg,mgfft,&
& nattyp,nfft,ngfft,ngfftf,nhat,nkxc,n3xccc,pawtab,ph1d,psps,rhor,rprimd,&
& usexcnhat,vxc,xccc3d,xred)

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
 logical :: nmxc
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
 nmxc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)

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
   call redgr (work,workgr,mpi_enreg,nfft,ngfft)
   do ifft=1,nfft
     rho0_redgr(:,ifft,1)=workgr(ifft,:)
   end do
   if(dtset%nspden==2) then
     work(:)=2.0_dp*(rhor_(:,1)-rhor_(:,2))
     if(n1xccc/=0) then
       work(:)=work(:)+xccc3d(:)
     end if
     call redgr(work,workgr,mpi_enreg,nfft,ngfft)
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
&     dummy_in,nkxc,nmxc,dtset%nspden,n3xccc_loc,option,qphon,rhor,rhor,&
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
&     dummy_in,nkxc,nmxc,dtset%nspden,n3xccc_loc,option,qphon,rhor,rhor,&
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
       call fourdp(1,vxc10_coreg,vxc10_core,-1,mpi_enreg,nfft,1, ngfft,0)
       call fourdp(1,vxc_coreg,vxc_core,-1,mpi_enreg,nfft,1, ngfft, 0)
       call fourdp(1,vxc1is_coreg,vxc1is_core,-1,mpi_enreg,nfft,1, ngfft, 0)

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
         call fourdp(1,vxc10_coreg,vxc10_core,-1,mpi_enreg,nfft,1, ngfft, 0)
         call fourdp(1,vxc_coreg,vxc_core,-1,mpi_enreg,nfft,1, ngfft,0)
         call fourdp(1,vxc1is_coreg,vxc1is_core,-1,mpi_enreg,nfft,1, ngfft, 0)
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

!!****f* ABINIT/eltxccore
!! NAME
!! eltxccore
!!
!! FUNCTION
!! Compute the core charge contributions to the 2nd derivatives
!! of the exchange-correlation energy with respect to all pairs of
!! strain or strain and atomic displacement for the frozen wavefunction
!! contribution to the elastic tensor. 1st-order potentials representing
!! the perturbation by one strain are supplied, and the routine loops
!! over the second strain and over all atomic displacements.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (DRH, DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  nfft=number of fft grid points
!!  ntypat=number of types of atoms in cell.
!!  n1,n2,n3=fft grid dimensions.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  typat(natom)=integer type for each atom in cell
!!  ucvol=unit cell volume (bohr**3).
!!  vxc_core(nfft)=spin-averaged xc potential
!!  vxc10_core(nfft)=spin-averaged 1st-order xc potential for elastic tensor
!!  vxc1is_core(nfft)=spin-averaged 1st-order xc potential for internal strain
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc,6,ntypat)=1D core charge function and five derivatives,
!!   for each type of atom, from psp
!!  xred(3,natom)=reduced coordinates for atoms in unit cell
!!
!! OUTPUT
!!  eltfrxc(6+3*natom,6) = xc frozen wavefunction contribution to the
!!   elastic tensor
!!
!! SIDE EFFECTS
!!  eltfrxc(6+3*natom,6) = xc frozen wavefunction contribution to the
!!   elastic and internal-strain tensor.  One column is incremented
!!   by the core contribution.
!!
!! NOTES
!! Note that this routine is related to the mkcore.f routine
!!
!! PARENTS
!!      dfpt_eltfrxc
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,timab,xmpi_sum
!!
!! SOURCE

subroutine eltxccore(eltfrxc,is2_in,my_natom,natom,nfft,ntypat,&
& n1,n1xccc,n2,n3,rprimd,typat,ucvol,vxc_core,vxc10_core,vxc1is_core,&
& xcccrc,xccc1d,xred, &
& mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: is2_in,n1,n1xccc,n2,n3,my_natom,natom,nfft,ntypat
 integer,optional,intent(in) :: comm_atom
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: rprimd(3,3),vxc10_core(nfft),vxc1is_core(nfft)
 real(dp),intent(in) :: vxc_core(nfft),xccc1d(n1xccc,6,ntypat)
 real(dp),intent(in) :: xcccrc(ntypat),xred(3,natom)
 real(dp),intent(inout) :: eltfrxc(6+3*natom,6)

!Local variables-------------------------------
!scalars
 integer,parameter :: mshift=401
 integer :: i1,i2,i3,iat,iatom,ierr,ifft,is1,is2,ishift,ishift1,ishift2
 integer :: ishift3,ixp,jj,js,ka,kb,kd,kg,mu,my_comm_atom,nu
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: aa,bb,cc,d2rss,dd,delta,delta2div6,deltam1,diff
 real(dp) :: difmag,difmag2,difmag2_fact,difmag2_part,drss1,drss2,func1
 real(dp) :: func2,range,range2,rangem1,rdiff1,rdiff2,rdiff3
 real(dp) :: term1,term2,yy
 character(len=500) :: message
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer :: igrid(3),ii(mshift,3),irange(3),ngfft(3)
 integer,pointer :: my_atmtab(:)
 real(dp) :: drm(3,3,6),eltfrxc_core(6+3*natom,6),lencp(3),rmet(3,3),rrdiff(mshift,3)
 real(dp) :: scale(3),tau(3),ts2(3),tsec(2),tt(3)
 real(dp),allocatable :: d2rm(:,:,:,:)

! *************************************************************************

!Compute lengths of cross products for pairs of primitive
!translation vectors (used in setting index search range below)
 lencp(1)=cross_elt(rprimd(1,2),rprimd(2,2),rprimd(3,2),&
& rprimd(1,3),rprimd(2,3),rprimd(3,3))
 lencp(2)=cross_elt(rprimd(1,3),rprimd(2,3),rprimd(3,3),&
& rprimd(1,1),rprimd(2,1),rprimd(3,1))
 lencp(3)=cross_elt(rprimd(1,1),rprimd(2,1),rprimd(3,1),&
& rprimd(1,2),rprimd(2,2),rprimd(3,2))

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Compute factor R1.(R2xR3)/|R2xR3| etc for 1, 2, 3
!(recall ucvol=R1.(R2xR3))
 scale(:)=ucvol/lencp(:)

!Compute metric tensor in real space rmet
 do nu=1,3
   rmet(:,nu)=rprimd(1,:)*rprimd(1,nu)+rprimd(2,:)*rprimd(2,nu)+&
&   rprimd(3,:)*rprimd(3,nu)
 end do

!Compute 1st and 2nd derivatives of metric tensor wrt all strain components
!and store for use in inner loop below.

 ABI_ALLOCATE(d2rm,(3,3,6,6))

!Loop over 2nd strain index
 do is2=1,6
   kg=idx(2*is2-1);kd=idx(2*is2)
   do jj = 1,3
     drm(:,jj,is2)=rprimd(kg,:)*rprimd(kd,jj)+rprimd(kd,:)*rprimd(kg,jj)
   end do

!  Loop over 1st strain index
   do is1=1,6

     ka=idx(2*is1-1);kb=idx(2*is1)
     d2rm(:,:,is1,is2)=0._dp
     do jj = 1,3
       if(ka==kg) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(kb,:)*rprimd(kd,jj)+rprimd(kd,:)*rprimd(kb,jj)
       if(ka==kd) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(kb,:)*rprimd(kg,jj)+rprimd(kg,:)*rprimd(kb,jj)
       if(kb==kg) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(ka,:)*rprimd(kd,jj)+rprimd(kd,:)*rprimd(ka,jj)
       if(kb==kd) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(ka,:)*rprimd(kg,jj)+rprimd(kg,:)*rprimd(ka,jj)
     end do
   end do !is1
 end do !is2

 ngfft(1)=n1
 ngfft(2)=n2
 ngfft(3)=n3
 delta=1.0_dp/(n1xccc-1)
 deltam1=n1xccc-1
 delta2div6=delta**2/6.0_dp

!Loop over atoms in unit cell
 eltfrxc_core(:,:)=zero

 do iat=1,my_natom
   iatom=iat;if (paral_atom) iatom=my_atmtab(iat)
   js=7+3*(iatom-1)
!  Set search range (density cuts off perfectly beyond range)
!  Cycle if no range.
   range=0.0_dp
   range=xcccrc(typat(iatom))
   if(range<1.d-16) cycle

   range2=range**2
   rangem1=1.0_dp/range

!  Consider each component in turn
   do mu=1,3
     tau(mu)=mod(xred(mu,iatom)+1._dp-aint(xred(mu,iatom)),1._dp)

!    Use tau to find nearest grid point along R(mu)
!    (igrid=0 is the origin; shift by 1 to agree with usual index)
     igrid(mu)=nint(tau(mu)*dble(ngfft(mu)))

!    Use range to compute an index range along R(mu)
!    (add 1 to make sure it covers full range)
     irange(mu)=1+nint((range/scale(mu))*dble(ngfft(mu)))

!    Check that the largest range is smallest than the maximum
!    allowed one
     if(2*irange(mu)+1 > mshift)then
       write(message, '(a,i0,a)' )' The range around atom',iatom,' is too large.'
       MSG_BUG(message)
     end if

!    Set up a counter that explore the relevant range
!    of points around the atom
     ishift=0
     do ixp=igrid(mu)-irange(mu),igrid(mu)+irange(mu)
       ishift=ishift+1
       ii(ishift,mu)=1+mod(ngfft(mu)+mod(ixp,ngfft(mu)),ngfft(mu))
       rrdiff(ishift,mu)=dble(ixp)/dble(ngfft(mu))-tau(mu)
     end do

!    End loop on mu
   end do

!  Conduct triple loop over restricted range of grid points for iatom

   do ishift3=1,1+2*irange(3)
!    map back to [1,ngfft(3)] for usual fortran index in unit cell
     i3=ii(ishift3,3)
!    find vector from atom location to grid point (reduced)
     rdiff3=rrdiff(ishift3,3)

     do ishift2=1,1+2*irange(2)
       i2=ii(ishift2,2)
       rdiff2=rrdiff(ishift2,2)
!      Prepare the computation of difmag2
       difmag2_part=rmet(3,3)*rdiff3**2+rmet(2,2)*rdiff2**2&
&       +2.0_dp*rmet(3,2)*rdiff3*rdiff2
       difmag2_fact=2.0_dp*(rmet(3,1)*rdiff3+rmet(2,1)*rdiff2)

       do ishift1=1,1+2*irange(1)
         rdiff1=rrdiff(ishift1,1)

!        Compute (rgrid-tau-Rprim)**2
         difmag2= difmag2_part+rdiff1*(difmag2_fact+rmet(1,1)*rdiff1)

!        Only accept contribution inside defined range
         if (difmag2<range2) then

!          Prepare computation of core charge function and derivatives,
!          using splines
           difmag=sqrt(difmag2)
           if (difmag>=1.0d-10) then
             i1=ii(ishift1,1)
             yy=difmag*rangem1

!            Compute index of yy over 1 to n1xccc scale
             jj=1+int(yy*(n1xccc-1))
             diff=yy-(jj-1)*delta

!            Will evaluate spline fit (p. 86 Numerical Recipes, Press et al;
!            NOTE error in book for sign of "aa" term in derivative;
!            also see splfit routine).
             bb = diff*deltam1
             aa = 1.0_dp-bb
             cc = aa*(aa**2-1.0_dp)*delta2div6
             dd = bb*(bb**2-1.0_dp)*delta2div6

!            Evaluate spline fit of 1st der of core charge density
!            from xccc1d(:,2,:) and (:,4,:)
             func1=aa*xccc1d(jj,2,typat(iatom))+bb*xccc1d(jj+1,2,typat(iatom)) +&
&             cc*xccc1d(jj,4,typat(iatom))+dd*xccc1d(jj+1,4,typat(iatom))
             term1=func1*rangem1
!            Evaluate spline fit of 2nd der of core charge density
!            from xccc1d(:,3,:) and (:,5,:)
             func2=aa*xccc1d(jj,3,typat(iatom))+bb*xccc1d(jj+1,3,typat(iatom)) +&
&             cc*xccc1d(jj,5,typat(iatom))+dd*xccc1d(jj+1,5,typat(iatom))
             term2=func2*rangem1**2

             ifft=i1+n1*(i2-1+n2*(i3-1))
             tt(:)=rmet(:,1)*rdiff1+rmet(:,2)*rdiff2+rmet(:,3)*rdiff3

!            Add contributions to 2nd derivative tensor
             drss2=&
&             (rdiff1*(drm(1,1,is2_in)*rdiff1+drm(1,2,is2_in)*rdiff2&
&             +drm(1,3,is2_in)*rdiff3)&
&             +rdiff2*(drm(2,1,is2_in)*rdiff1+drm(2,2,is2_in)*rdiff2&
&             +drm(2,3,is2_in)*rdiff3)&
&             +rdiff3*(drm(3,1,is2_in)*rdiff1+drm(3,2,is2_in)*rdiff2&
&             +drm(3,3,is2_in)*rdiff3))

!            Loop over 1st strain index
             do is1=1,6

               drss1=&
&               (rdiff1*(drm(1,1,is1)*rdiff1+drm(1,2,is1)*rdiff2&
&               +drm(1,3,is1)*rdiff3)&
&               +rdiff2*(drm(2,1,is1)*rdiff1+drm(2,2,is1)*rdiff2&
&               +drm(2,3,is1)*rdiff3)&
&               +rdiff3*(drm(3,1,is1)*rdiff1+drm(3,2,is1)*rdiff2&
&               +drm(3,3,is1)*rdiff3))

               d2rss=&
&               (rdiff1*(d2rm(1,1,is1,is2_in)*rdiff1+d2rm(1,2,is1,is2_in)*rdiff2&
&               +d2rm(1,3,is1,is2_in)*rdiff3)&
&               +rdiff2*(d2rm(2,1,is1,is2_in)*rdiff1+d2rm(2,2,is1,is2_in)*rdiff2&
&               +d2rm(2,3,is1,is2_in)*rdiff3)&
&               +rdiff3*(d2rm(3,1,is1,is2_in)*rdiff1+d2rm(3,2,is1,is2_in)*rdiff2&
&               +d2rm(3,3,is1,is2_in)*rdiff3))

!              Vall(0) X Rhocore(2) term
               eltfrxc_core(is1,is2_in)=eltfrxc_core(is1,is2_in)+0.25_dp*&
&               (vxc_core(ifft)*(term1*(d2rss/difmag&
&               -drss1*drss2/difmag**3)&
&               +term2*drss1*drss2/difmag**2))

!              Vall(1) X Rhocore(1) term
               eltfrxc_core(is1,is2_in)=eltfrxc_core(is1,is2_in)+0.25_dp*&
&               vxc10_core(ifft)*drss1*term1/difmag
               eltfrxc_core(is2_in,is1)=eltfrxc_core(is2_in,is1)+0.25_dp*&
&               vxc10_core(ifft)*drss1*term1/difmag

!              End loop in is1
             end do
!            Internal strain contributions
             ts2(:)=drm(:,1,is2_in)*rdiff1+drm(:,2,is2_in)*rdiff2&
&             +drm(:,3,is2_in)*rdiff3

             eltfrxc_core(js:js+2,is2_in)=eltfrxc_core(js:js+2,is2_in)&
&             -(vxc1is_core(ifft)*term1/difmag&
&             +0.5_dp*vxc_core(ifft)*(term2-term1/difmag)*drss2/difmag**2)*tt(:)&
&             -(vxc_core(ifft)*term1/difmag)*ts2(:)

!            End of the condition for the distance not to vanish
           end if

!          End of condition to be inside the range
         end if

!        End loop on ishift1
       end do

!      End loop on ishift2
     end do

!    End loop on ishift3
   end do

!  End loop on atoms
 end do

!In case of parallelism over atoms: communicate
 if (paral_atom) then
   call timab(48,1,tsec)
   call xmpi_sum(eltfrxc_core,my_comm_atom,ierr)
   call timab(48,2,tsec)
 end if

!Add core contribution to XC elastic tensor
 eltfrxc(:,:)=eltfrxc(:,:)+eltfrxc_core(:,:)

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 ABI_DEALLOCATE(d2rm)

 contains

   function cross_elt(xx,yy,zz,aa,bb,cc)
!Define magnitude of cross product of two vectors
   real(dp) :: cross_elt
   real(dp),intent(in) :: xx,yy,zz,aa,bb,cc
   cross_elt=sqrt((yy*cc-zz*bb)**2+(zz*aa-xx*cc)**2+(xx*bb-yy*aa)**2)
 end function cross_elt

end subroutine eltxccore
!!***

!!****f* ABINIT/dfpt_eltfrloc
!! NAME
!! dfpt_eltfrloc
!!
!! FUNCTION
!! Compute the frozen-wavefunction local pseudopotential contribution
!! to the elastic tensor and the internal strain (derivative wrt one
!! cartesian strain component and one reduced-coordinate atomic displacement).
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space (bohr**-1)
!!  gsqcut=cutoff on G^2 based on ecut
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  mqgrid=dimensioned number of q grid points for local psp spline
!!  natom=number of atoms in unit cell
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/variables/vargs.htm#ngfft
!!  ntypat=number of types of atoms
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  qgrid(mqgrid)=q point array for local psp spline fits
!!  rhog(2,nfft)=electron density in G space
!!  vlspl(mqgrid,2,ntypat)=q^2 v(q) spline for each type of atom.
!!
!! OUTPUT
!!  eltfrloc(6+3*natom,6)=non-symmetrized local pseudopotenial contribution
!!   to the elastic tensor and internal strain.
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine dfpt_eltfrloc(atindx,eltfrloc,gmet,gprimd,gsqcut,mgfft,&
&  mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,ntypat,ph1d,qgrid,rhog,vlspl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,mqgrid,natom,nfft,ntypat
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),rhog(2,nfft),vlspl(mqgrid,2,ntypat)
 real(dp),intent(out) :: eltfrloc(6+3*natom,6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,ia,ia1,ia2,id1,id2,id3,ielt,ieltx,ierr,ig1,ig2,ig3,ii
 integer :: is1,is2,itypat,jj,ka,kb,kd,kg,me_fft,n1,n2,n3,nproc_fft
 real(dp),parameter :: tolfix=1.0000001_dp
 real(dp) :: aa,bb,cc,cutoff,d2g
 real(dp) :: dd,dg1,dg2,diff,dq
 real(dp) :: dq2div6,dqdiv6,dqm1,ee,ff,gmag,gsquar
 real(dp) :: sfi,sfr,term,term1
!real(dp) :: ph1_elt,ph2_elt,ph3_elt,phi_elt,phr_elt
 real(dp) :: term2,term3,term4,term5,vion1,vion2,vion3
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: dgm(3,3,6),tsec(2)
 real(dp),allocatable :: d2gm(:,:,:,:),elt_work(:,:)

! *************************************************************************

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 me_fft=ngfft(11) ; nproc_fft=ngfft(10)

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!-----
!Compute 1st and 2nd derivatives of metric tensor wrt all strain components
!and store for use in inner loop below.
 ABI_ALLOCATE(d2gm,(3,3,6,6))

!Loop over 2nd strain index
 do is2=1,6
   kg=idx(2*is2-1);kd=idx(2*is2)
   do jj = 1,3
     dgm(:,jj,is2)=-(gprimd(kg,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(kg,jj))
   end do
!  Loop over 1st strain index, upper triangle only
   do is1=1,is2
     ka=idx(2*is1-1);kb=idx(2*is1)
     d2gm(:,:,is1,is2)=0._dp
     do jj = 1,3
       if(ka==kg) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(kb,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(kb,jj)
       if(ka==kd) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(kb,:)*gprimd(kg,jj)+gprimd(kg,:)*gprimd(kb,jj)
       if(kb==kg) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(ka,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(ka,jj)
       if(kb==kd) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(ka,:)*gprimd(kg,jj)+gprimd(kg,:)*gprimd(ka,jj)
     end do
   end do !is1
 end do !is2

!Zero out array to permit accumulation over atom types below:
 eltfrloc(:,:)=0.0_dp

 dq=(qgrid(mqgrid)-qgrid(1))/dble(mqgrid-1)
 dqm1=1.0_dp/dq
 dqdiv6=dq/6.0_dp
 dq2div6=dq**2/6.0_dp
 cutoff=gsqcut*tolfix
 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2

 ia1=1
 do itypat=1,ntypat
!  ia1,ia2 sets range of loop over atoms:
   ia2=ia1+nattyp(itypat)-1
   ii=0
   do i3=1,n3
     ig3=i3-(i3/id3)*n3-1
     do i2=1,n2
       if (fftn2_distrib(i2)==me_fft) then
         ig2=i2-(i2/id2)*n2-1
         do i1=1,n1
           ig1=i1-(i1/id1)*n1-1

           ii=ii+1
!          Skip G=0:
           if (ig1==0 .and. ig2==0 .and. ig3==0) cycle

!          Skip G**2 outside cutoff:
           gsquar=gsq_elt(ig1,ig2,ig3)
           if (gsquar<=cutoff) then
             gmag=sqrt(gsquar)

!            Compute vion(G) for given type of atom
             jj=1+int(gmag*dqm1)
             diff=gmag-qgrid(jj)

!            Evaluate spline fit from q^2 V(q) to get V(q):
!            (p. 86 Numerical Recipes, Press et al; NOTE error in book for sign
!             of "aa" term in derivative; also see splfit routine).
             bb = diff*dqm1
             aa = 1.0_dp-bb
             cc = aa*(aa**2-1.0_dp)*dq2div6
             dd = bb*(bb**2-1.0_dp)*dq2div6
             term1 = (aa*vlspl(jj,1,itypat)+bb*vlspl(jj+1,1,itypat) +&
&             cc*vlspl(jj,2,itypat)+dd*vlspl(jj+1,2,itypat))
             vion1=term1 / gsquar

!            Also get dV(q)/dq:
!            (note correction of Numerical Recipes sign error
!             before (3._dp*aa**2-1._dp)
             ee= vlspl(jj+1,1,itypat)-vlspl(jj,1,itypat)
             ff=  (3._dp*bb**2-1._dp)*vlspl(jj+1,2,itypat) &
&             - (3._dp*aa**2-1._dp)*vlspl(jj,2,itypat)
             term2 = ee*dqm1 + ff*dqdiv6
             vion2 = term2/gsquar - 2._dp*term1/(gsquar*gmag)

!            Also get V''(q)
             term3=aa*vlspl(jj,2,itypat)+bb*vlspl(jj+1,2,itypat)
             vion3 = (term3 - 4.0_dp*term2/gmag + 6._dp*term1/gsquar)/gsquar

!            Assemble structure factor over all atoms of given type:
             sfr=zero;sfi=zero
             do ia=ia1,ia2
               sfr=sfr+phre_elt(ig1,ig2,ig3,ia)
               sfi=sfi-phimag_elt(ig1,ig2,ig3,ia)
             end do
             term=(rhog(re,ii)*sfr+rhog(im,ii)*sfi)

!            Loop over 2nd strain index
             do is2=1,6
               dg2=0.5_dp*dgsqds_elt(ig1,ig2,ig3,is2)/gmag
!              Loop over 1st strain index, upper triangle only
               do is1=1,is2
                 dg1=0.5_dp*dgsqds_elt(ig1,ig2,ig3,is1)/gmag
                 d2g=(0.25_dp*d2gsqds_elt(ig1,ig2,ig3,is1,is2)-dg1*dg2)/gmag
                 eltfrloc(is1,is2)=eltfrloc(is1,is2)+&
&                 term*(vion3*dg1*dg2+vion2*d2g)
                 if(is2<=3)&
&                 eltfrloc(is1,is2)=eltfrloc(is1,is2)-term*vion2*dg1
                 if(is1<=3)&
&                 eltfrloc(is1,is2)=eltfrloc(is1,is2)-term*vion2*dg2
                 if(is1<=3 .and. is2<=3)&
&                 eltfrloc(is1,is2)=eltfrloc(is1,is2)+term*vion1
               end do !is1

!              Internal strain section - loop over current atoms
               do ia=ia1,ia2
                 if(is2 <=3) then
                   term4=vion2*dg2-vion1
                 else
                   term4=vion2*dg2
                 end if
                 term5=-two_pi*(rhog(re,ii)*phimag_elt(ig1,ig2,ig3,ia)&
&                 +rhog(im,ii)*phre_elt(ig1,ig2,ig3,ia))*term4
                 eltfrloc(7+3*(ia-1),is2)=eltfrloc(7+3*(ia-1),is2)+term5*dble(ig1)
                 eltfrloc(8+3*(ia-1),is2)=eltfrloc(8+3*(ia-1),is2)+term5*dble(ig2)
                 eltfrloc(9+3*(ia-1),is2)=eltfrloc(9+3*(ia-1),is2)+term5*dble(ig3)
               end do

             end do !is2

!            End skip G**2 outside cutoff:
           end if

!          End loop on n1, n2, n3. There is a "cycle" inside the loop
         end do
       end if
     end do
   end do

!  End loop on type of atoms
   ia1=ia2+1
 end do
!Init mpi_comm
 call timab(48,1,tsec)
 call xmpi_sum(eltfrloc,mpi_enreg%comm_fft,ierr)
 call timab(48,2,tsec)

!Fill in lower triangle
 do is2=2,6
   do is1=1,is2-1
     eltfrloc(is2,is1)=eltfrloc(is1,is2)
   end do
 end do

!The indexing array atindx is used to reestablish the correct
!order of atoms
 ABI_ALLOCATE(elt_work,(6+3*natom,6))
 elt_work(1:6,1:6)=eltfrloc(1:6,1:6)
 do ia=1,natom
   ielt=7+3*(ia-1)
   ieltx=7+3*(atindx(ia)-1)
   elt_work(ielt:ielt+2,1:6)=eltfrloc(ieltx:ieltx+2,1:6)
 end do
 eltfrloc(:,:)=elt_work(:,:)

 ABI_DEALLOCATE(d2gm)
 ABI_DEALLOCATE(elt_work)

 contains

!Real and imaginary parts of phase.
   function phr_elt(x1,y1,x2,y2,x3,y3)

   real(dp) :: phr_elt,x1,x2,x3,y1,y2,y3
   phr_elt=(x1*x2-y1*y2)*x3-(y1*x2+x1*y2)*y3
 end function phr_elt

   function phi_elt(x1,y1,x2,y2,x3,y3)

   real(dp):: phi_elt,x1,x2,x3,y1,y2,y3
   phi_elt=(x1*x2-y1*y2)*y3+(y1*x2+x1*y2)*x3
 end function phi_elt

   function ph1_elt(nri,ig1,ia)

   real(dp):: ph1_elt
   integer :: nri,ig1,ia
   ph1_elt=ph1d(nri,ig1+1+n1+(ia-1)*(2*n1+1))
 end function ph1_elt

   function ph2_elt(nri,ig2,ia)

   real(dp):: ph2_elt
   integer :: nri,ig2,ia
   ph2_elt=ph1d(nri,ig2+1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1))
 end function ph2_elt

   function ph3_elt(nri,ig3,ia)

   real(dp):: ph3_elt
   integer :: nri,ig3,ia
   ph3_elt=ph1d(nri,ig3+1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1))
 end function ph3_elt

   function phre_elt(ig1,ig2,ig3,ia)

   real(dp):: phre_elt
   integer :: ig1,ig2,ig3,ia
   phre_elt=phr_elt(ph1_elt(re,ig1,ia),ph1_elt(im,ig1,ia),&
&   ph2_elt(re,ig2,ia),ph2_elt(im,ig2,ia),ph3_elt(re,ig3,ia),ph3_elt(im,ig3,ia))
 end function phre_elt

   function phimag_elt(ig1,ig2,ig3,ia)

   real(dp) :: phimag_elt
   integer :: ig1,ig2,ig3,ia
   phimag_elt=phi_elt(ph1_elt(re,ig1,ia),ph1_elt(im,ig1,ia),&
&   ph2_elt(re,ig2,ia),ph2_elt(im,ig2,ia),ph3_elt(re,ig3,ia),ph3_elt(im,ig3,ia))
 end function phimag_elt

   function gsq_elt(i1,i2,i3)

   real(dp) :: gsq_elt
   integer :: i1,i2,i3
!Define G^2 based on G space metric gmet.
   gsq_elt=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+&
&   dble(i3*i3)*gmet(3,3)+dble(2*i1*i2)*gmet(1,2)+&
&   dble(2*i2*i3)*gmet(2,3)+dble(2*i3*i1)*gmet(3,1)
 end function gsq_elt

   function dgsqds_elt(i1,i2,i3,is)

   real(dp) :: dgsqds_elt
   integer :: i1,i2,i3,is
!Define dG^2/ds based on G space metric derivative
   dgsqds_elt=dble(i1*i1)*dgm(1,1,is)+dble(i2*i2)*dgm(2,2,is)+&
&   dble(i3*i3)*dgm(3,3,is)+&
&   dble(i1*i2)*(dgm(1,2,is)+dgm(2,1,is))+&
&   dble(i1*i3)*(dgm(1,3,is)+dgm(3,1,is))+&
&   dble(i2*i3)*(dgm(2,3,is)+dgm(3,2,is))
 end function dgsqds_elt

   function d2gsqds_elt(i1,i2,i3,is1,is2)

   real(dp) :: d2gsqds_elt
   integer :: i1,i2,i3,is1,is2
!Define 2dG^2/ds1ds2  based on G space metric derivative
   d2gsqds_elt=dble(i1*i1)*d2gm(1,1,is1,is2)+&
&   dble(i2*i2)*d2gm(2,2,is1,is2)+dble(i3*i3)*d2gm(3,3,is1,is2)+&
&   dble(i1*i2)*(d2gm(1,2,is1,is2)+d2gm(2,1,is1,is2))+&
&   dble(i1*i3)*(d2gm(1,3,is1,is2)+d2gm(3,1,is1,is2))+&
&   dble(i2*i3)*(d2gm(2,3,is1,is2)+d2gm(3,2,is1,is2))
 end function d2gsqds_elt

end subroutine dfpt_eltfrloc
!!***

!!****f* ABINIT/dfpt_eltfrkin
!! NAME
!! dfpt_eltfrkin
!!
!! FUNCTION
!! Compute the frozen-wavefunction kinetic enegy contribution to the
!! elastic tensor
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=<G|Cnk>=Fourier coefficients of wavefunction
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha) (NOT NEEDED !)
!!  effmass_free=effective mass for electrons (1. in common case)
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=work array for coordinates of G vectors in basis
!!  kptns(3,nkpt)=coordinates of k points in terms of reciprocal space
!!   primitive translations
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimension for number of planewaves
!!  nband(nkpt*nsppol)=number of bands being considered per k point
!!  nkpt=number of k points
!!  ngfft(18)=contain all needed information about 3D FFT, i
!!    see ~abinit/doc/variables/vargs.htm#ngfft
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for polarized
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2)
!!    at each k point
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  wtk(nkpt)=k point weights
!!
!! OUTPUT
!!  eltfrkin(6,6)=non-symmetrized kinetic energy contribution to the
!!                    elastic tensor
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      d2kindstr2,metric,sphereboundary,timab,xmpi_sum
!!
!! SOURCE

subroutine dfpt_eltfrkin(cg,eltfrkin,ecut,ecutsm,effmass_free,&
&  istwfk,kg,kptns,mband,mgfft,mkmem,mpi_enreg,&
&  mpw,nband,nkpt,ngfft,npwarr,nspinor,nsppol,occ,rprimd,wtk)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mgfft,mkmem,mpw,nkpt,nspinor,nsppol
 real(dp),intent(in) :: ecut,ecutsm,effmass_free
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: istwfk(nkpt),kg(3,mpw*mkmem),nband(nkpt*nsppol)
 integer,intent(in) :: ngfft(18),npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),kptns(3,nkpt)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),rprimd(3,3),wtk(nkpt)
 real(dp),intent(out) :: eltfrkin(6,6)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,iband,icg,ierr,ii,ikg
 integer :: ikpt,index,ipw,isppol,istwf_k,jj,master,me,n1,n2
 integer :: n3,nband_k,nkinout,npw_k,spaceComm
 real(dp) :: ucvol
!arrays
 integer,allocatable :: gbound(:,:),kg_k(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),rmet(3,3),tsec(2)
 real(dp),allocatable :: cwavef(:,:),ekinout(:)
 real(dp),allocatable :: eltfrkink(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!Default for sequential use
 master=0
!Init mpi_comm
 spaceComm=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt

!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 eltfrkin(:,:)=0.0_dp
 bdtot_index=0
 icg=0

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(cwavef,(2,mpw*nspinor))
 ABI_ALLOCATE(eltfrkink,(6,6))

!Define k-points distribution

!LOOP OVER SPINS
 do isppol=1,nsppol
   ikg=0

!  Loop over k points
   do ikpt=1,nkpt

     nband_k=nband(ikpt+(isppol-1)*nkpt)
     istwf_k=istwfk(ikpt)
     npw_k=npwarr(ikpt)

!    Skip this k-point if not the proper processor
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       bdtot_index=bdtot_index+nband_k
       cycle
     end if

     ABI_ALLOCATE(gbound,(2*mgfft+8,2))
     kpoint(:)=kptns(:,ikpt)

     kg_k(:,:) = 0

!$OMP PARALLEL DO PRIVATE(ipw) SHARED(ikg,kg,kg_k,npw_k)
     do ipw=1,npw_k
       kg_k(1,ipw)=kg(1,ipw+ikg)
       kg_k(2,ipw)=kg(2,ipw+ikg)
       kg_k(3,ipw)=kg(3,ipw+ikg)
     end do

     call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)

     index=1+icg

     eltfrkink(:,:)=0.0_dp

     nkinout=6*6
     ABI_ALLOCATE(ekinout,(nkinout))
     ekinout(:)=zero

     do iband=1,nband_k

       if(mpi_enreg%proc_distrb(ikpt,iband,isppol) /= me) cycle

       cwavef(:,1:npw_k*nspinor)=cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)

       call d2kindstr2(cwavef,ecut,ecutsm,effmass_free,ekinout,gmet,gprimd,&
&       istwf_k,kg_k,kpoint,npw_k,nspinor)

       eltfrkink(:,:)=eltfrkink(:,:)+ occ(iband+bdtot_index)* reshape(ekinout(:), (/6,6/) )

     end do !iband

     ABI_DEALLOCATE(ekinout)

     eltfrkin(:,:)=eltfrkin(:,:)+wtk(ikpt)*eltfrkink(:,:)

     ABI_DEALLOCATE(gbound)

     bdtot_index=bdtot_index+nband_k

     if (mkmem/=0) then
!      Handle case in which kg, cg, are kept in core
       icg=icg+npw_k*nspinor*nband_k
       ikg=ikg+npw_k
     end if

   end do
 end do  ! End loops on isppol and ikpt

!Fill in lower triangle
 do jj=2,6
   do ii=1,jj-1
     eltfrkin(jj,ii)=eltfrkin(ii,jj)
   end do
 end do

!Accumulate eltfrkin on all proc.
 call timab(48,1,tsec)
 call xmpi_sum(eltfrkin,spaceComm,ierr)
 call timab(48,2,tsec)

 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(eltfrkink)
 ABI_DEALLOCATE(kg_k)

 DBG_EXIT("COLL")

  contains
!!***

!!****f* ABINIT/d2kindstr2
!! NAME
!! d2kindstr2
!!
!! FUNCTION
!! compute expectation value of the second derivatives of the kinetic energy
!! wrt strain for one band and kpoint
!!
!! INPUTS
!!  cwavef(2,npw*nspinor)=wavefunction for current band
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!  effmass_free=effective mass for electrons (1. in common case)
!!  gmet(3,3)=reciprocal lattice metric tensor ($\textrm{Bohr}^{-2}$)
!!  gprimd(3,3)=primitive vectors in reciprocal space
!!  istwfk=information about wavefunction storage
!!  kg_k(3,npw)=integer coordinates of planewaves in basis sphere.
!!  kpt(3)=reduced coordinates of k point
!!  npw=number of plane waves at kpt.
!!  nspinor=number of spinorial components of the wavefunction
!!
!! OUTPUT
!!  ekinout(36)=expectation values of the second strain derivatives
!!   of the (modified) kinetic energy
!!
!! NOTES
!! Usually, the kinetic energy expression is $(1/2) (2 \pi)^2 (k+G)^2 $
!! However, the present implementation allows for a modification
!! of this kinetic energy, in order to obtain smooth total energy
!! curves with respect to the cut-off energy or the cell size and shape.
!! Thus the usual expression is kept if it is lower then ecut-ecutsm,
!! zero is returned beyond ecut, and in between, the kinetic
!! energy is DIVIDED by a smearing factor (to make it infinite at the
!! cut-off energy). The smearing factor is $x^2 (3-2x)$, where
!! x = (ecut- unmodified energy)/ecutsm.
!!
!! PARENTS
!!      dfpt_eltfrkin
!!
!! CHILDREN
!!
!! SOURCE

subroutine d2kindstr2(cwavef,ecut,ecutsm,effmass_free,ekinout,gmet,gprimd,&
&            istwfk,kg_k,kpt,npw,nspinor)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwfk,npw,nspinor
 real(dp),intent(in) :: ecut,ecutsm,effmass_free
!arrays
 integer,intent(in) :: kg_k(3,npw)
 real(dp),intent(in) :: cwavef(2,npw*nspinor),gmet(3,3),gprimd(3,3),kpt(3)
 real(dp),intent(inout) :: ekinout(36) !vz_i

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: ig,igs,ii,ispinor,istr1,istr2,ka,kb,kd,kg
 real(dp) :: d2fkin,d2fsm,d2kinacc,d2kpg2,dfkin,dfsm,dkpg21,dkpg22,ecutsm_inv
 real(dp) :: fsm,gpk1,gpk2,gpk3,htpisq,kpg2,term,xx
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 real(dp) :: d2gm(3,3),dgm01(3,3),dgm10(3,3)

! *************************************************************************
!
!htpisq is (1/2) (2 Pi) **2:
   htpisq=0.5_dp*(two_pi)**2

   ecutsm_inv=0.0_dp
   if(ecutsm>1.0d-20)ecutsm_inv=1/ecutsm

!Loop over 2nd strain index
   do istr2=1,6
!  Loop over 1st strain index, upper triangle only
     do istr1=1,istr2

       ka=idx(2*istr1-1);kb=idx(2*istr1);kg=idx(2*istr2-1);kd=idx(2*istr2)

       do ii = 1,3
         dgm01(:,ii)=-(gprimd(ka,:)*gprimd(kb,ii)+gprimd(kb,:)*gprimd(ka,ii))
         dgm10(:,ii)=-(gprimd(kg,:)*gprimd(kd,ii)+gprimd(kd,:)*gprimd(kg,ii))
       end do

       d2gm(:,:)=0._dp
       do ii = 1,3
         if(ka==kg) d2gm(:,ii)=d2gm(:,ii)&
&         +gprimd(kb,:)*gprimd(kd,ii)+gprimd(kd,:)*gprimd(kb,ii)
         if(ka==kd) d2gm(:,ii)=d2gm(:,ii)&
&         +gprimd(kb,:)*gprimd(kg,ii)+gprimd(kg,:)*gprimd(kb,ii)
         if(kb==kg) d2gm(:,ii)=d2gm(:,ii)&
&         +gprimd(ka,:)*gprimd(kd,ii)+gprimd(kd,:)*gprimd(ka,ii)
         if(kb==kd) d2gm(:,ii)=d2gm(:,ii)&
&         +gprimd(ka,:)*gprimd(kg,ii)+gprimd(kg,:)*gprimd(ka,ii)
       end do
       d2gm(:,:)=0.5_dp*d2gm(:,:)

       d2kinacc=0._dp

!    loop on spinor index
       do ispinor=1,nspinor
         igs=(ispinor-1)*npw
!      loop on plane waves
         do ig=1,npw
           gpk1=dble(kg_k(1,ig))+kpt(1)
           gpk2=dble(kg_k(2,ig))+kpt(2)
           gpk3=dble(kg_k(3,ig))+kpt(3)
           kpg2=htpisq*&
&           ( gmet(1,1)*gpk1**2+         &
&           gmet(2,2)*gpk2**2+         &
&           gmet(3,3)*gpk3**2          &
&           +2.0_dp*(gpk1*gmet(1,2)*gpk2+  &
&           gpk1*gmet(1,3)*gpk3+  &
&           gpk2*gmet(2,3)*gpk3 )  )
           dkpg21=htpisq*&
&           ( dgm01(1,1)*gpk1**2+         &
&           dgm01(2,2)*gpk2**2+         &
&           dgm01(3,3)*gpk3**2          &
&           +2.0_dp*(gpk1*dgm01(1,2)*gpk2+  &
&           gpk1*dgm01(1,3)*gpk3+  &
&           gpk2*dgm01(2,3)*gpk3 )  )
           dkpg22=htpisq*&
&           ( dgm10(1,1)*gpk1**2+         &
&           dgm10(2,2)*gpk2**2+         &
&           dgm10(3,3)*gpk3**2          &
&           +2.0_dp*(gpk1*dgm10(1,2)*gpk2+  &
&           gpk1*dgm10(1,3)*gpk3+  &
&           gpk2*dgm10(2,3)*gpk3 )  )
           d2kpg2=htpisq*&
&           ( d2gm(1,1)*gpk1**2+         &
&           d2gm(2,2)*gpk2**2+         &
&           d2gm(3,3)*gpk3**2          &
&           +2.0_dp*(gpk1*d2gm(1,2)*gpk2+  &
&           gpk1*d2gm(1,3)*gpk3+  &
&           gpk2*d2gm(2,3)*gpk3 )  )

           if(kpg2>ecut-tol12)then
             dfkin=0._dp
             d2fkin=0._dp
           elseif(kpg2>ecut-ecutsm)then
!          This kinetic cutoff smoothing function and its xx derivatives
!          were produced with Mathematica and the fortran code has been
!          numerically checked against Mathematica.
             xx=(ecut-kpg2)*ecutsm_inv
             fsm=1.0_dp/(xx**2*(3+xx*(1+xx*(-6+3*xx))))
             dfsm=-3.0_dp*(-1+xx)**2*xx*(2+5*xx)*fsm**2
             d2fsm=6.0_dp*xx**2*(9+xx*(8+xx*(-52+xx*(-3+xx*(137+xx*&
&             (-144+45*xx))))))*fsm**3
             dfkin=fsm-ecutsm_inv*kpg2*dfsm
             d2fkin=ecutsm_inv*(-2.0_dp*dfsm+ecutsm_inv*kpg2*d2fsm)
           else
             dfkin=1._dp
             d2fkin=0._dp
           end if

!        accumulate kinetic energy 2nd derivative with wavefunction components
           term=d2fkin*dkpg21*dkpg22 + dfkin*d2kpg2
           if(istwfk==2 .and. ig/=1)term=2.0_dp*term
           if(istwfk>2)term=2.0_dp*term
           d2kinacc=d2kinacc + term*(cwavef(re,ig+igs)**2 + cwavef(im,ig+igs)**2)

         end do  !ig
       end do !ispinor

       ekinout(istr1+6*(istr2-1))=d2kinacc/effmass_free

     end do !istr1
   end do !istr2

  end subroutine d2kindstr2
!!***

end subroutine dfpt_eltfrkin
!!***

!!****f* ABINIT/dfpt_eltfrhar
!! NAME
!! dfpt_eltfrhar
!!
!! FUNCTION
!! Compute the frozen-wavefunction hartree enegy contribution to the elastic tensor
!!
!! INPUTS
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  gsqcut =Fourier cutoff on G^2 for "large sphere" of radius double
!!   that of the basis sphere--appropriate for charge density rho(G),
!!   Hartree potential, and pseudopotentials
!!  mpi_enreg=informations about MPI parallelization
!!  nfft =(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/variables/vargs.htm#ngfft
!!  rhog(2,nfft)=total electron density in G space
!!
!! OUTPUT
!!  eltfrhar(6,6)=non-symmetrized kinetic energy contribution to the
!!                    elastic tensor
!! NOTES
!! *based largely on hartre.f
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      metric,ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine dfpt_eltfrhar(eltfrhar,rprimd,gsqcut,mpi_enreg,nfft,ngfft,rhog)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhog(2,nfft),rprimd(3,3)
 real(dp),intent(out) :: eltfrhar(6,6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i23,i3,id2,id3,ierr,ig,ig2,ig3,ii,ii1,ing,istr1,istr2,jj
 integer :: ka,kb,kd,kg,me_fft,n1,n2,n3,nproc_fft
 real(dp),parameter :: tolfix=1.000000001_dp
 real(dp) :: cutoff,d2eacc,d2etot,d2gs,deacc01,deacc10,dgs01,dgs10,eacc,fact,gs
 real(dp) :: term,ucvol
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer :: id(3)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: d2gm(3,3),dgm01(3,3),dgm10(3,3),gmet(3,3),gprimd(3,3),gqr(3)
 real(dp) :: rmet(3,3),tsec(2)
 real(dp),allocatable :: gq(:,:)

! *************************************************************************

!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 eltfrhar(:,:)=0.0_dp

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 me_fft=ngfft(11)
 nproc_fft=ngfft(10)

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Initialize a few quantities
 fact=0.5_dp*ucvol/pi
 cutoff=gsqcut*tolfix

!In order to speed the routine, precompute the components of g+q
!Also check if the booked space was large enough...
 ABI_ALLOCATE(gq,(3,max(n1,n2,n3)))
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id(ii))*ngfft(ii)-1
     gq(ii,ing)=ig
   end do
 end do

!Loop over 2nd strain index
 do istr2=1,6
!  Loop over 1st strain index, upper triangle only
   do istr1=1,istr2

     ka=idx(2*istr1-1);kb=idx(2*istr1);kg=idx(2*istr2-1);kd=idx(2*istr2)

     do ii = 1,3
       dgm01(:,ii)=-(gprimd(ka,:)*gprimd(kb,ii)+gprimd(kb,:)*gprimd(ka,ii))
       dgm10(:,ii)=-(gprimd(kg,:)*gprimd(kd,ii)+gprimd(kd,:)*gprimd(kg,ii))
     end do

     d2gm(:,:)=0._dp
     do ii = 1,3
       if(ka==kg) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(kb,:)*gprimd(kd,ii)+gprimd(kd,:)*gprimd(kb,ii)
       if(ka==kd) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(kb,:)*gprimd(kg,ii)+gprimd(kg,:)*gprimd(kb,ii)
       if(kb==kg) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(ka,:)*gprimd(kd,ii)+gprimd(kd,:)*gprimd(ka,ii)
       if(kb==kd) d2gm(:,ii)=d2gm(:,ii)&
&       +gprimd(ka,:)*gprimd(kg,ii)+gprimd(kg,:)*gprimd(ka,ii)
     end do
     d2gm(:,:)=0.5_dp*d2gm(:,:)

!    initialize energy accumulator
     eacc=0._dp
     deacc01=0._dp
     deacc10=0._dp
     d2eacc=0._dp

     id2=n2/2+2
     id3=n3/2+2
!    Triple loop on each dimension
     do i3=1,n3
       ig3=i3-(i3/id3)*n3-1
       gqr(3)=gq(3,i3)
       do i2=1,n2
         if (fftn2_distrib(i2)==me_fft) then
           gqr(2)=gq(2,i2)
           ig2=i2-(i2/id2)*n2-1
           i23=n1*(ffti2_local(i2)-1 +(n2/nproc_fft)*(i3-1))
!          Do the test that eliminates the Gamma point outside
!          of the inner loop
           ii1=1
           if(i23==0 .and. ig2==0 .and. ig3==0)then
             ii1=2
           end if

!          Final inner loop on the first dimension
!          (note the lower limit)
           do i1=ii1,n1
             gqr(1)=gq(1,i1)
             gs=(gmet(1,1)*gqr(1)*gqr(1)+gmet(2,2)*gqr(2)*gqr(2)+&
&             gmet(3,3)*gqr(3)*gqr(3)+2._dp*&
&             (gmet(1,2)*gqr(1)*gqr(2) + gmet(1,3)*gqr(1)*gqr(3)+&
&             gmet(2,3)*gqr(2)*gqr(3)) )
             ii=i1+i23
             if(gs<=cutoff)then
               dgs01=(dgm01(1,1)*gqr(1)*gqr(1)+dgm01(2,2)*gqr(2)*gqr(2)+&
&               dgm01(3,3)*gqr(3)*gqr(3)+2._dp*&
&               (dgm01(1,2)*gqr(1)*gqr(2) + dgm01(1,3)*gqr(1)*gqr(3)+&
&               dgm01(2,3)*gqr(2)*gqr(3)) )
               dgs10=(dgm10(1,1)*gqr(1)*gqr(1)+dgm10(2,2)*gqr(2)*gqr(2)+&
&               dgm10(3,3)*gqr(3)*gqr(3)+2._dp*&
&               (dgm10(1,2)*gqr(1)*gqr(2) + dgm10(1,3)*gqr(1)*gqr(3)+&
&               dgm10(2,3)*gqr(2)*gqr(3)) )
               d2gs =(d2gm(1,1)*gqr(1)*gqr(1)+d2gm(2,2)*gqr(2)*gqr(2)+&
&               d2gm(3,3)*gqr(3)*gqr(3)+2._dp*&
&               (d2gm(1,2)*gqr(1)*gqr(2) + d2gm(1,3)*gqr(1)*gqr(3)+&
&               d2gm(2,3)*gqr(2)*gqr(3)) )

               term=(rhog(re,ii)**2+rhog(im,ii)**2)/gs
               eacc=eacc+term
               deacc01=deacc01+dgs01*term/gs
               deacc10=deacc10+dgs10*term/gs
               d2eacc=d2eacc+(-d2gs+2._dp*dgs01*dgs10/gs)*term/gs
             end if

!            End loop on i1
           end do
         end if
!        End loop on i2
       end do
!      End loop on i3
     end do

!    Add contributions taking account diagonal strain terms (from ucvol
!    derivatives)
     d2etot=d2eacc
     if(istr1<=3) d2etot=d2etot+deacc10
     if(istr2<=3) d2etot=d2etot+deacc01
     if(istr1<=3 .and. istr2<=3) d2etot=d2etot+eacc

     eltfrhar(istr1,istr2)=fact*d2etot

!    End loop on istr1
   end do
!  End loop in istr2
 end do

 ABI_DEALLOCATE(gq)

!Init mpi_comm
 call timab(48,1,tsec)
 call xmpi_sum(eltfrhar,mpi_enreg%comm_fft,ierr)
 call timab(48,2,tsec)

!Fill in lower triangle
 do jj=2,6
   do ii=1,jj-1
     eltfrhar(jj,ii)=eltfrhar(ii,jj)
   end do
 end do
end subroutine dfpt_eltfrhar
!!***

!!****f* ABINIT/elt_ewald
!!
!! NAME
!! elt_ewald
!!
!! FUNCTION
!! Compute 2nd derivatives of Ewald energy wrt strain for frozen wavefunction
!! contributions to elastic tensor
!!
!! INPUTS
!! gmet(3,3)=metric tensor in reciprocal space (bohr^-2)
!! gprimd(3,3)=dimensional primitive translations for reciprocal space (bohr^-1)
!! mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!! comm_atom=--optional-- MPI communicator over atoms
!! my_natom=number of atoms treated by current processor
!! natom=number of atoms in unit cell
!! ntypat=numbe of type of atoms
!! rmet(3,3)=metric tensor in real space (bohr^2)
!! rprimd(3,3)=dimensional primitive translation vectors (bohr)
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume (bohr^3)
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!! zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!! elteew(6+3*natom,6)=2nd derivatives of Ewald energy wrt strain
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,timab,wrtout,xmpi_sum
!!
!! SOURCE

subroutine elt_ewald(elteew,gmet,gprimd,my_natom,natom,ntypat,rmet,rprimd,&
&                 typat,ucvol,xred,zion, &
&                 mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom,natom,ntypat
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: typat(natom)
 integer,optional,intent(in) :: comm_atom
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom),zion(ntypat)
 real(dp),intent(out) :: elteew(6+3*natom,6)

!Local variables-------------------------------
!scalars
 integer :: ia,ia0,ib,ierr,ig1,ig2,ig3,ir1,ir2,ir3,is1,is2,jj,js,ka,kb,kd,kg,my_comm_atom,newg,newr,ng,nr
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: arg,ch,chsq,cos_term,d2derfc,d2gss,d2r,d2rs,dderfc,derfc_arg
 real(dp) :: dgss1,dgss2,direct,dr1,dr2,drs1,drs2,eew,eta,fac,fraca1,fraca2
 real(dp) :: fraca3,fracb1,fracb2,fracb3,gsq,gsum,r1,r2,r3,recip,reta
 real(dp) :: rmagn,rsq,sin_term,sumg,summi,summr,sumr,t1,term
 character(len=500) :: message
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer,pointer :: my_atmtab(:)
 real(dp) :: d2gm(3,3,6,6),d2ris(3),d2rm(3,3,6,6),dgm(3,3,6),dris(3),drm(3,3,6)
 real(dp) :: t2(3),ts2(3),tsec(2),tt(3)
 real(dp),allocatable :: d2sumg(:,:),d2sumr(:,:),drhoisi(:,:),drhoisr(:,:)
 real(dp),allocatable :: mpibuf(:)

! *************************************************************************

!DEBUG
!write(std_out,*)' elt_ewald : enter '
!stop
!ENDDEBUG

!Compute 1st and 2nd derivatives of metric tensor wrt all strain components
!and store for use in inner loop below.

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Loop over 2nd strain index
 do is2=1,6
   kg=idx(2*is2-1);kd=idx(2*is2)
   do jj = 1,3
     drm(:,jj,is2)=rprimd(kg,:)*rprimd(kd,jj)+rprimd(kd,:)*rprimd(kg,jj)
     dgm(:,jj,is2)=-(gprimd(kg,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(kg,jj))
   end do

!  Loop over 1st strain index, upper triangle only
   do is1=1,is2

     ka=idx(2*is1-1);kb=idx(2*is1)
     d2rm(:,:,is1,is2)=zero
     d2gm(:,:,is1,is2)=zero
     do jj = 1,3
       if(ka==kg) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(kb,:)*rprimd(kd,jj)+rprimd(kd,:)*rprimd(kb,jj)
       if(ka==kd) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(kb,:)*rprimd(kg,jj)+rprimd(kg,:)*rprimd(kb,jj)
       if(kb==kg) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(ka,:)*rprimd(kd,jj)+rprimd(kd,:)*rprimd(ka,jj)
       if(kb==kd) d2rm(:,jj,is1,is2)=d2rm(:,jj,is1,is2)&
&       +rprimd(ka,:)*rprimd(kg,jj)+rprimd(kg,:)*rprimd(ka,jj)

       if(ka==kg) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(kb,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(kb,jj)
       if(ka==kd) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(kb,:)*gprimd(kg,jj)+gprimd(kg,:)*gprimd(kb,jj)
       if(kb==kg) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(ka,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(ka,jj)
       if(kb==kd) d2gm(:,jj,is1,is2)=d2gm(:,jj,is1,is2)&
&       +gprimd(ka,:)*gprimd(kg,jj)+gprimd(kg,:)*gprimd(ka,jj)
     end do
   end do !is1
 end do !is2

!Add up total charge and sum of $charge^2$ in cell
 chsq=zero
 ch=zero
 do ia=1,natom
   ch=ch+zion(typat(ia))
   chsq=chsq+zion(typat(ia))**2
 end do

!Compute eta, the Ewald summation convergence parameter,
!for approximately optimized summations:
 direct=rmet(1,1)+rmet(1,2)+rmet(1,3)+rmet(2,1)+&
& rmet(2,2)+rmet(2,3)+rmet(3,1)+rmet(3,2)+rmet(3,3)
 recip=gmet(1,1)+gmet(1,2)+gmet(1,3)+gmet(2,1)+&
& gmet(2,2)+gmet(2,3)+gmet(3,1)+gmet(3,2)+gmet(3,3)
!Here, a bias is introduced, because G-space summation scales
!better than r space summation ! Note : debugging is the most
!easier at fixed eta.
 eta=pi*200._dp/33.0_dp*sqrt(1.69_dp*recip/direct)

!Conduct reciprocal space summations
 fac=pi**2/eta ; gsum=zero
 ABI_ALLOCATE(d2sumg,(6+3*natom,6))
 ABI_ALLOCATE(drhoisr,(3,natom))
 ABI_ALLOCATE(drhoisi,(3,natom))
 d2sumg(:,:)=zero

!Sum over G space, done shell after shell until all
!contributions are too small.
 ng=0
 do
   ng=ng+1
   newg=0

   do ig3=-ng,ng
     do ig2=-ng,ng
       do ig1=-ng,ng

!        Exclude shells previously summed over
         if(abs(ig1)==ng .or. abs(ig2)==ng .or. abs(ig3)==ng&
&         .or. ng==1 ) then

!          gsq is G dot G = |G|^2
           gsq=gmet(1,1)*dble(ig1*ig1)+gmet(2,2)*dble(ig2*ig2)+&
&           gmet(3,3)*dble(ig3*ig3)+2._dp*(gmet(2,1)*dble(ig1*ig2)+&
&           gmet(3,1)*dble(ig1*ig3)+gmet(3,2)*dble(ig3*ig2))

!          Skip g=0:
           if (gsq>1.0d-20) then
             arg=fac*gsq

!            Larger arg gives 0 contribution because of exp(-arg)
             if (arg <= 80._dp) then
!              When any term contributes then include next shell
               newg=1
               term=exp(-arg)/gsq
               summr = zero
               summi = zero
!              Note that if reduced atomic coordinates xred drift outside
!              of unit cell (outside [0,1)) it is irrelevant in the following
!              term, which only computes a phase.
               do ia=1,natom
                 arg=two_pi*(ig1*xred(1,ia)+ig2*xred(2,ia)+ig3*xred(3,ia))
!                Sum real and imaginary parts (avoid complex variables)
                 cos_term=cos(arg)
                 sin_term=sin(arg)
                 summr=summr+zion(typat(ia))*cos_term
                 summi=summi+zion(typat(ia))*sin_term
                 drhoisr(1,ia)=-two_pi*zion(typat(ia))*sin_term*dble(ig1)
                 drhoisi(1,ia)= two_pi*zion(typat(ia))*cos_term*dble(ig1)
                 drhoisr(2,ia)=-two_pi*zion(typat(ia))*sin_term*dble(ig2)
                 drhoisi(2,ia)= two_pi*zion(typat(ia))*cos_term*dble(ig2)
                 drhoisr(3,ia)=-two_pi*zion(typat(ia))*sin_term*dble(ig3)
                 drhoisi(3,ia)= two_pi*zion(typat(ia))*cos_term*dble(ig3)
               end do

!              The following two checks avoid an annoying
!              underflow error message
               if (abs(summr)<1.d-16) summr=zero
               if (abs(summi)<1.d-16) summi=zero

!              The product of term and summr**2 or summi**2 below
!              can underflow if not for checks above
               t1=term*(summr*summr+summi*summi)
               gsum=gsum+t1
!              Loop over 2nd strain index
               do is2=1,6
                 dgss2=dgm(1,1,is2)*dble(ig1*ig1)+dgm(2,2,is2)*dble(ig2*ig2)+&
&                 dgm(3,3,is2)*dble(ig3*ig3)+2._dp*(dgm(2,1,is2)*dble(ig1*ig2)+&
&                 dgm(3,1,is2)*dble(ig1*ig3)+dgm(3,2,is2)*dble(ig3*ig2))
!                Loop over 1st strain index, upper triangle only
                 do is1=1,is2
                   dgss1=dgm(1,1,is1)*dble(ig1*ig1)+dgm(2,2,is1)*dble(ig2*ig2)+&
&                   dgm(3,3,is1)*dble(ig3*ig3)+2._dp*(dgm(2,1,is1)*dble(ig1*ig2)+&
&                   dgm(3,1,is1)*dble(ig1*ig3)+dgm(3,2,is1)*dble(ig3*ig2))

                   d2gss=d2gm(1,1,is1,is2)*dble(ig1*ig1)+&
&                   d2gm(2,2,is1,is2)*dble(ig2*ig2)+&
&                   d2gm(3,3,is1,is2)*dble(ig3*ig3)+2._dp*&
&                   (d2gm(2,1,is1,is2)*dble(ig1*ig2)+&
&                   d2gm(3,1,is1,is2)*dble(ig1*ig3)+&
&                   d2gm(3,2,is1,is2)*dble(ig3*ig2))

                   d2sumg(is1,is2)=d2sumg(is1,is2)+&
&                   t1*((fac**2 + 2.0_dp*fac/gsq + 2.0_dp/(gsq**2))*dgss1*dgss2 -&
&                   0.5_dp*(fac + 1.0_dp/gsq)*d2gss)
                   if(is1<=3) d2sumg(is1,is2)=d2sumg(is1,is2)+&
&                   t1*(fac + 1.0_dp/gsq)*dgss2
                   if(is2<=3) d2sumg(is1,is2)=d2sumg(is1,is2)+&
&                   t1*(fac + 1.0_dp/gsq)*dgss1
                   if(is1<=3 .and. is2<=3) d2sumg(is1,is2)=d2sumg(is1,is2)+t1

                 end do !is1

!                Internal strain contributions
                 do ia=1,natom
                   js=7+3*(ia-1)
                   t2(:)=2.0_dp*term*(summr*drhoisr(:,ia)+summi*drhoisi(:,ia))
                   d2sumg(js:js+2,is2)=d2sumg(js:js+2,is2)-&
&                   (fac + 1.0_dp/gsq)*dgss2*t2(:)
                   if(is2<=3) d2sumg(js:js+2,is2)=d2sumg(js:js+2,is2)-t2(:)
                 end do
               end do !is2

             end if ! End condition of not larger than 80.0
           end if ! End skip g=0
         end if ! End triple loop over G s and associated new shell condition
       end do
     end do
   end do

!  Check if new shell must be calculated
   if (newg==0) exit
 end do !  End the loop on ng (new shells). Note that there is one exit from this loop.

 sumg=gsum/(two_pi*ucvol)
 d2sumg(:,:)=d2sumg(:,:)/(two_pi*ucvol)

 ABI_DEALLOCATE(drhoisr)
 ABI_DEALLOCATE(drhoisi)
!Stress tensor is now computed elsewhere (ewald2) hence do not need
!length scale gradients (used to compute them here).

!Conduct real space summations
 reta=sqrt(eta)
 fac=2._dp*sqrt(eta/pi)
 ABI_ALLOCATE(d2sumr,(6+3*natom,6))
 sumr=zero;d2sumr(:,:)=zero

!In the following a summation is being conducted over all
!unit cells (ir1, ir2, ir3) so it is appropriate to map all
!reduced coordinates xred back into [0,1).
!
!Loop on shells in r-space as was done in g-space
 nr=0
 do
   nr=nr+1
   newr=0
!
   do ir3=-nr,nr
     do ir2=-nr,nr
       do ir1=-nr,nr
         if( abs(ir3)==nr .or. abs(ir2)==nr .or. abs(ir1)==nr .or. nr==1 )then

           do ia0=1,my_natom
             ia=ia0;if(paral_atom)ia=my_atmtab(ia0)
             js=7+3*(ia-1)
!            Map reduced coordinate xred(mu,ia) into [0,1)
             fraca1=xred(1,ia)-aint(xred(1,ia))+0.5_dp-sign(0.5_dp,xred(1,ia))
             fraca2=xred(2,ia)-aint(xred(2,ia))+0.5_dp-sign(0.5_dp,xred(2,ia))
             fraca3=xred(3,ia)-aint(xred(3,ia))+0.5_dp-sign(0.5_dp,xred(3,ia))
             do ib=1,natom
               fracb1=xred(1,ib)-aint(xred(1,ib))+0.5_dp-sign(0.5_dp,xred(1,ib))
               fracb2=xred(2,ib)-aint(xred(2,ib))+0.5_dp-sign(0.5_dp,xred(2,ib))
               fracb3=xred(3,ib)-aint(xred(3,ib))+0.5_dp-sign(0.5_dp,xred(3,ib))
               r1=dble(ir1)+fracb1-fraca1
               r2=dble(ir2)+fracb2-fraca2
               r3=dble(ir3)+fracb3-fraca3
               rsq=rmet(1,1)*r1*r1+rmet(2,2)*r2*r2+rmet(3,3)*r3*r3+&
&               2.0_dp*(rmet(2,1)*r2*r1+rmet(3,2)*r3*r2+rmet(3,1)*r1*r3)

!              Avoid zero denominators in 'term':
               if (rsq>=1.0d-24) then

!                Note: erfc(8) is about 1.1e-29,
!                so do not bother with larger arg.
!                Also: exp(-64) is about 1.6e-28,
!                so do not bother with larger arg**2 in exp.
                 term=zero
                 if (eta*rsq<64.0_dp) then
                   newr=1
                   rmagn=sqrt(rsq)
                   arg=reta*rmagn
!                  derfc computes the complementary error function
!                  dderfc is the derivative of the complementary error function
!                  d2derfc is the 2nd derivative of the complementary error function
                   dderfc=-fac*exp(-eta*rsq)
                   d2derfc=-2._dp*eta*rmagn*dderfc
                   derfc_arg = abi_derfc(arg)
                   term=derfc_arg/rmagn
                   sumr=sumr+zion(typat(ia))*zion(typat(ib))*term
                   tt(:)=rmet(:,1)*r1+rmet(:,2)*r2+rmet(:,3)*r3
                   dris(:)=tt(:)/rmagn
!                  Loop over 2nd strain index
                   do is2=1,6
                     drs2=drm(1,1,is2)*r1*r1+drm(2,2,is2)*r2*r2+&
&                     drm(3,3,is2)*r3*r3+&
&                     2.0_dp*(drm(2,1,is2)*r2*r1+drm(3,2,is2)*r3*r2+&
&                     drm(3,1,is2)*r1*r3)
                     dr2=0.5_dp*drs2/rmagn
!                    Loop over 1st strain index, upper triangle only
                     do is1=1,is2
                       drs1=drm(1,1,is1)*r1*r1+drm(2,2,is1)*r2*r2+&
&                       drm(3,3,is1)*r3*r3+&
&                       2.0_dp*(drm(2,1,is1)*r2*r1+drm(3,2,is1)*r3*r2+&
&                       drm(3,1,is1)*r1*r3)
                       dr1=0.5_dp*drs1/rmagn
                       d2rs=d2rm(1,1,is1,is2)*r1*r1+d2rm(2,2,is1,is2)*r2*r2+&
&                       d2rm(3,3,is1,is2)*r3*r3+&
&                       2.0_dp*(d2rm(2,1,is1,is2)*r2*r1+d2rm(3,2,is1,is2)*r3*r2+&
&                       d2rm(3,1,is1,is2)*r1*r3)
                       d2r=(0.25_dp*d2rs-dr1*dr2)/rmagn
                       d2sumr(is1,is2)=d2sumr(is1,is2)+&
&                       zion(typat(ia))*zion(typat(ib))*&
&                       ((d2derfc-2.0_dp*dderfc/rmagn+2.0_dp*derfc_arg/rsq)*dr1*dr2+&
&                       (dderfc-derfc_arg/rmagn)*d2r)/rmagn
                     end do !is1
!                    Internal strain contribution
                     ts2(:)=drm(:,1,is2)*r1+drm(:,2,is2)*r2+drm(:,3,is2)*r3
                     d2ris(:)=ts2(:)/rmagn-0.5_dp*drs2*tt(:)/(rsq*rmagn)

                     d2sumr(js:js+2,is2)=d2sumr(js:js+2,is2)-&
&                     2.0_dp*zion(typat(ia))*zion(typat(ib))*&
&                     ((d2derfc-2.0_dp*dderfc/rmagn+2.0_dp*derfc_arg/rsq)*dr2*dris(:)+&
&                     (dderfc-derfc_arg/rmagn)*d2ris(:))/rmagn
                   end do !is2
                 end if
               end if ! End avoid zero denominators in'term'

             end do ! end loop over ib:
           end do ! end loop over ia:
         end if ! end triple loop over real space points and associated condition of new shell
       end do
     end do
   end do

!  Check if new shell must be calculated
   if(newr==0) exit
 end do !  End loop on nr (new shells). Note that there is an exit within the loop

!In case of parallelism over atoms: communicate
 if (paral_atom) then
   call timab(48,1,tsec)
   ABI_ALLOCATE(mpibuf,((6+3*natom)*6+1))
   mpibuf(1:(6+3*natom)*6)=reshape(d2sumr(:,:),shape=(/((6+3*natom)*6)/))
   mpibuf((6+3*natom)*6+1)=sumr
   call xmpi_sum(mpibuf,my_comm_atom,ierr)
   sumr=mpibuf((6+3*natom)*6+1)
   d2sumr(:,:)=reshape(mpibuf(1:(6+3*natom)*6),shape=(/(6+3*natom),6/))
   ABI_DEALLOCATE(mpibuf)
   call timab(48,2,tsec)
 end if

 sumr=0.5_dp*sumr
 d2sumr(:,:)=0.5_dp*d2sumr(:,:)
 fac=pi*ch**2/(2.0_dp*eta*ucvol)

!Finally assemble Ewald energy, eew
 eew=sumg+sumr-chsq*reta/sqrt(pi)-fac

 elteew(:,:)=d2sumg(:,:)+d2sumr(:,:)

!Additional term for all strains diagonal (from "fac" term in eew)
 elteew(1:3,1:3)=elteew(1:3,1:3)-fac

!Fill in lower triangle
 do is2=2,6
   do is1=1,is2-1
     elteew(is2,is1)=elteew(is1,is2)
   end do
 end do

 ABI_DEALLOCATE(d2sumg)
 ABI_DEALLOCATE(d2sumr)

!Output the final values of ng and nr
 write(message, '(a,i4,a,i4)' )' elt_ewald : nr and ng are ',nr,' and ',ng
 call wrtout(std_out,message,'COLL')

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine elt_ewald
!!***

!!****f* ABINIT/dfpt_ewald
!!
!! NAME
!! dfpt_ewald
!!
!! FUNCTION
!! Compute ewald contribution to the dynamical matrix, at a given q wavevector.
!! Note: the q=0 part should be subtracted, by another call to
!! the present routine, with q=0. The present routine correspond
!! to the quantity C_bar defined in Eq.(24) or (27) in Phys. Rev. B 55, 10355 (1997) [[cite:Gonze1997a]].
!! The two calls correspond to Eq.(23) of the same paper.
!! If q=0 is asked, sumg0 should be put to 0. Otherwise, it should be put to 1.
!!
!! INPUTS
!! gmet(3,3)=metric tensor in reciprocal space (length units **-2)
!! mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!! comm_atom=--optional-- MPI communicator over atoms
!! my_natom=number of atoms treated by current processor
!! natom=number of atoms in unit cell
!! qphon(3)=phonon wavevector (same system of coordinates as the
!!          reciprocal lattice vectors)
!! rmet(3,3)=metric tensor in real space (length units squared)
!! sumg0: if=1, the sum in reciprocal space must include g=0,
!!   if=0, this contribution must be skipped (q=0 singularity)
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume in (whatever length scale units)**3
!! xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!! zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!! dyew(2,3,natom,3,natom)= Ewald part of the dynamical matrix,
!!    second energy derivative wrt xred(3,natom), Hartrees.
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,timab,xmpi_sum
!!
!! SOURCE

subroutine dfpt_ewald(dyew,gmet,my_natom,natom,qphon,rmet,sumg0,typat,ucvol,xred,zion, &
&                 mpi_atmtab,comm_atom ) ! optional arguments (parallelism))

!Arguments -------------------------------
!scalars
 integer,intent(in) :: my_natom,natom,sumg0
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: typat(natom)
 integer,optional,intent(in) :: comm_atom
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gmet(3,3),qphon(3),rmet(3,3),xred(3,natom),zion(*)
 real(dp),intent(out) :: dyew(2,3,natom,3,natom)

!Local variables -------------------------
!nr, ng affect convergence of sums (nr=3,ng=5 is not good enough):
!scalars
 integer,parameter :: im=2,ng=10,nr=6,re=1
 integer :: ia,ia0,ib,ierr,ig1,ig2,ig3,ii,ir1,ir2,ir3,mu,my_comm_atom,nu
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: arg,arga,argb,c1i,c1r,da1,da2,da3,derfc_arg
 real(dp) :: direct,dot1,dot2,dot3,dotr1,dotr2,dotr3
 real(dp) :: eta,fac,gdot12,gdot13,gdot23,gsq,gsum,norm1
 real(dp) :: r1,r2,r3,rdot12,rdot13,rdot23,recip,reta
 real(dp) :: reta3m,rmagn,rsq,term,term1,term2
 real(dp) :: term3
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 integer,pointer :: my_atmtab(:)
 real(dp) :: gpq(3),rq(3)

! *************************************************************************

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

!Compute eta for approximately optimized summations:
 direct=rmet(1,1)+rmet(1,2)+rmet(1,3)+rmet(2,1)+&
& rmet(2,2)+rmet(2,3)+rmet(3,1)+rmet(3,2)+rmet(3,3)
 recip=gmet(1,1)+gmet(1,2)+gmet(1,3)+gmet(2,1)+&
& gmet(2,2)+gmet(2,3)+gmet(3,1)+gmet(3,2)+gmet(3,3)
 eta=pi*(dble(ng)/dble(nr))*sqrt(1.69_dp*recip/direct)

!Test Ewald s summation
!eta=1.2_dp*eta

!Sum terms over g space:
 fac=pi**2/eta
 gsum=zero
 da1=zero
 da2=zero
 da3=zero
 dyew(:,:,:,:,:)=zero
 do ig3=-ng,ng
   do ig2=-ng,ng
     do ig1=-ng,ng
       gpq(1)=dble(ig1)+qphon(1)
       gpq(2)=dble(ig2)+qphon(2)
       gpq(3)=dble(ig3)+qphon(3)
       gdot12=gmet(2,1)*gpq(1)*gpq(2)
       gdot13=gmet(3,1)*gpq(1)*gpq(3)
       gdot23=gmet(3,2)*gpq(2)*gpq(3)
       dot1=gmet(1,1)*gpq(1)**2+gdot12+gdot13
       dot2=gmet(2,2)*gpq(2)**2+gdot12+gdot23
       dot3=gmet(3,3)*gpq(3)**2+gdot13+gdot23
       gsq=dot1+dot2+dot3
!      Skip q=0:
       if (gsq<1.0d-20) then
         if (sumg0==1) then
           write(message,'(5a)')&
&           'The phonon wavelength should not be zero : ',ch10,&
&           'there are non-analytical terms that the code cannot handle.',ch10,&
&           'Action : subtract this wavelength from the input.'
           MSG_ERROR(message)
         end if
       else
         arg=fac*gsq
!        Larger arg gives 0 contribution:
         if (arg <= 80._dp) then
           term=exp(-arg)/gsq
           do ia0=1,my_natom
             ia=ia0;if(paral_atom)ia=my_atmtab(ia0)
             arga=two_pi*(gpq(1)*xred(1,ia)+gpq(2)*xred(2,ia)+gpq(3)*xred(3,ia))
             do ib=1,ia
               argb=two_pi*(gpq(1)*xred(1,ib)+gpq(2)*xred(2,ib)+gpq(3)*xred(3,ib))
               arg=arga-argb
               c1r=cos(arg)*term
               c1i=sin(arg)*term

               do mu=1,3
                 do nu=1,mu
                   dyew(re,mu,ia,nu,ib)=dyew(re,mu,ia,nu,ib)+gpq(mu)*gpq(nu)*c1r
                   dyew(im,mu,ia,nu,ib)=dyew(im,mu,ia,nu,ib)+gpq(mu)*gpq(nu)*c1i
                 end do
               end do

             end do
           end do
         end if
!        Endif g/=0 :
       end if
!      End triple loop over G s:
     end do
   end do
 end do

!End G summation by accounting for some common factors.
!(for the charges:see end of routine)
 norm1=4.0_dp*pi/ucvol
 do ia0=1,my_natom
   ia=ia0;if(paral_atom)ia=my_atmtab(ia0)
   do ib=1,ia
     do mu=1,3
       do nu=1,mu
         dyew(:,mu,ia,nu,ib)=dyew(:,mu,ia,nu,ib)*norm1
       end do
     end do
   end do
 end do


!Do sums over real space:
 reta=sqrt(eta)
 reta3m=-eta*reta
 fac=4._dp/3.0_dp/sqrt(pi)
 do ir3=-nr,nr
   do ir2=-nr,nr
     do ir1=-nr,nr
       arg=-two_pi*(qphon(1)*ir1+qphon(2)*ir2+qphon(3)*ir3)
       c1r=cos(arg)*reta3m
       c1i=sin(arg)*reta3m
       do ia0=1,my_natom
         ia=ia0;if(paral_atom)ia=my_atmtab(ia0)
         do ib=1,ia
           r1=dble(ir1)+xred(1,ia)-xred(1,ib)
           r2=dble(ir2)+xred(2,ia)-xred(2,ib)
           r3=dble(ir3)+xred(3,ia)-xred(3,ib)
           rdot12=rmet(2,1)*r1*r2
           rdot13=rmet(3,1)*r1*r3
           rdot23=rmet(3,2)*r2*r3
           dotr1=rmet(1,1)*r1**2+rdot12+rdot13
           dotr2=rmet(2,2)*r2**2+rdot12+rdot23
           dotr3=rmet(3,3)*r3**2+rdot13+rdot23
           rsq=dotr1+dotr2+dotr3
           rmagn=sqrt(rsq)
!          Avoid zero denominators in term :
           if (rmagn>=1.0d-12) then
             arg=reta*rmagn
             term=zero
             if (arg<8.0_dp) then
!              Note: erfc(8) is about 1.1e-29,
!              so don t bother with larger arg.
!              Also: exp(-64) is about 1.6e-28,
!              so don t bother with larger arg**2 in exp.
               derfc_arg = abi_derfc(arg)
               term=derfc_arg/arg**3
               term1=2.0_dp/sqrt(pi)*exp(-arg**2)/arg**2
               term2=-(term+term1)
               term3=(3*term+term1*(3.0_dp+2.0_dp*arg**2))/rsq
               rq(1)=rmet(1,1)*r1+rmet(1,2)*r2+rmet(1,3)*r3
               rq(2)=rmet(2,1)*r1+rmet(2,2)*r2+rmet(2,3)*r3
               rq(3)=rmet(3,1)*r1+rmet(3,2)*r2+rmet(3,3)*r3
               do mu=1,3
                 do nu=1,mu
                   dyew(re,mu,ia,nu,ib)=dyew(re,mu,ia,nu,ib)+&
&                   c1r*(rq(mu)*rq(nu)*term3+rmet(mu,nu)*term2)
                   dyew(im,mu,ia,nu,ib)=dyew(im,mu,ia,nu,ib)+&
&                   c1i*(rq(mu)*rq(nu)*term3+rmet(mu,nu)*term2)
                 end do
               end do
             end if
           else
             if (ia/=ib)then
               write(message,'(a,a,a,a,a,i5,a,i5,a)')&
&               'The distance between two atoms vanishes.',ch10,&
&               'This is not allowed.',ch10,&
&               'Action: check the input for the atoms number',ia,' and',ib,'.'
               MSG_ERROR(message)
             else
               do mu=1,3
                 do nu=1,mu
                   dyew(re,mu,ia,nu,ib)=dyew(re,mu,ia,nu,ib)+&
&                   fac*reta3m*rmet(mu,nu)
                 end do
               end do
             end if
           end if

         end do ! End loop over ib:
       end do ! End loop over ia:
     end do ! End triple loop over real space points:
   end do
 end do

!Take account of the charges
!write(std_out,*)' '
 do ia0=1,my_natom
   ia=ia0;if(paral_atom)ia=my_atmtab(ia0)
   do ib=1,ia
     do mu=1,3
       do nu=1,mu
         do ii=1,2
!          write(std_out,*)dyew(ii,mu,ia,nu,ib)
           dyew(ii,mu,ia,nu,ib)=dyew(ii,mu,ia,nu,ib)*&
&           zion(typat(ia))*zion(typat(ib))
         end do
       end do
     end do
   end do
 end do

!Symmetrize with respect to the directions
 do ia0=1,my_natom
   ia=ia0;if(paral_atom)ia=my_atmtab(ia0)
   do ib=1,ia
     do mu=1,3
       do nu=1,mu
         dyew(re,nu,ia,mu,ib)=dyew(re,mu,ia,nu,ib)
         dyew(im,nu,ia,mu,ib)=dyew(im,mu,ia,nu,ib)
       end do
     end do
   end do
 end do

!In case of parallelism over atoms: communicate
 if (paral_atom) then
   call timab(48,1,tsec)
   call xmpi_sum(dyew,my_comm_atom,ierr)
   call timab(48,2,tsec)
 end if

!Fill the upper part of the matrix, with the hermitian conjugate
 do ia=1,natom
   do ib=1,ia
     do nu=1,3
       do mu=1,3
         dyew(re,mu,ib,nu,ia)=dyew(re,mu,ia,nu,ib)
         dyew(im,mu,ib,nu,ia)=-dyew(im,mu,ia,nu,ib)
       end do
     end do
   end do
 end do

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine dfpt_ewald
!!***

end module m_dfpt_elt
!!***
