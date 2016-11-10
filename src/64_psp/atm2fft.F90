!{\src2tex{textfont=tt}}
!!****f* ABINIT/atm2fft
!! NAME
!! atm2fft
!!
!! FUNCTION
!! This routine sums atomic functions (density or potential) defined
!! (in rec. space) on a radial grid to get global quantities on the
!! fine FFT grid. It can also compute contribution to energy derivatives
!! of these atomic functions.
!!
!! Possible options:
!!   optn=1: compute a sum of local atomic densities or contrib. to energy derivatives
!!   optv=1: compute a sum of local atomic potentials or contrib. to energy derivatives
!!
!!   optatm =1: computes sum of atomic potentials/densities
!!   optgr  =1: computes contribution of atomic pot./dens. to forces
!!   optstr =1: computes contribution of atomic pot./dens. to stress tensor
!!   optdyfr=1: computes contribution of atomic pot./dens. to frozen part of dyn. matrix
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  distribfft<type(distribfft_type)>=--optional-- contains infos related to FFT parallelism
!!  gauss(2,ntypat)= params for gaussian atm density (optn2=3) for each atom type
!!  gmet(3,3)=reciprocal space metric
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  gsqcut=cutoff on |G|^2: see setup1 for definition (doubled sphere)
!!  kxc(2,nfft)=exchange and correlation kernel
!!  me_g0=--optional-- 1 if the current process treat the g=0 plane-wave (only needed when comm_fft is present)
!!  mgfft=maximum size of 1D FFTs
!!  comm_fft=--optional-- MPI communicator over FFT components
!!  mqgrid=number of grid pts in q array for f(q) spline.
!!  natom=number of atoms in unit cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT
!!  ntypat=number of types of atoms.
!!  optatm,optdyfr,optgr,optn,optn2,optstr,optv= (see NOTES below)
!!  paral_kgb=--optional-- 1 if "band-FFT" parallelism is activated (only needed when comm_fft is present)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  qgrid(mqgrid)=q grid for spline from 0 to qmax.
!!  qprtrb(3)= integer wavevector of possible perturbing potential
!!             in basis of reciprocal lattice translations
!!  rhog(2,nfft)=electron density rho(G) in reciprocal space
!!               (used only if optv=1 and (optgr=1 or optstr=1 or optdyfr=1))
!!  ucvol=unit cell volume
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vspl(mqgrid,2,ntypat)=q^2 v(q) spline of an atomic potential
!!                        (used only if optv=1)
!!  vprtrb(2)=complex amplitude of possible perturbing potential; if nonzero,
!!            perturbing potential is added of the form V(G)=(vprtrb(1)+I*vprtrb(2))/2
!!            at the values G=qprtrb and (vprtrb(1)-I*vprtrb(2))/2 at G=-qprtrb
!!  vg(2,nfft)= potential V(G) in reciprocal space
!!              (used only if optn=1 and (optgr=1 or optstr=1 or optdyfr=1))
!!  vg1(2,nfft)= 1st-order potential V(G) in reciprocal space
!!               (used only if opteltfr==1)
!!  vg1_core(2,nfft)= 1-order potential V(G) in reciprocal space with only core contribution
!!                    (used only if opteltfr==1)
!! OUTPUT
!!  ======= if optv==1 =======
!!  ============================
!!   --- if optatm==1
!!    atmvloc(nfft)=sum of local atomic potentials in real space
!!   --- if optgr==1
!!    grv(3,natom)=contribution of atomic potentials to forces
!!   --- if optstr==1
!!    strv(6)=contribution of atomic potentials to stress tensor
!!            cart. coordinates, symmetric tensor, 6 comp. in order 11,22,33,32,31,21
!!   --- if optdyfr==1
!!    dyfrv(3,3,natom)=contribution of atomic potentials to frozen part of dyn. matrix
!!
!!  ======= if optn==1 =======
!!  ============================
!!   --- if optatm==1
!!    atmrho(nfft)=sum of atomic densities in real space
!!   --- if optgr==1
!!    grn(3,natom)=contribution of atomic densities to forces
!!   --- if optstr==1
!!    strn(6)=contribution of atomic densities to stress tensor
!!            cart. coordinates, symmetric tensor, 6 comp. in order 11,22,33,32,31,21
!!   --- if optdyfr==1
!!    dyfrn(3,3,natom)=contribution of atomic densities to frozen part of dyn. matrix
!!   --- if opteltfr==1
!!    eltfrn(6+3*natom,6)=contribution of atomic density to frozen part of stress tensor
!!
!! NOTES
!! Details on possible options:
!! ============================
!! optv: controls the computation of a local potential as sum of atomic potentials
!!          Vloc(r)=Sum_R[V^AT(r-R)]
!!       or its contributions to energy derivatives, i.e. derivatives of Int[Vloc(r).rho(r).dr]
!!          rho(r) is stored in reciprocal space in array rhog()
!!          V^AT is stored in reciprocal space in array vspl (in practice vspl(q)=q^2.V^AT(q))
!!
!! optn: controls the computation of a density as sum of atomic densities
!!          n(r)=Sum_R[n^AT(r-R)]
!!       or its contributions to energy derivatives, i.e. derivatives of Int[n(r).V(r).dr]
!!          V(r) is stored in reciprocal space in array vg()
!!          n^AT is stored in reciprocal space:
!!          if optn2=1: n^AT is the atomic PAW PS core density stored in array pawtab%tcorespl()
!!                   2: n^AT is the atomic PAW PS valence density stored in array pawtab%tvalespl()
!!                   3: n^AT is a gaussian density: n(g)=gauss(1,ityp)*exp[-(gauss(2,ityp)*G)^2]
!! Note: optv and optn can be activated together
!!
!! Options controlling which contrib. to Etot derivatives are computed:
!!   optatm  =1: computes Vloc(r) or n(r) as sum of atomic potentials/densities
!!   optgr   =1: computes contribution of atomic Vloc(r) or n(r) to forces
!!   optstr  =1: computes contribution of atomic Vloc(r) or n(r) to stress tensor
!!   optdyfr =1: computes contribution of atomic Vloc(r) or n(r) to fr part of dyn. matrix
!!   opteltfr=1: computes contribution of atomic Vloc(r) or n(r) to elastic tensor
!! Note: optatm, optgr, optstr, optelfr and optdyfr can be activated together
!!
!! Typical uses:
!! =============
!! Computation of:
!!  - local potential: optv=1, optatm=1
!!  - contrib. of local potential to Etot derivatives: optv=1, rhog=total valence density
!!                                                     optgr=1 or optstr=1 or optdyfr=1
!!  - PS core density: optn=1, optn2=1, optatm=1
!!  - contrib. of NLCC to Etot derivatives: optn=1, optn2=1, vg=XC potential
!!                                          optgr=1 or optstr=1 or optdyfr=1
!!  - sum of atomic valence densities: optn=1, optn2=2 or 3, optatm=1
!!  - correction of forces due to potential residual: optn=1, optn2=2 or 3, optgr=1
!!                                                    vg=potential residual
!!    etc...
!!
!! PARENTS
!!      dfpt_dyfro,dfpt_eltfrxc,extraprho,forces,fresidrsp,prcref,prcref_PMA
!!      respfn,setvtr,stress
!!
!! CHILDREN
!!      destroy_distribfft,fourdp,init_distribfft_seq,initmpi_seq
!!      set_mpi_enreg_fft,timab,unset_mpi_enreg_fft,wrtout,xmpi_sum,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine atm2fft(atindx1,atmrho,atmvloc,dyfrn,dyfrv,eltfrn,gauss,gmet,gprimd,&
&                  grn,grv,gsqcut,mgfft,mqgrid,natom,nattyp,nfft,ngfft,ntypat,&
&                  optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,&
&                  psps,pawtab,ph1d,qgrid,qprtrb,rhog,strn,strv,ucvol,usepaw,vg,vg1,vg1_core,vprtrb,vspl,&
&                  is2_in,comm_fft,me_g0,paral_kgb,distribfft) ! optional arguments

 use m_profiling_abi
 use m_errors
 use defs_basis
 use defs_abitypes
 use m_errors
 use m_xmpi

 use defs_datatypes,only : pseudopotential_type
 use m_pawtab,      only : pawtab_type
 use m_distribfft,  only : distribfft_type
 use m_mpinfo,      only : set_mpi_enreg_fft, unset_mpi_enreg_fft

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'atm2fft'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,mqgrid,natom,nfft,ntypat,optatm,optdyfr,opteltfr
 integer,intent(in) :: optgr,optn,optn2,optstr,optv,usepaw
 integer,optional,intent(in) :: is2_in,me_g0,comm_fft,paral_kgb
 real(dp),intent(in) :: gsqcut,ucvol
 type(pseudopotential_type),target,intent(in) :: psps
 type(distribfft_type),optional,intent(in),target :: distribfft
!arrays
 integer,intent(in) :: atindx1(natom),nattyp(ntypat),ngfft(18),qprtrb(3)
 real(dp),intent(in) :: gauss(2,ntypat*(optn2/3)),gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qgrid(mqgrid)
 real(dp),intent(in) :: rhog(2,nfft*optv*max(optgr,optstr,optdyfr,opteltfr))
 real(dp),intent(in) :: vg(2,nfft*optn*max(optgr,optstr,optdyfr,opteltfr))
 real(dp),intent(in) :: vg1(2,nfft*optn*opteltfr),vg1_core(2,nfft*optn*opteltfr)
 real(dp),intent(in) :: vprtrb(2),vspl(mqgrid,2,ntypat*optv)
 real(dp),intent(out) :: atmrho(nfft*optn)  !vz_i
 real(dp),intent(inout) :: atmvloc(nfft*optv) !vz_i
 real(dp),intent(out) :: dyfrn(3,3,natom*optn*optdyfr),dyfrv(3,3,natom*optv*optdyfr)
 real(dp),intent(out) :: eltfrn(6+3*natom,6)
 real(dp),intent(inout) :: grn(3,natom*optn*optgr) !vz_i
 real(dp),intent(out) :: grv(3,natom*optv*optgr),strn(6*optn*optstr)
 real(dp),intent(out) :: strv(6*optv*optstr)
 type(pawtab_type),target,intent(in) :: pawtab(ntypat*usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,ia,ia1,ia2,id1,id2,id3,ierr,ig1,ig1_,ig2,ig2_,ig3,ig3_,ii,is,is1,is2
 integer :: itypat,jj,js,ka,kb,kd,kg,me_fft,my_comm_fft,ndir,n1,n2,n3,nproc_fft,paral_kgb_fft
 integer :: shift1,shift2,shift3
 logical :: have_g0
 real(dp),parameter :: tolfix=1.0000001_dp
 real(dp) :: aa,alf2pi2,bb,cc,cutoff,dbl_ig1,dbl_ig2,dbl_ig3,dd,dg1,dg2,d2g,diff
 real(dp) :: dn_at,d2n_at,d2n_at2,dq,dq2div6,dqdiv6,dqm1,dv_at,ee,ff,gauss1,gauss2,gg,gmag,gsquar,n_at
 real(dp) :: ph12i,ph12r,ph1i,ph1r,ph2i,ph2r,ph3i,ph3r,sfi,sfr,term,term1,term2,tmpni,tmpnr
 real(dp) :: tmpvi,tmpvr,v_at,xnorm
 character(len=500) :: message
 type(distribfft_type),pointer :: my_distribfft
 type(mpi_type) :: mpi_enreg_fft
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 real(dp), ABI_CONTIGUOUS pointer :: tvalespl(:,:),tcorespl(:,:)
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer  :: delta(6)=(/1,1,1,0,0,0/)
 real(dp) :: dgm(3,3,6),d2gm(3,3,6,6),gcart(3),tsec(2)
 real(dp),allocatable :: dyfrn_indx(:,:,:),dyfrv_indx(:,:,:),grn_indx(:,:)
 real(dp),allocatable :: grv_indx(:,:),phim_igia(:),phre_igia(:),workn(:,:)
 real(dp),allocatable :: workv(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!Check optional arguments
 if (present(comm_fft)) then
   if ((.not.present(paral_kgb)).or.(.not.present(me_g0))) then
     MSG_BUG(' Need paral_kgb and me_g0 with comm_fft !')
   end if
 end if

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 
 me_fft=ngfft(11)
 nproc_fft=ngfft(10)


!Get the distrib associated with this fft_grid
 if (present(distribfft)) then
   my_distribfft => distribfft
 else
   ABI_DATATYPE_ALLOCATE(my_distribfft,)
   call init_distribfft_seq(my_distribfft,'f',n2,n3,'fourdp')
 end if
 if (n2==my_distribfft%n2_coarse) then
   fftn2_distrib => my_distribfft%tab_fftdp2_distrib
   ffti2_local => my_distribfft%tab_fftdp2_local
 else if (n2 == my_distribfft%n2_fine) then
   fftn2_distrib => my_distribfft%tab_fftdp2dg_distrib
   ffti2_local => my_distribfft%tab_fftdp2dg_local
 else
   MSG_BUG("Unable to find an allocated distrib for this fft grid")
 end if

 if (present(is2_in)) then
   if(is2_in<1.or.is2_in>6) then
     MSG_BUG("is2_in must be between 1 and 6")
   else
     ndir = 1
   end if
   ndir = -1
 end if

!Zero out arrays to permit accumulation over atom types
 if (optv==1.and.optatm==1) then
   ABI_ALLOCATE(workv,(2,nfft))
   workv(:,:)=zero
 end if
 if (optn==1.and.optatm==1) then
   ABI_ALLOCATE(workn,(2,nfft))
   workn(:,:)=zero
 end if
 if (optv==1.and.optgr==1) then
   ABI_ALLOCATE(grv_indx,(3,natom))
   grv_indx(:,:)=zero
 end if
 if (optn==1.and.optgr==1) then
   ABI_ALLOCATE(grn_indx,(3,natom))
   grn_indx(:,:)=zero
 end if
 if (optv==1.and.optdyfr==1) then
   ABI_ALLOCATE(dyfrv_indx,(3,3,natom))
   dyfrv_indx(:,:,:)=zero
 end if
 if (optn==1.and.optdyfr==1) then
   ABI_ALLOCATE(dyfrn_indx,(3,3,natom))
   dyfrn_indx(:,:,:)=zero
 end if
 if (optv==1.and.optstr==1) strv(:)=zero
 if (optn==1.and.optstr==1) strn(:)=zero
 if (opteltfr==1) eltfrn(:,:) = zero

!Compute 1st and 2nd derivatives of metric tensor wrt all strain components
!and store for use in inner loop below for elastic tensor.
 if (opteltfr==1) then
   dgm(:,:,:)=zero
   d2gm(:,:,:,:)=zero
!  Loop over 2nd strain index
   do is2=1,6
     kg=idx(2*is2-1);kd=idx(2*is2)
     do jj = 1,3
       dgm(:,jj,is2)=-(gprimd(kg,:)*gprimd(kd,jj)+gprimd(kd,:)*gprimd(kg,jj))
     end do

!    Loop over 1st strain index, upper triangle only
     do is1=1,is2
       ka=idx(2*is1-1);kb=idx(2*is1)
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
     end do !is1
   end do !is2
 end if


 dq=(qgrid(mqgrid)-qgrid(1))/dble(mqgrid-1)
 dqm1=1.0_dp/dq
 dqdiv6=dq/6.0_dp
 dq2div6=dq**2/6.0_dp
 cutoff=gsqcut*tolfix
 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2

 ABI_ALLOCATE(phre_igia,(natom))
 ABI_ALLOCATE(phim_igia,(natom))

 ia1=1
 do itypat=1,ntypat
!  ia1,ia2 sets range of loop over atoms:
   ia2=ia1+nattyp(itypat)-1
   ii=0

   if (optn2==3)then
     gauss1=gauss(1,itypat)
     gauss2=gauss(2,itypat)
     alf2pi2=(two_pi*gauss2)**2
   end if

   if (usepaw == 1) then
     tcorespl => pawtab(itypat)%tcorespl
     tvalespl => pawtab(itypat)%tvalespl
   else
     tcorespl => psps%nctab(itypat)%tcorespl
     tvalespl => psps%nctab(itypat)%tvalespl
   end if
   
   do i3=1,n3
     ig3=i3-(i3/id3)*n3-1
     ig3_=ig3;if (ig3_==(n3/2+1)) ig3_=0
     do i2=1,n2
       ig2=i2-(i2/id2)*n2-1
       ig2_=ig2;if (ig2_==(n2/2+1)) ig2_=0
       if(fftn2_distrib(i2)==me_fft) then
         do i1=1,n1
           ig1=i1-(i1/id1)*n1-1
           ig1_=ig1;if (ig1_==(n1/2+1)) ig1_=0
           ii=ii+1
           gsquar=gsq_atm(ig1,ig2,ig3)

!          Skip G**2 outside cutoff:
           if (gsquar<=cutoff) then

             gmag=sqrt(gsquar)
             have_g0=(ig1==0.and.ig2==0.and.ig3==0)

             jj=1+int(gmag*dqm1)
             diff=gmag-qgrid(jj)

!            Compute structure factor for all atoms of given type:
             do ia=ia1,ia2
               shift1=1+n1+(ia-1)*(2*n1+1)
               shift2=1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1)
               shift3=1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
               ph1r=ph1d(1,ig1+shift1);ph1i=ph1d(2,ig1+shift1)
               ph2r=ph1d(1,ig2+shift2);ph2i=ph1d(2,ig2+shift2)
               ph3r=ph1d(1,ig3+shift3);ph3i=ph1d(2,ig3+shift3)
               ph12r=ph1r*ph2r-ph1i*ph2i
               ph12i=ph1r*ph2i+ph1i*ph2r
               phre_igia(ia)=ph12r*ph3r-ph12i*ph3i
               phim_igia(ia)=ph12r*ph3i+ph12i*ph3r
             end do

!            Assemble structure factors for this type of atom= sum[exp(-i.piG.R)]
             if (optatm==1.or.optstr==1.or.opteltfr==1) then
               sfr=zero;sfi=zero
               do ia=ia1,ia2
                 sfr=sfr+phre_igia(ia)
                 sfi=sfi-phim_igia(ia)
               end do
             end if

!            Compute V^AT(G) and/or n^AT(G) for given type of atom
!            Evaluate spline fit: p. 86 Numerical Recipes, Press et al;
!            NOTE: error in book for sign of "aa" term in derivative;
!            !           also see splfit routine.
             if (optv==1.or.optn2/=3) then
               bb = diff*dqm1
               aa = 1.0_dp-bb
               cc = aa*(aa**2-1.0_dp)*dq2div6
               dd = bb*(bb**2-1.0_dp)*dq2div6
             end if
             if (optv==1) then
               if (have_g0) then
                 v_at=zero
               else
                 v_at=(aa*vspl(jj,1,itypat)+bb*vspl(jj+1,1,itypat)+&
&                 cc*vspl(jj,2,itypat)+dd*vspl(jj+1,2,itypat)) &
&                 /gsquar
               end if
             end if
             if (optn==1) then
               if (optn2==1) then
                 n_at=(aa*tcorespl(jj,1)+bb*tcorespl(jj+1,1)+cc*tcorespl(jj,2)+dd*tcorespl(jj+1,2))
               else if (optn2==2) then
                 n_at=(aa*tvalespl(jj,1)+bb*tvalespl(jj+1,1)+cc*tvalespl(jj,2)+dd*tvalespl(jj+1,2))
               else if (optn2==3) then
                 n_at=gauss1*exp(-gsquar*alf2pi2)
               else
                 n_at=zero
               end if
             end if

!            Compute sum of local atomic potentials or densities
!            ---------------------------------------------------
             if(optatm==1) then
!              Accumulate V^AT(G)*SF(G) or n^AT(G)*SF(G)
               if (optv==1) then
                 workv(re,ii)=workv(re,ii)+sfr*v_at
                 workv(im,ii)=workv(im,ii)+sfi*v_at
               end if
               if (optn==1) then
                 workn(re,ii)=workn(re,ii)+sfr*n_at
                 workn(im,ii)=workn(im,ii)+sfi*n_at
               end if

!            Compute contrib. to forces and/or frozen part of dyn. matrix
!            -------------------------------------------------------------
             else if (optgr==1.or.optdyfr==1) then
               dbl_ig1=dble(ig1_);dbl_ig2=dble(ig2_);dbl_ig3=dble(ig3_)
!              Compute (2Pi)*V^AT(G)*rho(G) or (2Pi)*n^AT(G)*V(G)
               if (optv==1) then
                 tmpvr=(two_pi*v_at)*rhog(re,ii)
                 tmpvi=(two_pi*v_at)*rhog(im,ii)
               end if
               if (optn==1) then
                 tmpnr=(two_pi*n_at)*vg(re,ii)
                 tmpni=(two_pi*n_at)*vg(im,ii)
               end if
!              === contrib. to forces
               if (optgr==1) then
!                Accumulate -(2Pi.G)*V^AT(G)*rho(G)*SF(G)
!                or -(2Pi)*n^AT(G)*V(G)*SF(G) into forces
                 if (optv==1) then
                   do ia=ia1,ia2
                     term=tmpvi*phre_igia(ia)+tmpvr*phim_igia(ia)
                     grv_indx(1,ia)=grv_indx(1,ia)-dbl_ig1*term
                     grv_indx(2,ia)=grv_indx(2,ia)-dbl_ig2*term
                     grv_indx(3,ia)=grv_indx(3,ia)-dbl_ig3*term
                   end do
                 end if
                 if (optn==1) then
                   do ia=ia1,ia2
                     term=tmpni*phre_igia(ia)+tmpnr*phim_igia(ia)
                     grn_indx(1,ia)=grn_indx(1,ia)-dbl_ig1*term
                     grn_indx(2,ia)=grn_indx(2,ia)-dbl_ig2*term
                     grn_indx(3,ia)=grn_indx(3,ia)-dbl_ig3*term
                   end do
                 end if
               end if
!              === contrib. to frozen part of dyn. matrix
               if (optdyfr==1) then
!                Accumulate -(2Pi^2.Gi.Gj)*V^AT(G)*rho(G)*SF(G)
!                or -(2Pi^2.Gi.Gj)*n^AT(G)*V(G)*SF(G) into dyn. matrix
                 if (optv==1) then
                   do ia=ia1,ia2
                     term=two_pi*(tmpvr*phre_igia(ia)-tmpvi*phim_igia(ia))
                     dyfrv_indx(1,1,ia)=dyfrv_indx(1,1,ia)-dbl_ig1*dbl_ig1*term
                     dyfrv_indx(1,2,ia)=dyfrv_indx(1,2,ia)-dbl_ig1*dbl_ig2*term
                     dyfrv_indx(1,3,ia)=dyfrv_indx(1,3,ia)-dbl_ig1*dbl_ig3*term
                     dyfrv_indx(2,2,ia)=dyfrv_indx(2,2,ia)-dbl_ig2*dbl_ig2*term
                     dyfrv_indx(2,3,ia)=dyfrv_indx(2,3,ia)-dbl_ig2*dbl_ig3*term
                     dyfrv_indx(3,3,ia)=dyfrv_indx(3,3,ia)-dbl_ig3*dbl_ig3*term
                   end do
                 end if
                 if (optn==1) then
                   do ia=ia1,ia2
                     term=two_pi*(tmpnr*phre_igia(ia)-tmpni*phim_igia(ia))
                     dyfrn_indx(1,1,ia)=dyfrn_indx(1,1,ia)-dbl_ig1*dbl_ig1*term
                     dyfrn_indx(1,2,ia)=dyfrn_indx(1,2,ia)-dbl_ig1*dbl_ig2*term
                     dyfrn_indx(1,3,ia)=dyfrn_indx(1,3,ia)-dbl_ig1*dbl_ig3*term
                     dyfrn_indx(2,2,ia)=dyfrn_indx(2,2,ia)-dbl_ig2*dbl_ig2*term
                     dyfrn_indx(2,3,ia)=dyfrn_indx(2,3,ia)-dbl_ig2*dbl_ig3*term
                     dyfrn_indx(3,3,ia)=dyfrn_indx(3,3,ia)-dbl_ig3*dbl_ig3*term
                   end do
                 end if
               end if
             end if

!            Compute (dV^AT(q)/dq)/q and/or (dn^AT(q)/dq)/q
!            For stress tensor and/or elastic tensor
!            ---------------------------------
             if (optstr==1.or.opteltfr==1) then
!              Note: correction of Numerical Recipes sign error before (3._dp*aa**2-1._dp)
!              ee*dqm1 + ff*dqdiv6 is the best estimate of dV(q)/dq from splines
               if (optv==1) then
                 if (have_g0) then
                   dv_at=zero
                 else
                   ee=vspl(jj+1,1,itypat)-vspl(jj,1,itypat)
                   ff=(3._dp*bb**2-1._dp)*vspl(jj+1,2,itypat)&
&                   -(3._dp*aa**2-1._dp)*vspl(jj  ,2,itypat)
                   dv_at=((ee*dqm1+ff*dqdiv6)/gmag-2.0_dp*v_at)/gsquar
                 end if
               end if
               if (optn==1) then
                 if (have_g0) then
                   if (optn2==1) then
                     if (usepaw ==1) then
                       dn_at=pawtab(itypat)%dncdq0
                     else
                       dn_at=psps%nctab(itypat)%dncdq0
                     end if
                   else if (optn2==2) then
                     if (usepaw == 1) then
                       dn_at=pawtab(itypat)%dnvdq0
                     else
                       dn_at=psps%nctab(itypat)%dnvdq0
                     end if
                   else if (optn2==3) then
                     dn_at=-two*gauss1*alf2pi2
                   end if
                   if (opteltfr==1) then
                     d2n_at = 0
                   end if
                 else
                   if (optn2==1) then
                     ee=tcorespl(jj+1,1)-tcorespl(jj,1)
                     ff=(3._dp*bb**2-1._dp)*tcorespl(jj+1,2) &
&                     -(3._dp*aa**2-1._dp)*tcorespl(jj,2)
!                    Also get nc''(q)
                     if (opteltfr==1) then
                       gg=aa*tcorespl(jj,2)+bb*tcorespl(jj+1,2)
                     end if
                   else if (optn2==2) then
                     ee=tvalespl(jj+1,1)-tvalespl(jj,1)
                     ff=(3._dp*bb**2-1._dp)*tvalespl(jj+1,2) &
&                     -(3._dp*aa**2-1._dp)*tvalespl(jj,2)
!                    Also get nc''(q)
                     if (opteltfr==1) then
                       gg=aa*tvalespl(jj,2)+bb*tvalespl(jj+1,2)
                     end if
                   else if (optn2==3) then
                     dn_at=-two*gauss1*alf2pi2*exp(-gsquar*alf2pi2)
                   else
                   end if
                   dn_at  = (ee*dqm1+ff*dqdiv6)/gmag
                   if (opteltfr==1) then
                     d2n_at = (gg-dn_at)/gsquar
                     d2n_at2 = gg/gmag**3
                   end if
                 end if
               end if

!              Compute G in cartesian coordinates
               gcart(1)=gprimd(1,1)*dble(ig1)+gprimd(1,2)*dble(ig2)+&
&               gprimd(1,3)*dble(ig3)
               gcart(2)=gprimd(2,1)*dble(ig1)+gprimd(2,2)*dble(ig2)+&
&               gprimd(2,3)*dble(ig3)
               gcart(3)=gprimd(3,1)*dble(ig1)+gprimd(3,2)*dble(ig2)+&
&               gprimd(3,3)*dble(ig3)
             end if

!            Compute contrib. to stress tensor
!            ---------------------------------
             if (optstr==1)then
!              Accumulate -dV^AT/dG*rho(G)*SF(G)*Gi.Gj/G
!              or -dn^AT/dG*V(G)*SF(G)*Gi.Gj/G
!              into stress tensor
               if (optv==1) then
                 term=(rhog(re,ii)*sfr+rhog(im,ii)*sfi)
                 strv(1)=strv(1)-term*(dv_at*gcart(1)*gcart(1)+v_at)
                 strv(2)=strv(2)-term*(dv_at*gcart(2)*gcart(2)+v_at)
                 strv(3)=strv(3)-term*(dv_at*gcart(3)*gcart(3)+v_at)
                 strv(4)=strv(4)-term*dv_at*gcart(3)*gcart(2)
                 strv(5)=strv(5)-term*dv_at*gcart(3)*gcart(1)
                 strv(6)=strv(6)-term*dv_at*gcart(2)*gcart(1)
               end if
               if (optn==1) then
                 term=(vg(re,ii)*sfr+vg(im,ii)*sfi)*dn_at
                 strn(1)=strn(1)-term*gcart(1)*gcart(1)
                 strn(2)=strn(2)-term*gcart(2)*gcart(2)
                 strn(3)=strn(3)-term*gcart(3)*gcart(3)
                 strn(4)=strn(4)-term*gcart(3)*gcart(2)
                 strn(5)=strn(5)-term*gcart(3)*gcart(1)
                 strn(6)=strn(6)-term*gcart(2)*gcart(1)
               end if
             end if

!            Compute contrib. to elastic tensor
!            ---------------------------------
             if (opteltfr==1) then
               dbl_ig1=dble(ig1_);dbl_ig2=dble(ig2_);dbl_ig3=dble(ig3_)
!             if (ig1==0 .and. ig2==0 .and. ig3==0) cycle
!              Compute G*dG/Deps_{\gamme\delta}
               dg2=0.5_dp*dgsqds_atm(ig1,ig2,ig3,is2_in)

               term  = vg (re,ii)*sfr + vg (im,ii)*sfi
               term1 = vg1(re,ii)*sfr + vg1(im,ii)*sfi
               term2 = vg1_core(re,ii)*sfr + vg1_core(im,ii)*sfi
               
               do is1=1,6

!                Compute G*dG/Deps_{\alpha\beta}
                 dg1=0.5_dp*dgsqds_atm(ig1,ig2,ig3,is1)
!                Compute G^2*d2G/Deps_{alphabetagammadelta}
                 d2g=(0.25_dp*d2gsqds_atm(ig1,ig2,ig3,is1,is2_in))

                 eltfrn(is1,is2_in) = eltfrn(is1,is2_in) + (term*(d2n_at*dg1*dg2 + d2g*dn_at))
                 eltfrn(is1,is2_in) = eltfrn(is1,is2_in) + 0.5*(term1*dn_at*dg1)
                 eltfrn(is2_in,is1) = eltfrn(is2_in,is1) + 0.5*(term1*dn_at*dg1)

                 if(is2_in<=3)then
                   eltfrn(is1,is2_in) = eltfrn(is1,is2_in) - term*dn_at*dg1
                 end if
                 if(is1<=3)then
                   eltfrn(is1,is2_in) = eltfrn(is1,is2_in) - term*dn_at*dg2
                 end if
                 if(is2_in<=3.and.is1<=3)then
                   eltfrn(is1,is2_in) = eltfrn(is1,is2_in) - (term1-term)*n_at
                 end if
               end do

!              internal strain
               do ia=ia1,ia2
                 js=7+3*(ia-1)
!                Compute -2pi*G*i*vxcis_core(G)*nat(G)*(exp(-iGr))
                 term=(vg1_core(im,ii)*phre_igia(ia)+vg1_core(re,ii)*phim_igia(ia))*n_at*two_pi
                 eltfrn(js  ,is2_in) = eltfrn(js  ,is2_in) - dbl_ig1*term
                 eltfrn(js+1,is2_in) = eltfrn(js+1,is2_in) - dbl_ig2*term
                 eltfrn(js+2,is2_in) = eltfrn(js+2,is2_in) - dbl_ig3*term

!                Compute -2pi*G*i*vxc(G)*(dnat*dG/deps-delta*nat(G))*(exp(-iGr))
                 term=(vg(im,ii)*phre_igia(ia)+vg(re,ii)*phim_igia(ia))*&
&                 (dn_at*dg2-delta(is2_in)*n_at)*two_pi
                 eltfrn(js  ,is2_in) = eltfrn(js  ,is2_in) - dbl_ig1*term 
                 eltfrn(js+1,is2_in) = eltfrn(js+1,is2_in) - dbl_ig2*term 
                 eltfrn(js+2,is2_in) = eltfrn(js+2,is2_in) - dbl_ig3*term 
               end do
             end if
!            End skip G**2 outside cutoff:
           end if
!          End loop on n1, n2, n3
         end do
       end if ! this plane is for me_fft
     end do
   end do

!  Symmetrize the dynamical matrix with respect to indices
   if (optdyfr==1) then
     if (optv==1) then
       do ia=ia1,ia2
         dyfrv_indx(2,1,ia)=dyfrv_indx(1,2,ia)
         dyfrv_indx(3,1,ia)=dyfrv_indx(1,3,ia)
         dyfrv_indx(3,2,ia)=dyfrv_indx(2,3,ia)
       end do
     end if
     if (optn==1) then
       do ia=ia1,ia2
         dyfrn_indx(2,1,ia)=dyfrn_indx(1,2,ia)
         dyfrn_indx(3,1,ia)=dyfrn_indx(1,3,ia)
         dyfrn_indx(3,2,ia)=dyfrn_indx(2,3,ia)
       end do
     end if
   end if

   ia1=ia2+1

!  End loop on type of atoms
 end do

 ABI_DEALLOCATE(phre_igia)
 ABI_DEALLOCATE(phim_igia)

!Get local potential or density back to real space
 if(optatm==1)then
!  Allow for the addition of a perturbing potential
   if ((optv==1).and.(vprtrb(1)**2+vprtrb(2)**2) > 1.d-30) then
!    Find the linear indices which correspond with the input wavevector qprtrb
!    The double modulus handles both i>=n and i<0, mapping into [0,n-1];
!    then add 1 to get range [1,n] for each
     i3=1+mod(n3+mod(qprtrb(3),n3),n3)
     i2=1+mod(n2+mod(qprtrb(2),n2),n2)
     i1=1+mod(n1+mod(qprtrb(1),n1),n1)
!    Compute the linear index in the 3 dimensional array
     ii=i1+n1*((ffti2_local(i2)-1)+(n2/nproc_fft)*(i3-1))
!    Add in the perturbation at G=qprtrb
     workv(re,ii)=workv(re,ii)+0.5_dp*vprtrb(1)
     workv(im,ii)=workv(im,ii)+0.5_dp*vprtrb(2)
!    Same thing for G=-qprtrb
     i3=1+mod(n3+mod(-qprtrb(3),n3),n3)
     i2=1+mod(n2+mod(-qprtrb(2),n2),n2)
     i1=1+mod(n1+mod(-qprtrb(1),n1),n1)
!    ii=i1+n1*((i2-1)+n2*(i3-1))
     workv(re,ii)=workv(re,ii)+0.5_dp*vprtrb(1)
     workv(im,ii)=workv(im,ii)-0.5_dp*vprtrb(2)
     write(message, '(a,1p,2e12.4,a,0p,3i4,a)' )&
&     ' atm2fft: perturbation of vprtrb=', vprtrb,&
&     ' and q=',qprtrb,' has been added'
     call wrtout(std_out,message,'COLL')
   end if

   if (optv==1.or.optn==1) then
!    Create fake mpi_enreg to wrap fourdp
     call initmpi_seq(mpi_enreg_fft)
     ABI_DATATYPE_DEALLOCATE(mpi_enreg_fft%distribfft)
     if (present(comm_fft)) then
       call set_mpi_enreg_fft(mpi_enreg_fft,comm_fft,my_distribfft,me_g0,paral_kgb)
       my_comm_fft=comm_fft;paral_kgb_fft=paral_kgb
     else
       my_comm_fft=xmpi_comm_self;paral_kgb_fft=0;
       mpi_enreg_fft%distribfft => my_distribfft
     end if
!    Non-symetrized non-zero elements have to be nullified
!    Transform back to real space; divide by unit cell volume
     xnorm=one/ucvol
     if (optv==1) then
       call zerosym(workv,2,n1,n2,n3,comm_fft=my_comm_fft,distribfft=my_distribfft)
       call fourdp(1,workv,atmvloc,1,mpi_enreg_fft,nfft,ngfft,paral_kgb_fft,0)
       atmvloc(:)=atmvloc(:)*xnorm
       ABI_DEALLOCATE(workv)
     end if
     if (optn==1) then
       call zerosym(workn,2,n1,n2,n3,comm_fft=my_comm_fft,distribfft=my_distribfft)
       call fourdp(1,workn,atmrho,1,mpi_enreg_fft,nfft,ngfft,paral_kgb_fft,0)
       atmrho(:)=atmrho(:)*xnorm
       ABI_DEALLOCATE(workn)
     end if
!    Destroy fake mpi_enreg
     call unset_mpi_enreg_fft(mpi_enreg_fft)
   end if

 end if

!Additional treatment in case of parallelization
 if (present(comm_fft)) then
   if ((xmpi_comm_size(comm_fft)>1).and.&
&   (optgr==1.or.optstr==1.or.optdyfr==1.or.opteltfr==1)) then
     call timab(48,1,tsec)
     if (optv==1) then
       if (optgr==1)then
         call xmpi_sum(grv_indx,comm_fft,ierr)
       end if
       if (optstr==1)then
         call xmpi_sum(strv,comm_fft,ierr)
       end if
       if (optdyfr==1)then
         call xmpi_sum(dyfrv_indx,comm_fft,ierr)
       end if
     end if
     if (optn==1) then
       if (optgr==1)then
         call xmpi_sum(grn_indx,comm_fft,ierr)
       end if
       if (optstr==1)then
         call xmpi_sum(strn,comm_fft,ierr)
       end if
       if (optdyfr==1)then
         call xmpi_sum(dyfrn_indx,comm_fft,ierr)
       end if
     end if
     call timab(48,2,tsec)
   end if
 end if

!Forces: re-order atoms
 if (optgr==1) then
   if (optv==1) then
     do ia=1,natom
       grv(1:3,atindx1(ia))=grv_indx(1:3,ia)
     end do
     ABI_DEALLOCATE(grv_indx)
   end if
   if (optn==1) then
     do ia=1,natom
       grn(1:3,atindx1(ia))=grn_indx(1:3,ia)
     end do
     ABI_DEALLOCATE(grn_indx)
   end if
 end if

!Elastic tensor: Fill in lower triangle
! if (opteltfr==1) then
!   do is2=2,6
!     do is1=1,is2-1
!       eltfrn(is2,is1)=eltfrn(is1,is2)
!     end do
!   end do
! end if

!Normalize stress tensor: 
 if (optstr==1) then
   if (optv==1) then
     strv(:)=strv(:)/ucvol
   end if
   if (optn==1) strn(:)=strn(:)/ucvol
 end if

!Dynamical matrix: re-order atoms
 if (optdyfr==1) then
   if (optv==1) then
     do ia=1,natom
       dyfrv(1:3,1:3,atindx1(ia))=dyfrv_indx(1:3,1:3,ia)
     end do
     ABI_DEALLOCATE(dyfrv_indx)
   end if
   if (optn==1) then
     do ia=1,natom
       dyfrn(1:3,1:3,atindx1(ia))=dyfrn_indx(1:3,1:3,ia)
     end do
     ABI_DEALLOCATE(dyfrn_indx)
   end if
 end if

 if (.not.present(distribfft)) then
   call destroy_distribfft(my_distribfft)
   ABI_DATATYPE_DEALLOCATE(my_distribfft)
 end if

 DBG_EXIT("COLL")

 contains

   function gsq_atm(i1,i2,i3)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gsq_atm'
!End of the abilint section

   real(dp) :: gsq_atm
   integer,intent(in) :: i1,i2,i3
   gsq_atm=dble(i1*i1)*gmet(1,1)+dble(i2*i2)*gmet(2,2)+dble(i3*i3)*gmet(3,3) &
&   +two*(dble(i1*i2)*gmet(1,2)+dble(i2*i3)*gmet(2,3)+dble(i3*i1)*gmet(3,1))
 end function gsq_atm

 function dgsqds_atm(i1,i2,i3,is)
!Define dG^2/ds based on G space metric derivative

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dgsqds_atm'
!End of the abilint section

   real(dp) :: dgsqds_atm
   integer,intent(in) :: i1,i2,i3,is
   dgsqds_atm=dble(i1*i1)*dgm(1,1,is)+dble(i2*i2)*dgm(2,2,is)+&
&   dble(i3*i3)*dgm(3,3,is)+&
&   dble(i1*i2)*(dgm(1,2,is)+dgm(2,1,is))+&
&   dble(i1*i3)*(dgm(1,3,is)+dgm(3,1,is))+&
&   dble(i2*i3)*(dgm(2,3,is)+dgm(3,2,is))
 end function dgsqds_atm

 function d2gsqds_atm(i1,i2,i3,is1,is2)
!  Define 2dG^2/ds1ds2  based on G space metric derivative

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'd2gsqds_atm'
!End of the abilint section

   real(dp) :: d2gsqds_atm
   integer,intent(in) :: i1,i2,i3,is1,is2
   d2gsqds_atm=dble(i1*i1)*d2gm(1,1,is1,is2)+&
&   dble(i2*i2)*d2gm(2,2,is1,is2)+dble(i3*i3)*d2gm(3,3,is1,is2)+&
&   dble(i1*i2)*(d2gm(1,2,is1,is2)+d2gm(2,1,is1,is2))+&
&   dble(i1*i3)*(d2gm(1,3,is1,is2)+d2gm(3,1,is1,is2))+&
&   dble(i2*i3)*(d2gm(2,3,is1,is2)+d2gm(3,2,is1,is2))
 end function d2gsqds_atm

end subroutine atm2fft
!!***
