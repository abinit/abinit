!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_atm2fft
!! NAME
!! dfpt_atm2fft
!!
!! FUNCTION
!! This routine sums 1st-order atomic functions (density or potential)
!! defined (in rec. space) on a radial grid to get global 1st-order
!! quantities on the fine FFT grid.
!!
!! Possible options:
!!   optn=1: compute a sum of local 1st-order atomic densities
!!   optv=1: compute a sum of local 1st-order atomic potentials
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms ordered by type
!!  cplex: if 1, real space 1-order functions on FFT grid
!!  distribfft<type(distribfft_type)>=--optional-- contains infos related to FFT parallelism
!!  eei=local pseudopotential part of total energy
!!  gauss(2,ntypat)= params for gaussian atm density (optn2=3) for each atom type
!!  gmet(3,3)=reciprocal space metric
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  gsqcut=cutoff on |G|^2: see setup1 for definition (doubled sphere)
!!  idir=direction of atomic displacement (in case of phonons perturb.)
!!       used only if ndir=1 (see below)
!!  ipert=nindex of perturbation
!!  me_g0=--optional-- 1 if the current process treat the g=0 plane-wave (only needed when comm_fft is present)
!!  mgfft=maximum size of 1D FFTs
!!  comm_fft=--optional-- MPI communicator over FFT components
!!  mqgrid=number of grid pts in q array for f(q) spline.
!!  natom=number of atoms in unit cell.
!!  ndir=number of directions of atomic displacement (in case of phonon):
!!       can be 1 (idir direction in then used) or 3 (all directions)
!!       6 cartesian strain component 11,22,33,32,31,21 (in case strain perturbation)
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT
!!  nattyp(ntypat)=array describing how many atoms of each type in cell 
!!  ntypat=number of types of atoms.
!!  optn,optn2,optv= (see NOTES below)
!!  paral_kgb=--optional-- 1 if "band-FFT" parallelism is activated (only needed when comm_fft is present)
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  qgrid(mqgrid)=q grid for spline from 0 to qmax.
!!  qphon(3)=wavevector of the phonon
!!  typat(natom)=type of each atom
!!  ucvol=unit cell volume
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vspl(mqgrid,2,ntypat)=q^2 v(q) spline of an atomic potential
!!                        (used only if optv=1)
!!  xred(3,natom)=reduced atomic coordinates
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!! OUTPUT
!!  ======= if optv==1 =======
!!    atmvlocr1(cplex*nfft,ndir)=sum of local 1st-order atomic potentials in real space
!!    atmvlocg1(2,nfft,ndir)=sum of local 1st-order atomic potentials in G-space
!!  ======= if optn==1 =======
!!   --- if optatm==1
!!    atmrhor1(cplex*nfft,ndir)=sum of 1st-order atomic densities in real space
!!    atmrhog1(2,nfft,ndir)=sum of 1st-order atomic densities in G-space
!!
!! NOTES
!! Details on possible options:
!! ============================
!! optv: controls the computation of a local 1st-order potential as sum of atomic potentials
!!          Vloc(r)=Sum_R[V1^AT(r-R)]
!! optn: controls the computation of a 1st-order density as sum of atomic densities
!!          n(r)=Sum_R[n1^AT(r-R)]
!!          n^AT is stored in reciprocal space:
!!          if optn2=1: n^AT is the atomic PAW PS core density stored in array pawtab%tcorespl()
!!                   2: n^AT is the atomic PAW PS valence density stored in array pawtab%tvalespl()
!!                   3: n^AT is a gaussian density: n(g)=gauss(1,ityp)*exp[-(gauss(2,ityp)*G)^2]
!! Note: optv and optn can be activated together
!!
!! Typical uses:
!! =============
!! Computation of:
!!  - 1st-order local potential: optv=1
!!  - 1st-order PS core density: optn=1, optn2=1
!!
!! PARENTS
!!      dfpt_dyxc1,dfpt_eltfrxc,dfpt_looppert,dfpt_nstpaw,pawgrnl
!!
!! CHILDREN
!!      destroy_distribfft,fourdp,init_distribfft_seq,initmpi_seq
!!      set_mpi_enreg_fft,unset_mpi_enreg_fft,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dfpt_atm2fft(atindx,cplex,gmet,gprimd,gsqcut,idir,ipert,&
&                   mgfft,mqgrid,natom,ndir,nfft,ngfft,ntypat,&
&                   ph1d,qgrid,qphon,typat,ucvol,usepaw,xred,psps,pawtab,&
&                   atmrhor1,atmrhog1,atmvlocr1,atmvlocg1,distribfft,gauss,comm_fft,me_g0,optn_in,&
&                   optn2_in,optv_in,paral_kgb,vspl) ! optional arguments

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_xmpi,       only : xmpi_comm_self, xmpi_comm_size
 use defs_datatypes,only : pseudopotential_type
 use m_pawtab,     only : pawtab_type
 use m_distribfft, only : distribfft_type
 use m_mpinfo,     only : set_mpi_enreg_fft, unset_mpi_enreg_fft

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_atm2fft'
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mgfft,mqgrid,natom,ndir,nfft,ntypat,usepaw
 integer,optional,intent(in) :: optn_in,optn2_in,optv_in
 integer,optional,intent(in) :: me_g0,comm_fft,paral_kgb
 real(dp),intent(in) :: gsqcut,ucvol
 type(pseudopotential_type),target,intent(in) :: psps
 type(distribfft_type),optional,intent(in),target :: distribfft
!arrays
 integer,intent(in) :: atindx(natom),ngfft(18),typat(natom)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qgrid(mqgrid),qphon(3)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),optional,intent(in)  :: gauss(2,ntypat),vspl(mqgrid,2,ntypat)
 real(dp),optional,intent(out) :: atmrhor1(cplex*nfft,ndir)
 real(dp),optional,intent(out) :: atmrhog1(2,nfft,ndir)
 real(dp),optional,intent(out) :: atmvlocr1(cplex*nfft,ndir)
 real(dp),optional,intent(out) :: atmvlocg1(2,nfft,ndir)
 type(pawtab_type),target,intent(in) :: pawtab(ntypat*usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,ia,ia1,ia2,iatm,iatom,id,id1,id2,id3
 integer :: ig1,ig1max,ig1min,ig2,ig2max,ig2min,ig3,ig3max,ig3min
 integer :: ii,itypat,jj,me_fft,my_comm_fft,n1,n2,n3,nattyp,nproc_fft,ntype,paral_kgb_fft
 integer :: optn,optv,optn2,shift1,shift2,shift3,type1,type2
 logical :: have_g0,qeq0,qeq05
 real(dp),parameter :: tolfix=1.0000001_dp
 real(dp) :: aa,alf2pi2,bb,cc,cutoff,dd,diff,dq,dq2div6,dqdiv6,dqm1,ee,ff
 real(dp) :: gauss1,gauss2,gmag,gq1,gq2,gq3,gsquar,n_at,dn_at,ph12i,ph12r,ph1i
 real(dp) :: ph1r,ph2i,ph2r,ph3i,ph3r,phqim,phqre,qxred2pi
 real(dp) :: sfi,sfqi,sfqr,sfr,term_n,term_v,v_at,dv_at,xnorm
 type(distribfft_type),pointer :: my_distribfft
 type(mpi_type) :: mpi_enreg_fft
!arrays
 integer :: eps1(6)=(/1,2,3,2,3,1/),eps2(6)=(/1,2,3,3,1,2/),jdir(ndir)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:)
 real(dp), ABI_CONTIGUOUS pointer :: tvalespl(:,:),tcorespl(:,:)
 real(dp) ::  gq(6),gcart(3)
 real(dp),allocatable :: phim_igia(:),phre_igia(:),workn(:,:,:),workv(:,:,:)

!no_abirules
!Define G^2 based on G space metric gmet.
! gsq(g1,g2,g3)=g1*g1*gmet(1,1)+g2*g2*gmet(2,2)+g3*g3*gmet(3,3) &
! &       +two*(g1*g2*gmet(1,2)+g2*g3*gmet(2,3)+g3*g1*gmet(3,1))
! *************************************************************************

 DBG_ENTER("COLL")

!  Check optional arguments
 if (present(comm_fft)) then
   if ((.not.present(paral_kgb)).or.(.not.present(me_g0))) then
     MSG_BUG('Need paral_kgb and me_g0 with comm_fft !')
   end if
 end if

 if (present(gauss))then
   if (.not.present(optn2_in)) then   
     optn2 = 3
   else
     if(optn2_in/=3)then
       MSG_BUG('optn2 must be set to 3!')
     else
       optn2 = optn2_in
     end if
   end if     
 end if

 if (present(atmrhor1))then
   if(.not.present(optn_in))then
     optn = 1
   else
     optn = optn_in   
   end if
   if (.not.present(optn2_in)) then   
     MSG_BUG('rho1 calculation need optn2 !')
   else
     optn2 = optn2_in
   end if
 else
   optn  = 0
   optn2 = 0
 end if
 
 if (present(atmvlocr1))then
   if(.not.present(optv_in))then
     optv = 1
   else
     if(optv_in/=1)then
       MSG_BUG('optv_in must be set to 1!')
     else
       optv = optv_in
     end if
   end if
   if(.not.present(vspl))then
     MSG_BUG('vloc1 calculation need vspl!')
   end if
 else
   optv  = 0 
 end if
 
 if(ipert==natom+1.or.ipert==natom+2.or.ipert==natom+10.or.ipert==natom+11) then

!  (In case of d/dk or an electric/magnetic field)
   if (optn==1) then
     atmrhor1(1:cplex*nfft,1:ndir)=zero
     if (present(atmrhog1)) atmrhog1 = zero
   end if
   if (optv==1) then 
     atmvlocr1(1:cplex*nfft,1:ndir)=zero
     if (present(atmvlocg1)) atmvlocg1 = zero
   end if

 else

!   Useful quantities
   if (ipert/=natom+3.and.ipert/=natom+4) then
     iatom=ipert;iatm=atindx(iatom)
     itypat=typat(iatom)
   else
      !sum of all (strain pertubation)
     iatom  = 1
     iatm   = 1
     itypat = 1
   end if

   n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
   me_fft=ngfft(11)
   nproc_fft=ngfft(10)
   if (ndir==1) then
     jdir(1)=idir
   else
     do id=1,ndir
       jdir(id)=id
     end do
   end if

!   Get the distrib associated with this fft_grid
   if (present(distribfft)) then
     my_distribfft => distribfft
   else
     ABI_DATATYPE_ALLOCATE(my_distribfft,)
     call init_distribfft_seq(my_distribfft,'f',n2,n3,'fourdp')
   end if
   if (n2==my_distribfft%n2_coarse) then
     fftn2_distrib => my_distribfft%tab_fftdp2_distrib
   else if (n2 == my_distribfft%n2_fine) then
     fftn2_distrib => my_distribfft%tab_fftdp2dg_distrib
   else
     MSG_BUG("Unable to find an allocated distrib for this fft grid")
   end if

   qeq0=(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15)
   qeq05=(abs(abs(qphon(1))-half)<tol12.or. &
&   abs(abs(qphon(2))-half)<tol12.or. &
&   abs(abs(qphon(3))-half)<tol12)

   if (nproc_fft>1.and.qeq05) then
     MSG_ERROR('not compatible with FFT parallelism')
   end if

   if (optn2==3)then
     gauss1=gauss(1,itypat)
     gauss2=gauss(2,itypat)
     alf2pi2=(two_pi*gauss2)**2
   end if

   dq=(qgrid(mqgrid)-qgrid(1))/dble(mqgrid-1)
   dqm1=one/dq
   dqdiv6=dq/six
   dq2div6=dq**2/six
   cutoff=gsqcut*tolfix
   id1=n1/2+2;id2=n2/2+2;id3=n3/2+2
   ig1max=-1;ig2max=-1;ig3max=-1
   ig1min=n1;ig2min=n2;ig3min=n3

!   Determination of phase qxred*
   qxred2pi=two_pi*(qphon(1)*xred(1,iatom)+ &
&   qphon(2)*xred(2,iatom)+ &
&   qphon(3)*xred(3,iatom) )
   phqre=cos(qxred2pi)
   phqim=sin(qxred2pi)

!   Zero out temporary arrays
   if (optv==1) then
     ABI_ALLOCATE(workv,(2,nfft,ndir))
     workv(:,:,:)=zero
   end if
   if (optn==1) then
     ABI_ALLOCATE(workn,(2,nfft,ndir))
     workn(:,:,:)=zero
   end if

   if (ipert==natom+3.or.ipert==natom+4) then
     ntype = ntypat
     ia1=1
     type1 = 1
     type2 = ntype
     ABI_ALLOCATE(phre_igia,(natom))
     ABI_ALLOCATE(phim_igia,(natom))
   else
     type1 = itypat
     type2 = itypat
     ntype = 1
     ABI_ALLOCATE(phre_igia,(iatm:iatm))
     ABI_ALLOCATE(phim_igia,(iatm:iatm))
   end if


   do itypat=type1,type2
!    ia1,ia2 sets range of loop over atoms:
     if (ipert==natom+3.or.ipert==natom+4) then
       nattyp = count(typat(:)==itypat)
       ia2=ia1+nattyp-1
     else
       ia1 = iatm
       ia2 = iatm
     end if

     if (usepaw == 1) then
       tcorespl => pawtab(itypat)%tcorespl
       tvalespl => pawtab(itypat)%tvalespl
     else
       tcorespl => psps%nctab(itypat)%tcorespl
       tvalespl => psps%nctab(itypat)%tvalespl
     end if

     ii=0
     do i3=1,n3
       ig3=i3-(i3/id3)*n3-1
       gq3=dble(ig3)+qphon(3)
       gq(3)=gq3

       do i2=1,n2
         if (fftn2_distrib(i2)==me_fft) then
           ig2=i2-(i2/id2)*n2-1
           gq2=dble(ig2)+qphon(2)
           gq(2)=gq2
           
           do i1=1,n1
             ig1=i1-(i1/id1)*n1-1
             gq1=dble(ig1)+qphon(1)
             gq(1)=gq1

             ii=ii+1
!            gsquar=gsq(gq1,gq2,gq3)
             gsquar=gq1*gq1*gmet(1,1)+gq2*gq2*gmet(2,2)+gq3*gq3*gmet(3,3) &
&             +two*(gq1*gq2*gmet(1,2)+gq2*gq3*gmet(2,3)+gq3*gq1*gmet(3,1))

             
!             Skip G**2 outside cutoff:
             if (gsquar<=cutoff) then
               
!               Identify min/max indexes (to cancel unbalanced contributions later)
               if (qeq05) then
                 ig1max=max(ig1max,ig1);ig1min=min(ig1min,ig1)
                 ig2max=max(ig2max,ig2);ig2min=min(ig2min,ig2)
                 ig3max=max(ig3max,ig3);ig3min=min(ig3min,ig3)
               end if
               
               gmag=sqrt(gsquar)
               have_g0=(ig1==0.and.ig2==0.and.ig3==0.and.qeq0)

               jj=1+int(gmag*dqm1)
               diff=gmag-qgrid(jj)
               
!               Compute structure factor
               phre_igia(:) = zero
               phim_igia(:) = zero
               
               do ia=ia1,ia2
                 shift1=1+n1+(ia-1)*(2*n1+1)
                 shift2=1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1)
                 shift3=1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
                 ph1r=ph1d(re,ig1+shift1);ph1i=ph1d(im,ig1+shift1)
                 ph2r=ph1d(re,ig2+shift2);ph2i=ph1d(im,ig2+shift2)
                 ph3r=ph1d(re,ig3+shift3);ph3i=ph1d(im,ig3+shift3)
                 ph12r=ph1r*ph2r-ph1i*ph2i
                 ph12i=ph1r*ph2i+ph1i*ph2r
                 phre_igia(ia)=ph12r*ph3r-ph12i*ph3i
                 phim_igia(ia)=ph12r*ph3i+ph12i*ph3r
               end do
!              Compute V^AT(g+q) and/or n^AT(g+q) for given type of atom
!              Evaluate spline fit: p. 86 Numerical Recipes, Press et al;
!              Note the error in book for sign of "aa" term in derivative.
               if (optv==1.or.optn2/=3) then
                 bb = diff*dqm1
                 aa = one-bb
                 cc = aa*(aa**2-one)*dq2div6
                 dd = bb*(bb**2-one)*dq2div6
               end if

               if (optv==1) then
                 if (have_g0) then
                   v_at=zero
                   dv_at=zero
                 else
                   v_at=(aa*vspl(jj,1,itypat)+bb*vspl(jj+1,1,itypat)+&
&                   cc*vspl(jj,2,itypat)+dd*vspl(jj+1,2,itypat))/gsquar

!                   Also get (dV(q)/dq)/q:
!                   (note correction of Numerical Recipes sign error
!                   before (3._dp*aa**2-1._dp)
                   ee= vspl(jj+1,1,itypat)-vspl(jj,1,itypat)
                   ff=  (3._dp*bb**2-1._dp)*vspl(jj+1,2,itypat) &
&                   - (3._dp*aa**2-1._dp)*vspl(jj,2,itypat)
                   dv_at = ( ( ee*dqm1 + ff*dqdiv6 )/gmag&
&                   - 2.0_dp*v_at) / gsquar                   
                 end if
               end if ! end optv

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

!                Also get (dn^AT(q)/dq)/q:
                 if (have_g0) then
                   if (optn2==1) then
                     if (usepaw == 1) then
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
                 else
                   if (optn2==1) then
                     ee=tcorespl(jj+1,1)-tcorespl(jj,1)
                     ff=(3._dp*bb**2-1._dp)*tcorespl(jj+1,2) &
&                     -(3._dp*aa**2-1._dp)*tcorespl(jj,2)
                   else if (optn2==2) then
                     ee=tvalespl(jj+1,1)-tvalespl(jj,1)
                     ff=(3._dp*bb**2-1._dp)*tvalespl(jj+1,2) &
&                     -(3._dp*aa**2-1._dp)*tvalespl(jj,2)
                   else if (optn2==3) then
                     dn_at=-two*gauss1*alf2pi2*exp(-gsquar*alf2pi2)
                   else
                   end if
                   dn_at=(ee*dqm1+ff*dqdiv6)/gmag
                 end if
               end if ! end optn
               
               do id=1,ndir
                 sfr=zero;sfi=zero                
                 do ia=ia1,ia2
                   if (ipert==natom+3.or.ipert==natom+4) then
!                    sum[Exp(-i.2pi.g.xred)]  
                     sfr=sfr+phre_igia(ia)
                     sfi=sfi-phim_igia(ia)
                   else
!                    Exp(-i.2pi.g.xred)  * -i.2pi.(g+q)
                     sfr=-(two_pi*gq(jdir(id))*phim_igia(ia))
                     sfi=-(two_pi*gq(jdir(id))*phre_igia(ia))
                   end if
                 end do

                 if (ipert/=natom+3.and.ipert/=natom+4) then
!                  Phonons case

!                  Exp(-i.2pi.q.xred)       => -i.2pi.(g+q).Exp(-i.2pi.(g+q).xred)
                   sfqr= sfr*phqre+sfi*phqim
                   sfqi=-sfr*phqim+sfi*phqre

                   if (optv == 1) then
                     workv(re,ii,id) = workv(re,ii,id) + sfqr*v_at
                     workv(im,ii,id) = workv(im,ii,id) + sfqi*v_at
                   end if
                   if (optn == 1) then
                     workn(re,ii,id) = workn(re,ii,id) + sfqr*n_at
                     workn(im,ii,id) = workn(im,ii,id) + sfqi*n_at
                   end if

                 else
!                  Strain case
                   
!                  Compute G in cartesian coordinates
                   gcart(1)=gprimd(1,1)*dble(ig1)+gprimd(1,2)*dble(ig2)+&
&                   gprimd(1,3)*dble(ig3)
                   gcart(2)=gprimd(2,1)*dble(ig1)+gprimd(2,2)*dble(ig2)+&
&                   gprimd(2,3)*dble(ig3)
                   gcart(3)=gprimd(3,1)*dble(ig1)+gprimd(3,2)*dble(ig2)+&
&                   gprimd(3,3)*dble(ig3)

!                  Accumulate -dV^AT/dG*rho(G)*SF(G)*Gi.Gj/G
!                  or -dn^AT/dG*V(G)*SF(G)*Gi.Gj/G              
                   if (optv==1) then
                     if(jdir(id)<=3) then
                       term_v = dv_at*gcart(eps1(jdir(id)))*gcart(eps2(jdir(id))) + v_at
                     else
                       term_v = dv_at*gcart(eps1(jdir(id)))*gcart(eps2(jdir(id)))
                     end if
                     workv(re,ii,id) = workv(re,ii,id) - (sfr*term_v)
                     workv(im,ii,id) = workv(im,ii,id) - (sfi*term_v)
                   end if

                   if (optn==1) then
                     if(jdir(id)<=3) then
                       term_n = dn_at*gcart(eps1(jdir(id)))*gcart(eps2(jdir(id))) + n_at
                     else
                       term_n = dn_at*gcart(eps1(jdir(id)))*gcart(eps2(jdir(id)))
                     end if

                     workn(re,ii,id) = workn(re,ii,id) - (sfr*term_n)
                     workn(im,ii,id) = workn(im,ii,id) - (sfi*term_n)
                   end if

                 end if
!                End loop on ndir

               end do
!              End skip G**2 outside cutoff
             end if
!            End loop on n1, n2, n3
           end do
         end if ! this plane is selected
       end do
     end do
     ia1=ia2+1
   end do ! end loop itype

   ABI_DEALLOCATE(phre_igia)
   ABI_DEALLOCATE(phim_igia)

   if(ipert==natom+3.or.ipert==natom+4) then
!    Set Vloc(G=0)=0:
     if (optv==1) then
       workv(re,1,:)=zero
       workv(im,1,:)=zero
     end if
   end if
   
!  Identify unbalanced g-vectors
   if (qeq05) then  !This doesn't work in parallel
     ig1=-1;if (mod(n1,2)==0) ig1=1+n1/2
     ig2=-1;if (mod(n2,2)==0) ig2=1+n2/2
     ig3=-1;if (mod(n3,2)==0) ig3=1+n3/2
     if (abs(abs(qphon(1))-half)<tol12) then
       if (abs(ig1min)<abs(ig1max)) ig1=abs(ig1max)
       if (abs(ig1min)>abs(ig1max)) ig1=n1-abs(ig1min)
     end if
     if (abs(abs(qphon(2))-half)<tol12) then
       if (abs(ig2min)<abs(ig2max)) ig2=abs(ig2max)
       if (abs(ig2min)>abs(ig2max)) ig2=n2-abs(ig2min)
     end if
     if (abs(abs(qphon(3))-half)<tol12) then
       if (abs(ig3min)<abs(ig3max)) ig3=abs(ig3max)
       if (abs(ig3min)>abs(ig3max)) ig3=n3-abs(ig3min)
     end if
   end if


!  Get 1st-order potential/density back to real space
!  Non-symetrized non-zero elements have to be nullified
!  Divide by unit cell volume

   if (optv==1.or.optn==1) then

     xnorm=one/ucvol
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
     
     if (optv==1) then
       do id=1,ndir
!        Eliminate unbalanced g-vectors
         if (qeq0) then       !q=0
           call zerosym(workv(:,:,id),2,n1,n2,n3,comm_fft=my_comm_fft,distribfft=my_distribfft)
         else if (qeq05) then !q=1/2; this doesn't work in parallel
           call zerosym(workv(:,:,id),2,n1,n2,n3,ig1=ig1,ig2=ig2,ig3=ig3)
         end if
         call fourdp(cplex,workv(:,:,id),atmvlocr1(:,id),1,mpi_enreg_fft,nfft,ngfft,paral_kgb_fft,0)  
         atmvlocr1(:,id)=atmvlocr1(:,id)*xnorm
       end do

       !if (present(atmvlocg1)) atmvlocg1 = workv
       ABI_DEALLOCATE(workv)
     end if
     
     if (optn==1) then
       do id=1,ndir
!        Eliminate unbalanced g-vectors
         if (qeq0) then       !q=0
           call zerosym(workn(:,:,id),2,n1,n2,n3,comm_fft=my_comm_fft,distribfft=my_distribfft)
         else if (qeq05) then !q=1/2; this doesn't work in parallel
           call zerosym(workn(:,:,id),2,n1,n2,n3,ig1=ig1,ig2=ig2,ig3=ig3)
         end if
         call fourdp(cplex,workn(:,:,id),atmrhor1(:,id),1,mpi_enreg_fft,nfft,ngfft,paral_kgb_fft,0)
         atmrhor1(:,id)=atmrhor1(:,id)*xnorm
       end do
       !if (present(atmrhog1)) atmrhog1 = workn
       ABI_DEALLOCATE(workn)
     end if

!    Destroy fake mpi_enreg
     call unset_mpi_enreg_fft(mpi_enreg_fft)
   end if
   
   if (.not.present(distribfft)) then
     call destroy_distribfft(my_distribfft)
     ABI_DATATYPE_DEALLOCATE(my_distribfft)
   end if

!  End the condition of non-electric-field
 end if

 DBG_EXIT("COLL")
 
end subroutine dfpt_atm2fft
!!***
