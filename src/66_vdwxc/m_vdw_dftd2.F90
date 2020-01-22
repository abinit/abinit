!!****m* ABINIT/m_vdw_dftd2
!! NAME
!!  m_vdw_dftd2
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MT)
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

module m_vdw_dftd2

 use defs_basis
 use m_abicore
 use m_errors
 use m_atomdata

 use m_geometry,         only : metric

 implicit none

 private
!!***

 public :: vdw_dftd2
!!***

contains
!!***

!!****f* ABINIT/vdw_dftd2
!!
!! NAME
!! vdw_dftd2
!!
!! FUNCTION
!! Compute energy and derivatives with respect to dimensionless
!! reduced atom coordinates due to Van der Waals interaction.
!! The formalism here follows the DFT-D2 approach of Grimme
!! which consists in adding a semi-empirical dispersion potential
!! (pair-wise force field) to the conventional Kohn-Sham DFT energy.
!!
!! INPUTS
!!  ixc=choice of exchange-correlation functional
!!  natom=number of atoms
!!  ntypat=number of atom types
!!  prtvol=printing volume (if >0, print computation parameters)
!!  typat(natom)=type integer for each atom in cell
!!  rprimd(3,3)=real space primitive translations
!!  vdw_tol=tolerance use to converge the potential (a pair of atoms is included
!!          in potential if its contribution is larger than vdw_tol)
!!          vdw_tol<0 takes default value (10^-10)
!!  xred(3,natom)=reduced atomic coordinates
!!  znucl(ntypat)=atomic number of atom type
!!  === optional inputs ===
!!  [qphon(3)]=wavevector of the phonon;
!!             used only for dynamical matrix computation
!!
!! OUTPUT
!!  e_vdw_dftd2=contribution to energy from DFT-D2 dispersion potential
!!  === optional outputs ===
!!  [dyn_vdw_dftd2(2,3,natom,3,natom)]=contribution to dynamical matrix from DFT-D2 dispersion potential
!!  [elt_vdw_dftd2(6+3*natom,6)]=contribution to elastic tensor and internal strains from DFT-D2 disp. pot.
!!  [fred_vdw_dftd2(3,natom)]=contribution to forces from DFT-D2 dispersion potential
!!  [str_vdw_dftd2(6)]=contribution to stress tensor from DFT-D2 dispersion potential
!!
!! NOTES
!!  Ref.: S. Grimme, Semiempirical GGA-type density functional
!!        constructed with a long-range dispersion correction,
!!        J. Comp. Chem. 27, 1787 (2006) [[cite:Grimme2006]]
!!
!! PARENTS
!!      respfn,setvtr,stress
!!
!! CHILDREN
!!
!! SOURCE

subroutine vdw_dftd2(e_vdw_dftd2,ixc,natom,ntypat,prtvol,typat,rprimd,vdw_tol,xred,znucl,&
&          dyn_vdw_dftd2,elt_vdw_dftd2,fred_vdw_dftd2,str_vdw_dftd2,qphon) ! Optionals

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,natom,ntypat,prtvol
 real(dp),intent(in) :: vdw_tol
 real(dp),intent(out) :: e_vdw_dftd2
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom),znucl(ntypat)
 real(dp),intent(in),optional :: qphon(3)
 real(dp),intent(out),optional :: dyn_vdw_dftd2(2,3,natom,3,natom)
 real(dp),intent(out),optional :: elt_vdw_dftd2(6+3*natom,6)
 real(dp),intent(out),optional :: fred_vdw_dftd2(3,natom)
 real(dp),intent(out),optional :: str_vdw_dftd2(6)

!Local variables-------------------------------
!scalars
 integer,parameter :: vdw_nspecies=55
 integer :: ia,ia1,ii,is1,is2,is3,itypat,ja,ja1,jj,jtypat,kk,ll,mu,npairs,nshell,nu
 logical :: need_dynmat,need_elast,need_forces,need_intstr,need_stress
 logical :: need_gradient,need_gradient2,newshell,qeq0=.true.
 real(dp),parameter :: e_conv=(10/Bohr_Ang)**6/Ha_J/Avogadro ! 1 J.nm^6.mol^-1 in Ha.Bohr^6
 real(dp),parameter :: vdw_d=20._dp,vdw_tol_default=tol10
 real(dp),parameter :: vdw_s_pbe=0.75_dp, vdw_s_blyp=1.2_dp, vdw_s_b3lyp=1.05_dp
 real(dp),parameter :: vdw_s_bp86=1.05_dp, vdw_s_tpss=1.0_dp, vdw_s_b97d=1.25_dp
 real(dp) :: c6,c6r6,ex,fr,gr,gr2,grad,grad2,ph,ph1r,ph1i
 real(dp) :: r0,r1,r2,r3,rcut,rcut2,rsq,rr,sfact,ucvol,vdw_s
 character(len=500) :: msg
 type(atomdata_t) :: atom
!arrays
 integer,allocatable :: ivdw(:)
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 real(dp) :: gmet(3,3),gprimd(3,3),mat(3,3),rcart(3),rmet(3,3),vec(3)
 real(dp),allocatable :: vdw_c6(:,:),vdw_r0(:,:),xred01(:,:)
 real(dp),parameter :: vdw_c6_dftd2(vdw_nspecies)= &
&      (/ 0.14, 0.08, 1.61, 1.61, 3.13, 1.75, 1.23, 0.70, 0.75, 0.63,&
&         5.71, 5.71,10.79, 9.23, 7.84, 5.57, 5.07, 4.61,10.80,10.80,&
&        10.80,10.80,10.80,10.80,10.80,10.80,10.80,10.80,10.80,10.80,&
&        16.99,17.10,16.37,12.64,12.47,12.01,24.67,24.67,24.67,24.67,&
&        24.67,24.67,24.67,24.67,24.67,24.67,24.67,24.67,37.32,38.71,&
&        38.44,31.74,31.50,29.99, 0.00/)
 real(dp),parameter :: vdw_r0_dftd2(vdw_nspecies)= &
&      (/1.001,1.012,0.825,1.408,1.485,1.452,1.397,1.342,1.287,1.243,&
&        1.144,1.364,1.639,1.716,1.705,1.683,1.639,1.595,1.485,1.474,&
&        1.562,1.562,1.562,1.562,1.562,1.562,1.562,1.562,1.562,1.562,&
&        1.650,1.727,1.760,1.771,1.749,1.727,1.628,1.606,1.639,1.639,&
&        1.639,1.639,1.639,1.639,1.639,1.639,1.639,1.639,1.672,1.804,&
&        1.881,1.892,1.892,1.881,1.000/)
 character(len=2),parameter :: vdw_symb(vdw_nspecies)= &
&      (/' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',&
&        'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca',&
&        'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',&
&        'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr',&
&        'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',&
&        'Sb','Te',' I','Xe','no'/)

! *************************************************************************

 DBG_ENTER("COLL")

!Extract options
 need_forces=present(fred_vdw_dftd2)
 need_stress=present(str_vdw_dftd2)
 need_dynmat=present(dyn_vdw_dftd2)
 need_elast=present(elt_vdw_dftd2)
 need_intstr=present(elt_vdw_dftd2)
 need_gradient=(need_forces.or.need_stress)
 need_gradient2=(need_dynmat.or.need_elast.or.need_intstr)
 if (need_dynmat) then
   if (.not.present(qphon)) then
     msg='Dynamical matrix required without a q-vector'
     MSG_BUG(msg)
   end if
   qeq0=(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15)
 end if

!Identify type(s) of atoms
 ABI_ALLOCATE(ivdw,(ntypat))
 do itypat=1,ntypat
   ivdw(itypat)=-1;jtypat=0
   call atomdata_from_znucl(atom,znucl(itypat))
   do while ((ivdw(itypat)<0).and.(jtypat<vdw_nspecies))
     jtypat=jtypat+1;if (vdw_symb(jtypat)==atom%symbol) ivdw(itypat)=jtypat
   end do
   if (ivdw(itypat)<0) then
     write(msg,'(3a)') &
&     'Van der Waals DFT-D2 correction not available for atom type: ',atom%symbol,' !'
     MSG_ERROR(msg)
   end if
 end do

!Select DFT-D2 VdW parameters according to system data
 vdw_s=e_conv
 if (ixc==11.or.ixc==-101130.or.ixc==-130101) then
   vdw_s=vdw_s*vdw_s_pbe
 else if (ixc==18.or.ixc==-106131.or.ixc==-131106) then
   vdw_s=vdw_s*vdw_s_blyp
 else if (ixc==19.or.ixc==-106132.or.ixc==-132106) then
   vdw_s=vdw_s*vdw_s_bp86
 else if (ixc==-202231.or.ixc==-231202) then
   vdw_s=vdw_s*vdw_s_tpss
 else
   write(msg,'(a,i8,a)')'  Van der Waals DFT-D2 correction not compatible with ixc=',ixc,' !'
   MSG_ERROR(msg)
 end if
 ABI_ALLOCATE(vdw_c6,(ntypat,ntypat))
 ABI_ALLOCATE(vdw_r0,(ntypat,ntypat))
 do itypat=1,ntypat
   do jtypat=1,ntypat
     vdw_c6(itypat,jtypat)=sqrt(vdw_c6_dftd2(ivdw(itypat))*vdw_c6_dftd2(ivdw(jtypat)))
     vdw_r0(itypat,jtypat)=(vdw_r0_dftd2(ivdw(itypat))+vdw_r0_dftd2(ivdw(jtypat)))/Bohr_Ang
   end do
 end do

!Computation of cut-off radius according to tolerance
!We take: r_cut=(s6*max(C6)/tol)**(1/6) and rcut<=75 bohr
 if (vdw_tol<zero) then
   rcut=(vdw_s/vdw_tol_default*maxval(vdw_c6))**sixth
 else
   rcut=(vdw_s/vdw_tol*maxval(vdw_c6))**sixth
 end if
!rcut=min(rcut,100._dp)
 rcut2=rcut*rcut

!Retrieve cell geometry data
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Map reduced coordinates into [0,1]
 ABI_ALLOCATE(xred01,(3,natom))
 do ia=1,natom
   xred01(1,ia)=xred(1,ia)-aint(xred(1,ia))+half-sign(half,xred(1,ia))
   xred01(2,ia)=xred(2,ia)-aint(xred(2,ia))+half-sign(half,xred(2,ia))
   xred01(3,ia)=xred(3,ia)-aint(xred(3,ia))+half-sign(half,xred(3,ia))
 end do

!Set accumulated quantities to zero
 npairs=0
 e_vdw_dftd2=zero
 if (need_forces) fred_vdw_dftd2=zero
 if (need_stress) str_vdw_dftd2=zero
 if (need_dynmat) dyn_vdw_dftd2=zero
 if (need_elast)  elt_vdw_dftd2(1:6,1:6)=zero
 if (need_intstr) elt_vdw_dftd2(7:6+3*natom,1:6)=zero

!Loop over shells of cell replicas
 nshell=0
 do
   newshell=.false.;nshell=nshell+1

!  Loop over cell replicas in the shell
!  ns1=1+int(rcut*sqrt(SUM(gprimd(:,1)**2))
!  ns2=1+int(rcut*sqrt(SUM(gprimd(:,2)**2))
!  ns3=1+int(rcut*sqrt(SUM(gprimd(:,3)**2))
   do is3=-nshell,nshell
     do is2=-nshell,nshell
       do is1=-nshell,nshell
         if (nshell==1.or. &
&         abs(is3)==nshell.or.abs(is2)==nshell.or.abs(is1)==nshell) then

!          Phase for dynamical matrix
           if (need_dynmat) then
             ph1r=one;ph1i=zero  !ph1=exp(-iqL)
             if (.not.qeq0) then
               ph=-two_pi*(qphon(1)*is1+qphon(2)*is2+qphon(3)*is3)
               ph1r=cos(ph);ph1i=sin(ph)
             end if
           end if

!          Loops over atoms a and b
           do ia=1,natom
             itypat=typat(ia)
             do ja=1,ia
               jtypat=typat(ja)
               r1=xred01(1,ia)-xred01(1,ja)-dble(is1)
               r2=xred01(2,ia)-xred01(2,ja)-dble(is2)
               r3=xred01(3,ia)-xred01(3,ja)-dble(is3)
               rsq=rmet(1,1)*r1*r1+rmet(2,2)*r2*r2+rmet(3,3)*r3*r3 &
&               +two*(rmet(2,1)*r2*r1+rmet(3,2)*r3*r2+rmet(3,1)*r1*r3)

!              Select atomic pairs (a,b) and avoid atom_a=atom_b
               if (rsq>=tol16.and.rsq<=rcut2) then

!                Data for the selected pair
                 npairs=npairs+1;newshell=.true.
                 sfact=vdw_s;if (ia==ja) sfact=half*sfact
                 rr=sqrt(rsq)
                 c6=vdw_c6(itypat,jtypat)
                 r0=vdw_r0(itypat,jtypat)

!                Computation of pair-wise potential
                 ex=exp(-vdw_d*(rr/r0-one))
                 fr=one/(one+ex)
                 c6r6=c6/rr**6

!                Contribution to energy
                 e_vdw_dftd2=e_vdw_dftd2-sfact*fr*c6r6

                 if (need_gradient.or.need_gradient2) then
                   gr=(vdw_d/r0)*(fr**2)*ex
                   grad=-sfact*(gr-six*fr/rr)*c6r6/rr
                   rcart(1)=rprimd(1,1)*r1+rprimd(1,2)*r2+rprimd(1,3)*r3
                   rcart(2)=rprimd(2,1)*r1+rprimd(2,2)*r2+rprimd(2,3)*r3
                   rcart(3)=rprimd(3,1)*r1+rprimd(3,2)*r2+rprimd(3,3)*r3

!                  Contribution to forces
                   if (need_forces.and.ia/=ja) then
                     vec(1:3)=grad*rcart(1:3)
                     fred_vdw_dftd2(1:3,ia)=fred_vdw_dftd2(1:3,ia)+vec(1:3)
                     fred_vdw_dftd2(1:3,ja)=fred_vdw_dftd2(1:3,ja)-vec(1:3)
                   end if

!                  Contribution to stress tensor
                   if (need_stress) then
                     do mu=1,6
                       ii=alpha(mu);jj=beta(mu)
                       str_vdw_dftd2(mu)=str_vdw_dftd2(mu)+grad*rcart(ii)*rcart(jj)
                     end do
                   end if

                   if (need_gradient2) then
                     gr2=(vdw_d/r0)*gr*(2*fr*ex-one)
                     grad2=-sfact*(gr2-13._dp*gr/rr+48._dp*fr/rr**2)*c6r6/rr**2

!                    Contribution to dynamical matrix (phase factors are subtle!)
                     if (need_dynmat) then
                       mat(1:3,1)=grad2*rcart(1:3)*rcart(1) ; mat(1,1)=mat(1,1)+grad
                       mat(1:3,2)=grad2*rcart(1:3)*rcart(2) ; mat(2,2)=mat(2,2)+grad
                       mat(1:3,3)=grad2*rcart(1:3)*rcart(3) ; mat(3,3)=mat(3,3)+grad
                       if (ia/=ja) then
                         do ii=1,3
                           dyn_vdw_dftd2(1,1:3,ia,ii,ia)=dyn_vdw_dftd2(1,1:3,ia,ii,ia)+mat(1:3,ii)
                           dyn_vdw_dftd2(1,1:3,ja,ii,ja)=dyn_vdw_dftd2(1,1:3,ja,ii,ja)+mat(1:3,ii)
                           dyn_vdw_dftd2(1,1:3,ia,ii,ja)=dyn_vdw_dftd2(1,1:3,ia,ii,ja)-mat(1:3,ii)*ph1r
                           dyn_vdw_dftd2(2,1:3,ia,ii,ja)=dyn_vdw_dftd2(2,1:3,ia,ii,ja)-mat(1:3,ii)*ph1i
                           dyn_vdw_dftd2(1,1:3,ja,ii,ia)=dyn_vdw_dftd2(1,1:3,ja,ii,ia)-mat(1:3,ii)*ph1r
                           dyn_vdw_dftd2(2,1:3,ja,ii,ia)=dyn_vdw_dftd2(2,1:3,ja,ii,ia)+mat(1:3,ii)*ph1i
                         end do
                       else if (.not.qeq0) then
                         do ii=1,3
                           dyn_vdw_dftd2(1,1:3,ia,ii,ia)=dyn_vdw_dftd2(1,1:3,ia,ii,ia) &
&                           +two*mat(1:3,ii)*(one-ph1r)
                         end do
                       end if
                     end if

!                    Contribution to elastic tensor
                     if (need_elast) then
                       do mu=1,6
                         ii=alpha(mu);jj=beta(mu)
                         do nu=1,6
                           kk=alpha(nu);ll=beta(nu)
                           elt_vdw_dftd2(mu,nu)=elt_vdw_dftd2(mu,nu) &
&                           +grad2*rcart(ii)*rcart(jj)*rcart(kk)*rcart(ll)
                           if (ii==kk) elt_vdw_dftd2(mu,nu)=elt_vdw_dftd2(mu,nu) &
&                           +half*grad*rcart(jj)*rcart(ll)
                           if (ii==ll) elt_vdw_dftd2(mu,nu)=elt_vdw_dftd2(mu,nu) &
&                           +half*grad*rcart(jj)*rcart(kk)
                           if (jj==kk) elt_vdw_dftd2(mu,nu)=elt_vdw_dftd2(mu,nu) &
&                           +half*grad*rcart(ii)*rcart(ll)
                           if (jj==ll) elt_vdw_dftd2(mu,nu)=elt_vdw_dftd2(mu,nu) &
&                           +half*grad*rcart(ii)*rcart(kk)
                         end do
                       end do
                     end if

!                    Contribution to internal strains
                     if (need_intstr.and.ia/=ja) then
                       ia1=6+3*(ia-1);ja1=6+3*(ja-1)
                       do mu=1,6
                         ii=alpha(mu);jj=beta(mu)
                         vec(1:3)=grad2*rcart(ii)*rcart(jj)*rcart(1:3)
                         vec(ii)=vec(ii)+half*grad*rcart(jj)
                         vec(jj)=vec(jj)+half*grad*rcart(ii)
                         elt_vdw_dftd2(ia1+1:ia1+3,mu)=elt_vdw_dftd2(ia1+1:ia1+3,mu)+vec(1:3)
                         elt_vdw_dftd2(ja1+1:ja1+3,mu)=elt_vdw_dftd2(ja1+1:ja1+3,mu)-vec(1:3)
                       end do
                     end if

                   end if ! Computation of 2nd gradient
                 end if ! Computation of gradient
               end if   ! Pairs selection
             end do     ! Loop over atom b
           end do       ! Loop over atom a
         end if         ! Triple loop over cell replicas in shell
       end do
     end do
   end do
   if(.not.newshell) exit ! Check if new shell must be calculated
 end do ! Loop over shells

!Gradients: convert them from cartesian to reduced coordinates
 if (need_forces) then
   do ia=1,natom
     call grad_cart2red(fred_vdw_dftd2(:,ia))
   end do
 end if
 if (need_dynmat) then
   do ja=1,natom
     do ia=1,natom
       do kk=1,merge(2,1,qeq0)
         do ii=1,3
           vec(1:3)=dyn_vdw_dftd2(kk,1:3,ia,ii,ja)
           call grad_cart2red(vec)
           dyn_vdw_dftd2(kk,1:3,ia,ii,ja)=vec(1:3)
         end do
         do ii=1,3
           vec(1:3)=dyn_vdw_dftd2(kk,ii,ia,1:3,ja)
           call grad_cart2red(vec)
           dyn_vdw_dftd2(kk,ii,ia,1:3,ja)=vec(1:3)
         end do
       end do
     end do
   end do
 end if
 if (need_intstr) then
   do mu=1,6
     ia1=6
     do ia=1,natom
       call grad_cart2red(elt_vdw_dftd2(ia1+1:ia1+3,mu))
       ia1=ia1+3
     end do
   end do
 end if

!DEBUG
!write(77,*) "---------------"
!write(77,*) "E=",e_vdw_dftd2
!if (need_forces) then
! do ia=1,natom
!  write(77,*) "F=",ia,fred_vdw_dftd2(:,ia)
! end do
!end if
!if (need_stress) write(77,*) "S=",str_vdw_dftd2(:)
!if (need_dynmat) then
! do ia=1,natom
!  do ii=1,3
!   do ja=1,natom
!    write(77,*) "D=",ia,ii,ja,dyn_vdw_dftd2(:,:,ja,ii,ia)
!   end do
!  end do
! end do
!end if
!if (need_elast) then
! do ii=1,6
!  write(77,*) "e=",ii,elt_vdw_dftd2(1:6,ii)
! end do
!end if
!if (need_intstr) then
! do ii=1,6
!  do ia=1,natom
!   write(77,*) "I=",ii,ia,elt_vdw_dftd2(7+3*(ia-1):9+3*(ia-1),ii)
!  end do
! end do
!end if
!flush(77)
!DEBUG

!Stress tensor: divide by volume
 if (need_stress) str_vdw_dftd2=str_vdw_dftd2/ucvol

!Printing
 if (prtvol>0) then
   write(msg,'(10a)') ch10,&
&   '  --------------------------------------------------------------',ch10,&
&   '  Van der Waals DFT-D2 semi-empirical dispersion potential added',ch10,&
&   '      with following parameters:',ch10,&
&   '      Specie  C6 (J.nm^6.mol^-1)  R0 (Ang)',ch10,&
&   '      ------------------------------------'
   call wrtout(std_out,msg,'COLL')
   do itypat=1,ntypat
     write(msg,'(9X,a2,11X,f5.2,8X,f6.3)') &
&     vdw_symb(ivdw(itypat)),vdw_c6_dftd2(ivdw(itypat)),vdw_r0_dftd2(ivdw(itypat))
     call wrtout(std_out,msg,'COLL')
   end do
   write(msg,'(2a,f6.2,2a,f6.2,2a,f6.2,a)') ch10,&
&   '      Scaling factor   = ',vdw_s/e_conv,ch10,&
&   '      Damping parameter= ',vdw_d,ch10,&
&   '      Cut-off radius   = ',rcut,' bohr'
   call wrtout(std_out,msg,'COLL')
   write(msg,'(2a,i8,2a,es14.5,4a)') ch10,&
&   '      Number of pairs contributing = ',npairs,ch10,&
&   '      DFT-D2 energy contribution   = ',e_vdw_dftd2,' Ha',ch10,&
&   '  --------------------------------------------------------------',ch10
   call wrtout(std_out,msg,'COLL')
 end if

 ABI_DEALLOCATE(ivdw)
 ABI_DEALLOCATE(vdw_c6)
 ABI_DEALLOCATE(vdw_r0)
 ABI_DEALLOCATE(xred01)

 DBG_EXIT("COLL")

 contains
!!***

!!****f* vdw_dftd2/grad_cart2red
!!
!! NAME
!! grad_cart2red
!!
!! FUNCTION
!! Convert gradients from cartesian to reduced coordinates
!!
!! PARENTS
!!      vdw_dftd2
!!
!! CHILDREN
!!
!! SOURCE

subroutine grad_cart2red(grad)

implicit none

!Arguments ------------------------------------
 real(dp),intent(inout) :: grad(3)
!Local variables-------------------------------
 real(dp) :: tmp(3)

! *********************************************************************

   tmp(1)=rprimd(1,1)*grad(1)+rprimd(2,1)*grad(2)+rprimd(3,1)*grad(3)
   tmp(2)=rprimd(1,2)*grad(1)+rprimd(2,2)*grad(2)+rprimd(3,2)*grad(3)
   tmp(3)=rprimd(1,3)*grad(1)+rprimd(2,3)*grad(2)+rprimd(3,3)*grad(3)
   grad(1:3)=tmp(1:3)

 end subroutine grad_cart2red
!!***

end subroutine vdw_dftd2
!!***

end module m_vdw_dftd2
!!***
