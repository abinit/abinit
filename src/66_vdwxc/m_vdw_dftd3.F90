!!****m* ABINIT/m_vdw_dftd3
!! NAME
!!  m_vdw_dftd3
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2015-2024 ABINIT group (BVT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_vdw_dftd3

 use defs_basis
 use m_abicore
 use m_errors
 use m_atomdata

 use m_special_funcs,  only : abi_derfc
 use m_geometry,       only : metric
 use m_vdw_dftd3_data, only : vdw_dftd3_data

 implicit none

 private
!!***

 public :: vdw_dftd3
!!***

contains
!!***

!!****f* ABINIT/vdw_dftd3
!!
!! NAME
!! vdw_dftd3
!!
!! FUNCTION
!! Compute energy, forces, stress, interatomic force constant and elastic
!! contribution due to dispersion interaction as formulated by Grimme in
!! the DFT-D3 approach. The last cited adds a dispersion potential
!! (pair-wise force field, rij^6 and rij^8) to Kohn-Sham DFT energy.
!! It is also possible to include a three-body term and molecular
!! dispersion (using vdw_tol_3bt>0).
!! DFT-D3(Becke and Johnson), another formulation which avoids the use of a damping
!! function to remove the undesired short-range behaviour
!!  is also activable using vdw_xc=7
!!
!! INPUTS
!!  ixc=choice of exchange-correlation functional
!!  natom=number of atoms
!!  ntypat=number of atom types
!!  prtvol=printing volume (if >0, print computation parameters)
!!  typat(natom)=type integer for each atom in cell
!!  vdw_xc= select van-der-Waals correction
!!             if =6: DFT-D3 as in Grimme, J. Chem. Phys. 132, 154104 (2010) [[cite:Grimme2010]]
!!             if =7: DFT-D3(BJ) as in Grimme, Comput. Chem. 32, 1456 (2011) [[cite:Grimme2011]]
!!                    Only the use of R0 = a1 C8/C6 + a2 is available here
!!
!!  vdw_tol=tolerance use to converge the pair-wise potential
!!          (a pair of atoms is included in potential if its contribution
!!          is larger than vdw_tol) vdw_tol<0 takes default value (10^-10)
!!  vdw_tol_3bt= tolerance use to converge three body terms (only for vdw_xc=6)
!!               a triplet of atom contributes to the correction if its
!!               contribution is larger than vdw_tol_3bt
!!  xred(3,natom)=reduced atomic coordinates
!!  znucl(ntypat)=atomic number of atom type
!!  === optional input ===
!!  [qphon(3)]= reduced q-vector along which is computed the DFT-D3 contribution
!!  to the IFCs in reciprocal space
!!
!! OUTPUT
!!  e_vdw_dftd3=contribution to energy from DFT-D3 dispersion potential
!!  === optional outputs ===
!!  [elt_vdw_dftd3(6+3*natom,6)]= contribution to elastic constant and
!!  internal strains from DFT-D3 dispersion potential
!!  [gred_vdw_dftd3(3,natom)]=contribution to gradient w.r.to atomic displ.
!!  from DFT-D3 dispersion potential
!!  [str_vdw_dftd3(6)]=contribution to stress tensor from DFT-D3 dispersion potential
!!  [dyn_vdw_dftd3(2,3,natom,3,natom)]= contribution to the interatomic force
!!  constants (in reciprocal space) at given input q-vector
!!  from DFT-D3 dispersion potential
!!
!! NOTES
!!  Ref.:
!!  DFT-D3: S. Grimme, J. Antony, S. Ehrlich, and H. Krieg
!!  A consistent and accurate ab initio parametrization of density functional
!!  dispersion correction (DFT-D) for the 94 elements H-Pu
!!  J. Chem. Phys. 132, 154104 (2010) [[cite:Grimme2010]]
!!  DFT-D3(BJ) S. Grimme, S. Ehrlich and L. Goerigk
!!  Effect of the damping function in dispersion corrected density functional theory
!!  Comput. Chem. 32, 1456 (2011) [[cite:Grimme2011]]
!!
!! SOURCE

subroutine vdw_dftd3(e_vdw_dftd3,ixc,natom,ntypat,prtvol,typat,rprimd,vdw_xc,&
&          vdw_tol,vdw_tol_3bt,xred,znucl,dyn_vdw_dftd3,elt_vdw_dftd3,&
&          gred_vdw_dftd3,str_vdw_dftd3,qphon)

implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,natom,ntypat,prtvol,vdw_xc
 real(dp),intent(in) :: vdw_tol,vdw_tol_3bt
 real(dp),intent(out) :: e_vdw_dftd3
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) ::  rprimd(3,3),xred(3,natom),znucl(ntypat)
 real(dp),intent(in),optional :: qphon(3)
 real(dp),intent(out),optional :: dyn_vdw_dftd3(2,3,natom,3,natom)
 real(dp),intent(out),optional :: elt_vdw_dftd3(6+3*natom,6)
 real(dp),intent(out),optional :: gred_vdw_dftd3(3,natom)
 real(dp),intent(out),optional :: str_vdw_dftd3(6)

!Local variables-------------------------------
!scalars
! The maximal number of reference systems for c6 is 5 (for C)
 integer,parameter :: vdw_nspecies=94
 integer:: alpha,beta,ia,ii,indi,indj,index_ia,index_ja,index_ka
 integer :: is1,is2,is3,itypat,ja,jj,js1,js2,js3
 integer :: jtypat,ka,kk,ktypat,la,ll,ierr
 integer :: nline,npairs,nshell
 integer :: refi,refj,refmax
 logical :: bol_3bt,found
 logical :: need_dynmat,need_elast,need_forces,need_grad,need_hess,need_stress,newshell
 real(dp),parameter :: alpha6=14.0_dp, alpha8=16.0_dp
 real(dp),parameter:: k1=16.0_dp, k2=15.0_dp, k3=4.0_dp

! s6 parameters (BJ case)
 real(dp),parameter :: vdwbj_s6_b2gpplyp=0.560_dp, vdwbj_s6_ptpss=0.750_dp
 real(dp),parameter :: vdwbj_s6_b2plyp=0.640_dp,   vdwbj_s6_dsdblyp=0.500_dp
 real(dp),parameter :: vdwbj_s6_pwpb95=0.820_dp
! s8 parameters (BJ case)
 real(dp),parameter :: vdwbj_s8_b1b95=1.4507_dp,    vdwbj_s8_b2gpplyp=0.2597_dp
 real(dp),parameter :: vdwbj_s8_b3pw91=2.8524_dp,   vdwbj_s8_bhlyp=1.0354_dp
 real(dp),parameter :: vdwbj_s8_bmk=2.0860_dp,      vdwbj_s8_bop=3.295_dp
 real(dp),parameter :: vdwbj_s8_bpbe=4.0728_dp,     vdwbj_s8_camb3lyp=2.0674_dp
 real(dp),parameter :: vdwbj_s8_lcwpbe=1.8541_dp,   vdwbj_s8_mpw1b95=1.0508_dp
 real(dp),parameter :: vdwbj_s8_mpwb1k=0.9499_dp,   vdwbj_s8_mpwlyp=2.0077_dp
 real(dp),parameter :: vdwbj_s8_olyp=2.6205_dp,     vdwbj_s8_opbe=3.3816_dp
 real(dp),parameter :: vdwbj_s8_otpss=2.7495_dp,    vdwbj_s8_pbe38=1.4623_dp
 real(dp),parameter :: vdwbj_s8_pbesol=2.9491_dp,   vdwbj_s8_ptpss=0.2804_dp
 real(dp),parameter :: vdwbj_s8_pwb6k=0.9383_dp,    vdwbj_s8_revssb=0.4389_dp
 real(dp),parameter :: vdwbj_s8_ssb=-0.1744_dp,     vdwbj_s8_tpssh=2.2382_dp
 real(dp),parameter :: vdwbj_s8_hcth120=1.0821_dp,  vdwbj_s8_b2plyp=0.9147_dp
 real(dp),parameter :: vdwbj_s8_b3lyp=1.9889_dp,    vdwbj_s8_b97d=2.2609_dp
 real(dp),parameter :: vdwbj_s8_blyp=2.6996_dp,     vdwbj_s8_bp86=3.2822_dp
 real(dp),parameter :: vdwbj_s8_dsdblyp=0.2130_dp,  vdwbj_s8_pbe0=1.2177_dp
 real(dp),parameter :: vdwbj_s8_pbe=0.7875_dp,      vdwbj_s8_pw6b95=0.7257_dp
 real(dp),parameter :: vdwbj_s8_pwpb95=0.2904_dp,   vdwbj_s8_revpbe0=1.7588_dp
 real(dp),parameter :: vdwbj_s8_revpbe38=1.4760_dp, vdwbj_s8_revpbe=2.3550_dp
 real(dp),parameter :: vdwbj_s8_rpw86pbe=1.3845_dp, vdwbj_s8_tpss0=1.2576_dp
 real(dp),parameter :: vdwbj_s8_tpss=1.9435_dp
! a1 parameters (BJ only)
 real(dp),parameter :: vdwbj_a1_b1b95=0.2092_dp,    vdwbj_a1_b2gpplyp=0.0000_dp
 real(dp),parameter :: vdwbj_a1_b3pw91=0.4312_dp,   vdwbj_a1_bhlyp=0.2793_dp
 real(dp),parameter :: vdwbj_a1_bmk=0.1940_dp,      vdwbj_a1_bop=0.4870_dp
 real(dp),parameter :: vdwbj_a1_bpbe=0.4567_dp,     vdwbj_a1_camb3lyp=0.3708_dp
 real(dp),parameter :: vdwbj_a1_lcwpbe=0.3919_dp,   vdwbj_a1_mpw1b95=0.1955_dp
 real(dp),parameter :: vdwbj_a1_mpwb1k=0.1474_dp,   vdwbj_a1_mpwlyp=0.4831_dp
 real(dp),parameter :: vdwbj_a1_olyp=0.5299_dp,     vdwbj_a1_opbe=0.5512_dp
 real(dp),parameter :: vdwbj_a1_otpss=0.4634_dp,    vdwbj_a1_pbe38=0.3995_dp
 real(dp),parameter :: vdwbj_a1_pbesol=0.4466_dp,   vdwbj_a1_ptpss=0.000_dp
 real(dp),parameter :: vdwbj_a1_pwb6k=0.1805_dp,    vdwbj_a1_revssb=0.4720_dp
 real(dp),parameter :: vdwbj_a1_ssb=-0.0952_dp,     vdwbj_a1_tpssh=0.4529_dp
 real(dp),parameter :: vdwbj_a1_hcth120=0.3563_dp,  vdwbj_a1_b2plyp=0.3065_dp
 real(dp),parameter :: vdwbj_a1_b3lyp=0.3981_dp,    vdwbj_a1_b97d=0.5545_dp
 real(dp),parameter :: vdwbj_a1_blyp=0.4298_dp,     vdwbj_a1_bp86=0.3946_dp
 real(dp),parameter :: vdwbj_a1_dsdblyp=0.000_dp,   vdwbj_a1_pbe0=0.4145_dp
 real(dp),parameter :: vdwbj_a1_pbe=0.4289_dp,      vdwbj_a1_pw6b95=0.2076_dp
 real(dp),parameter :: vdwbj_a1_pwpb95=0.0000_dp,   vdwbj_a1_revpbe0=0.4679_dp
 real(dp),parameter :: vdwbj_a1_revpbe38=0.4309_dp, vdwbj_a1_revpbe=0.5238_dp
 real(dp),parameter :: vdwbj_a1_rpw86pbe=0.4613_dp, vdwbj_a1_tpss0=0.3768_dp
 real(dp),parameter :: vdwbj_a1_tpss=0.4535_dp
! a2 parameters (BJ only)
 real(dp),parameter :: vdwbj_a2_b1b95=5.5545_dp,    vdwbj_a2_b2gpplyp=6.3332_dp
 real(dp),parameter :: vdwbj_a2_b3pw91=4.4693_dp,   vdwbj_a2_bhlyp=4.9615_dp
 real(dp),parameter :: vdwbj_a2_bmk=5.9197_dp,      vdwbj_a2_bop=3.5043_dp
 real(dp),parameter :: vdwbj_a2_bpbe=4.3908_dp,     vdwbj_a2_camb3lyp=5.4743_dp
 real(dp),parameter :: vdwbj_a2_lcwpbe=5.0897_dp,   vdwbj_a2_mpw1b95=6.4177_dp
 real(dp),parameter :: vdwbj_a2_mpwb1k=6.6223_dp,   vdwbj_a2_mpwlyp=4.5323_dp
 real(dp),parameter :: vdwbj_a2_olyp=2.8065_dp,     vdwbj_a2_opbe=2.9444_dp
 real(dp),parameter :: vdwbj_a2_otpss=4.3153_dp,    vdwbj_a2_pbe38=5.1405_dp
 real(dp),parameter :: vdwbj_a2_pbesol=6.1742_dp,   vdwbj_a2_ptpss=6.5745_dp
 real(dp),parameter :: vdwbj_a2_pwb6k=7.7627_dp,    vdwbj_a2_revssb=4.0986_dp
 real(dp),parameter :: vdwbj_a2_ssb=5.2170_dp,      vdwbj_a2_tpssh=4.6550_dp
 real(dp),parameter :: vdwbj_a2_hcth120=4.3359_dp,  vdwbj_a2_b2plyp=5.0570_dp
 real(dp),parameter :: vdwbj_a2_b3lyp=4.4211_dp,    vdwbj_a2_b97d=3.2297_dp
 real(dp),parameter :: vdwbj_a2_blyp=4.2359_dp,     vdwbj_a2_bp86=4.8516_dp
 real(dp),parameter :: vdwbj_a2_dsdblyp=6.0519_dp,  vdwbj_a2_pbe0=4.8593_dp
 real(dp),parameter :: vdwbj_a2_pbe=4.4407_dp,      vdwbj_a2_pw6b95=6.3750_dp
 real(dp),parameter :: vdwbj_a2_pwpb95=7.3141_dp,   vdwbj_a2_revpbe0=3.7619_dp
 real(dp),parameter :: vdwbj_a2_revpbe38=3.9446_dp, vdwbj_a2_revpbe=3.5016_dp
 real(dp),parameter :: vdwbj_a2_rpw86pbe=4.5062_dp, vdwbj_a2_tpss0=4.5865_dp
 real(dp),parameter :: vdwbj_a2_tpss=4.4752_dp
! s6 parameters (zero damping)
 real(dp),parameter :: vdw_s6_b2gpplyp=0.56_dp, vdw_s6_b2plyp=0.64_dp
 real(dp),parameter :: vdw_s6_dsdblyp=0.50_dp,  vdw_s6_ptpss=0.75_dp
 real(dp),parameter :: vdw_s6_pwpb95=0.82_dp
! s8 parameters (zero damping)
 real(dp),parameter :: vdw_s8_b1b95=1.868_dp,   vdw_s8_b2gpplyp=0.760_dp
 real(dp),parameter :: vdw_s8_b3lyp=1.703_dp,   vdw_s8_b97d=0.909_dp
 real(dp),parameter :: vdw_s8_bhlyp=1.442_dp,   vdw_s8_blyp=1.682_dp
 real(dp),parameter :: vdw_s8_bp86=1.683_dp,    vdw_s8_bpbe=2.033_dp
 real(dp),parameter :: vdw_s8_mpwlyp=1.098_dp,  vdw_s8_pbe=0.722_dp
 real(dp),parameter :: vdw_s8_pbe0=0.928_dp,    vdw_s8_pw6b95=0.862_dp
 real(dp),parameter :: vdw_s8_pwb6k=0.550_dp,   vdw_s8_revpbe=1.010_dp
 real(dp),parameter :: vdw_s8_tpss=1.105_dp,    vdw_s8_tpss0=1.242_dp
 real(dp),parameter :: vdw_s8_tpssh=1.219_dp,   vdw_s8_bop=1.975_dp
 real(dp),parameter :: vdw_s8_mpw1b95=1.118_dp, vdw_s8_mpwb1k=1.061_dp
 real(dp),parameter :: vdw_s8_olyp=1.764_dp,    vdw_s8_opbe=2.055_dp
 real(dp),parameter :: vdw_s8_otpss=1.494_dp,   vdw_s8_pbe38=0.998_dp
 real(dp),parameter :: vdw_s8_pbesol=0.612_dp,  vdw_s8_revssb=0.560_dp
 real(dp),parameter :: vdw_s8_ssb=0.663_dp,     vdw_s8_b3pw91=1.775_dp
 real(dp),parameter :: vdw_s8_bmk=2.168_dp,     vdw_s8_camb3lyp=1.217_dp
 real(dp),parameter :: vdw_s8_lcwpbe=1.279_dp,  vdw_s8_m052x=0.00_dp
 real(dp),parameter :: vdw_s8_m05=0.595_dp,     vdw_s8_m062x=0.00_dp
 real(dp),parameter :: vdw_s8_m06hf=0.00_dp,    vdw_s8_m06l=0.00_dp
 real(dp),parameter :: vdw_s8_m06=0.00_dp,      vdw_s8_hcth120=1.206_dp
 real(dp),parameter :: vdw_s8_b2plyp=1.022_dp,  vdw_s8_dsdblyp=0.705_dp
 real(dp),parameter :: vdw_s8_ptpss=0.879_dp,   vdw_s8_pwpb95=0.705_dp
 real(dp),parameter :: vdw_s8_revpbe0=0.792_dp, vdw_s8_revpbe38=0.862_dp
 real(dp),parameter :: vdw_s8_rpw86pbe=0.901_dp
! sr6 parameters (zero damping)
 real(dp),parameter :: vdw_sr6_b1b95=1.613_dp,   vdw_sr6_b2gpplyp=1.586_dp
 real(dp),parameter :: vdw_sr6_b3lyp=1.261_dp,   vdw_sr6_b97d=0.892_dp
 real(dp),parameter :: vdw_sr6_bhlyp=1.370_dp,   vdw_sr6_blyp=1.094_dp
 real(dp),parameter :: vdw_sr6_bp86=1.139_dp,    vdw_sr6_bpbe=1.087_dp
 real(dp),parameter :: vdw_sr6_mpwlyp=1.239_dp,  vdw_sr6_pbe=1.217_dp
 real(dp),parameter :: vdw_sr6_pbe0=1.287_dp,    vdw_sr6_pw6b95=1.532_dp
 real(dp),parameter :: vdw_sr6_pwb6k=1.660_dp,   vdw_sr6_revpbe=0.923_dp
 real(dp),parameter :: vdw_sr6_tpss=1.166_dp,    vdw_sr6_tpss0=1.252_dp
 real(dp),parameter :: vdw_sr6_tpssh=1.223_dp,   vdw_sr6_bop=0.929_dp
 real(dp),parameter :: vdw_sr6_mpw1b95=1.605_dp, vdw_sr6_mpwb1k=1.671_dp
 real(dp),parameter :: vdw_sr6_olyp=0.806_dp,    vdw_sr6_opbe=0.837_dp
 real(dp),parameter :: vdw_sr6_otpss=1.128_dp,   vdw_sr6_pbe38=1.333_dp
 real(dp),parameter :: vdw_sr6_pbesol=1.345_dp,  vdw_sr6_revssb=1.221_dp
 real(dp),parameter :: vdw_sr6_ssb=1.215_dp,     vdw_sr6_b3pw91=1.176_dp
 real(dp),parameter :: vdw_sr6_bmk=1.931_dp,     vdw_sr6_camb3lyp=1.378_dp
 real(dp),parameter :: vdw_sr6_lcwpbe=1.355_dp,  vdw_sr6_m052x=1.417_dp
 real(dp),parameter :: vdw_sr6_m05=1.373_dp,     vdw_sr6_m062x=1.619_dp
 real(dp),parameter :: vdw_sr6_m06hf=1.446_dp,   vdw_sr6_m06l=1.581_dp
 real(dp),parameter :: vdw_sr6_m06=1.325_dp,     vdw_sr6_hcth120=1.221_dp
 real(dp),parameter :: vdw_sr6_b2plyp=1.427_dp,  vdw_sr6_dsdblyp=1.569_dp
 real(dp),parameter :: vdw_sr6_ptpss=1.541_dp,   vdw_sr6_pwpb95=1.557_dp
 real(dp),parameter :: vdw_sr6_revpbe0=0.949_dp, vdw_sr6_revpbe38=1.021_dp
 real(dp),parameter :: vdw_sr6_rpw86pbe=1.224_dp
! sr8 parameters (zero damping) = 1.000_dp

 real(dp),parameter :: vdw_sr9=3.0/4.0
 real(dp),parameter :: vdw_tol_default=tol10
 real(dp) :: ang,arg,cn_dmp,cosa,cosb,cosc,c6,c8
 real(dp) :: dcosa_r3drij,dcosa_r3drjk,dcosa_r3drki
 real(dp) :: dcosb_r3drij,dcosb_r3drjk,dcosb_r3drki
 real(dp) :: dcosc_r3drij,dcosc_r3drjk,dcosc_r3drki
 real(dp) :: dcn_dmp,dexp_cn,dfdmp,dfdmp_drij
 real(dp) :: dfdmp_drjk,dfdmp_drki,dlri,dlrj,dmp,dmp6,dmp8,dmp9,dr,d2lri,d2lrj,d2lrirj
 real(dp) :: dsysref,dsysref_a, dsysref_b
 real(dp) :: d1_r3drij,d1_r3drjk,d1_r3drki,d2cn_dmp,d2cn_exp,d2frac_cn
 real(dp) :: d_drij,d_drjk,d_drki
 real(dp) :: exp_cn,e_no_c6,e_no_c8,e_3bt,fdmp6,fdmp8,fdmp9,frac_cn
 real(dp) :: grad,grad_no_c,grad6,grad6_no_c6,grad8,grad8_no_c8,gr6,gr8
 real(dp) :: hess,hessij, hess6, hess8,im_arg,l,ltot
 real(dp) :: max_vdw_c6,min_dsys,re_arg,rcovij,rcut,rcutcn,rcut2,rcut9
 real(dp) :: rsq,rsqij,rsqjk,rsqki,rmean,rr,rrij,rrjk,rrki,rijk,r0,r6,r8
 real(dp) :: sfact6,sfact8,sfact9,sum_dlri,sum_dlrj,sum_dlc6ri,sum_dlc6rj
 real(dp) :: sum_d2lri,sum_d2lrj,sum_d2lrirj,sum_d2lc6ri,sum_d2lc6rj,sum_d2lc6rirj
 real(dp) :: temp,temp2
 real(dp) :: ucvol,vdw_s6,vdw_s8,vdw_sr6,vdw_sr8,vdw_a1,vdw_a2,vdw_q
 character(len=500) :: msg
 type(atomdata_t) :: atom1,atom2

!arrays

! Covalence radius of the different species for CN (coordination number)
real(dp),parameter:: rcov(vdw_nspecies)=&
&    (/0.80628308, 1.15903197, 3.02356173, 2.36845659, 1.94011865, &
&      1.88972601, 1.78894056, 1.58736983, 1.61256616, 1.68815527, &
&      3.52748848, 3.14954334, 2.84718717, 2.62041997, 2.77159820, &
&      2.57002732, 2.49443835, 2.41884923, 4.43455700, 3.88023730, &
&      3.35111422, 3.07395437, 3.04875805, 2.77159820, 2.69600923, &
&      2.62041997, 2.51963467, 2.49443835, 2.54483100, 2.74640188, &
&      2.82199085, 2.74640188, 2.89757982, 2.77159820, 2.87238349, &
&      2.94797246, 4.76210950, 4.20778980, 3.70386304, 3.50229216, &
&      3.32591790, 3.12434702, 2.89757982, 2.84718717, 2.84718717, &
&      2.72120556, 2.89757982, 3.09915070, 3.22513231, 3.17473967, &
&      3.17473967, 3.09915070, 3.32591790, 3.30072128, 5.26603625, &
&      4.43455700, 4.08180818, 3.70386304, 3.98102289, 3.95582657, &
&      3.93062995, 3.90543362, 3.80464833, 3.82984466, 3.80464833, &
&      3.77945201, 3.75425569, 3.75425569, 3.72905937, 3.85504098, &
&      3.67866672, 3.45189952, 3.30072128, 3.09915070, 2.97316878, &
&      2.92277614, 2.79679452, 2.82199085, 2.84718717, 3.32591790, &
&      3.27552496, 3.27552496, 3.42670319, 3.30072128, 3.47709584, &
&      3.57788113, 5.06446567, 4.56053862, 4.20778980, 3.98102289, &
&      3.82984466, 3.85504098, 3.88023730, 3.90543362 /)

! q = arrays of vdw_species elements containing the link between C6ij and C8ij:
! C8ij = 3sqrt(qi)sqrt(qj)C6ij
 real(dp),parameter :: vdw_q_dftd3(vdw_nspecies)= &
&   (/2.00734898,  1.56637132,  5.01986934,  3.85379032,  3.64446594, &
&     3.10492822,  2.71175247,  2.59361680,  2.38825250,  2.21522516, &
&     6.58585536,  5.46295967,  5.65216669,  4.88284902,  4.29727576, &
&     4.04108902,  3.72932356,  3.44677275,  7.97762753,  7.07623947, &
&     6.60844053,  6.28791364,  6.07728703,  5.54643096,  5.80491167, &
&     5.58415602,  5.41374528,  5.28497229,  5.22592821,  5.09817141, &
&     6.12149689,  5.54083734,  5.06696878,  4.87005108,  4.59089647, &
&     4.31176304,  9.55461698,  8.67396077,  7.97210197,  7.43439917, &
&     6.58711862,  6.19536215,  6.01517290,  5.81623410,  5.65710424, &
&     5.52640661,  5.44263305,  5.58285373,  7.02081898,  6.46815523, &
&     5.98089120,  5.81686657,  5.53321815,  5.25477007, 11.02204549, &
&    10.15679528,  9.35167836,  9.06926079,  8.97241155,  8.90092807, &
&     8.85984840,  8.81736827,  8.79317710,  7.89969626,  8.80588454, &
&     8.42439218,  8.54289262,  8.47583370,  8.45090888,  8.47339339, &
&     7.83525634,  8.20702843,  7.70559063,  7.32755997,  7.03887381, &
&     6.68978720,  6.05450052,  5.88752022,  5.70661499,  5.78450695, &
&     7.79780729,  7.26443867,  6.78151984,  6.67883169,  6.39024318, &
&     6.09527958, 11.79156076, 11.10997644,  9.51377795,  8.67197068, &
&     8.77140725,  8.65402716,  8.53923501,  8.85024712 /)

 integer  :: is(3), nshell_3bt(3)
 integer,allocatable :: ivdw(:)
 integer  :: jmin(3), jmax(3), js(3)
 integer,parameter :: voigt1(6)=(/1,2,3,2,1,1/),voigt2(6)=(/1,2,3,3,3,2/)
 real(dp),allocatable :: cn(:),cfgrad_no_c(:,:,:,:)
 real(dp),allocatable :: dcn(:,:,:),dcn_cart(:,:,:),dc6ri(:,:),dc6rj(:,:),dc9ijri(:,:),dc9ijrj(:,:)
 !real(dp),allocatable:: d2cn(:,:,:,:,:,:)
 real(dp),allocatable:: d2cn_iii(:,:,:,:)
 real(dp),allocatable:: d2cn_jji(:,:,:,:,:)
 real(dp),allocatable:: d2cn_iji(:,:,:,:,:)
 real(dp),allocatable:: d2cn_jii(:,:,:,:,:)
 real(dp),allocatable:: d2cn_tmp(:)
 real(dp),allocatable :: d2c6ri(:,:),d2c6rj(:,:),d2c6rirj(:,:)
 real(dp),allocatable:: elt_cn(:,:,:),e3bt_ij(:,:),e3bt_jk(:,:),e3bt_ki(:,:),e_no_c(:,:)
 real(dp),allocatable:: e_alpha1(:),e_alpha2(:),e_alpha3(:),e_alpha4(:)
 real(dp) :: gred(3),gredij(3),gredjk(3),gredki(3)
 real(dp),allocatable:: fe_no_c(:,:,:),cfdcn(:,:,:,:),fdcn(:,:,:,:),fgrad_no_c(:,:,:,:),gred_vdw_3bt(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3)
 real(dp),allocatable:: grad_no_cij(:,:,:)
 real(dp) :: mcart(3,3)
 real(dp) :: r(3),rcart(3),rcart2(3,3),rcartij(3),rcartjk(3),rcartki(3)
 real(dp):: rij(3), rjk(3), rki(3),rmet(3,3),rred(3)
 real(dp),allocatable:: r0ijk(:,:,:)
 real(dp),allocatable :: str_alpha1(:,:),str_alpha2(:,:),str_dcn(:,:),str_no_c(:,:,:)
 real(dp) :: str_3bt(6)
 real(dp) :: temp_comp(2),temp_comp2(2)
 real(dp),allocatable:: temp_prod(:,:)
 real(dp) :: vec(6),vecij(6), vecjk(6),vecki(6)
 real(dp),allocatable :: vdw_cnrefi(:,:,:,:),vdw_cnrefj(:,:,:,:)
 real(dp),allocatable :: vdw_c6(:,:),vdw_c6ref(:,:,:,:),vdw_c8(:,:),vdw_c9(:,:,:),vdw_r0(:,:)
 real(dp):: vdw_dftd3_r0(4465)
 real(dp):: vdw_dftd3_c6(32385)
 integer:: index_c6(254)
 real(dp):: vdw_dftd3_cni(27884)
 integer:: index_cni(27884)
 real(dp):: vdw_dftd3_cnj(13171)
 integer:: index_cnj(13171)
 real(dp),allocatable :: xred01(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 write(msg,'(1a)')&
& '====> STARTING DFT-D3 computation'
 call wrtout(std_out,msg,'COLL')

 call vdw_dftd3_data(vdw_dftd3_r0,vdw_dftd3_c6,index_c6,vdw_dftd3_cni,index_cni,&
& vdw_dftd3_cnj,index_cnj)
! Determine the properties which have to be studied
 bol_3bt = (vdw_tol_3bt>0)
 need_forces = present(gred_vdw_dftd3)
 need_stress= present(str_vdw_dftd3)
 need_dynmat= present(dyn_vdw_dftd3)
 need_elast= present(elt_vdw_dftd3)
 need_grad=(need_forces.or.need_stress.or.need_dynmat.or.need_elast)
 need_hess=(need_dynmat.or.need_elast)
 if (need_dynmat) then
   if (.not.present(qphon)) then
     msg='Dynamical matrix required without a q-vector'
     ABI_BUG(msg)
   end if
   dyn_vdw_dftd3=zero
 end if
 e_vdw_dftd3 = zero
 if (need_forces) gred_vdw_dftd3=zero
 if (need_stress) str_vdw_dftd3=zero
 if (need_elast) elt_vdw_dftd3=zero

!Identify type(s) of atoms
 ABI_MALLOC(ivdw,(ntypat))
 do itypat=1,ntypat
   call atomdata_from_znucl(atom1,znucl(itypat))
   if (znucl(itypat).gt.94.0_dp) then
     write(msg,'(3a,es14.2)') &
&     'Van der Waals DFT-D3 correction not available for atom type: ',znucl(itypat),' !'
     ABI_ERROR(msg)
   else
     ivdw(itypat) = znucl(itypat)
   end if
 end do

! Determination of coefficients that depend of the
! exchange-correlation functional
 vdw_s6 =one; vdw_s8 =one
 vdw_sr6=one; vdw_sr8=one
 vdw_a1 =one; vdw_a2 =one
! Case one : DFT-D3
 if (vdw_xc == 6) then
   select case (ixc)
   case(11, -130101, -101130)
     vdw_sr6=vdw_sr6_pbe; vdw_s8=vdw_s8_pbe
   case(14, -130102, -102130)
     vdw_sr6=vdw_sr6_revpbe; vdw_s8=vdw_s8_revpbe
   case(18, -131106, -106131)
     vdw_sr6=vdw_sr6_blyp; vdw_s8=vdw_s8_blyp
   case(19, -132106, -106132)
     vdw_sr6=vdw_sr6_bp86; vdw_s8=vdw_s8_bp86
   case(41, -406)
     vdw_sr6=vdw_sr6_pbe0; vdw_s8=vdw_s8_pbe0
   case(-440)
     vdw_sr6=vdw_sr6_b1b95; vdw_s8=vdw_s8_b1b95
   case(-402)
     vdw_sr6=vdw_sr6_b3lyp; vdw_s8=vdw_s8_b3lyp
   case(-170)
     vdw_sr6=vdw_sr6_b97d; vdw_s8=vdw_s8_b97d
   case(-435)
     vdw_sr6=vdw_sr6_bhlyp; vdw_s8=vdw_s8_bhlyp
   case(-130106, -106130)
     vdw_sr6=vdw_sr6_bpbe; vdw_s8=vdw_s8_bpbe
   case(-174)
     vdw_sr6=vdw_sr6_mpwlyp; vdw_s8=vdw_s8_mpwlyp
   case(-451)
     vdw_sr6=vdw_sr6_pw6b95; vdw_s8=vdw_s8_pw6b95
   case(-452)
     vdw_sr6=vdw_sr6_pwb6k; vdw_s8=vdw_s8_pwb6k
   case(-231202, -202231)
     vdw_sr6=vdw_sr6_tpss; vdw_s8=vdw_s8_tpss
   case(-396)
     vdw_sr6=vdw_sr6_tpss0; vdw_s8=vdw_s8_tpss0
   case(-457)
     vdw_sr6=vdw_sr6_tpssh; vdw_s8=vdw_s8_tpssh
   case(-636)
     vdw_sr6=vdw_sr6_bop; vdw_s8=vdw_s8_bop
   case(-445)
     vdw_sr6=vdw_sr6_mpw1b95; vdw_s8=vdw_s8_mpw1b95
   case(-446)
     vdw_sr6=vdw_sr6_mpwb1k; vdw_s8=vdw_s8_mpwb1k
   case(-67)
     vdw_sr6=vdw_sr6_olyp; vdw_s8=vdw_s8_olyp
   case(-65)
     vdw_sr6=vdw_sr6_opbe; vdw_s8=vdw_s8_opbe
   case(-64)
     vdw_sr6=vdw_sr6_otpss; vdw_s8=vdw_s8_otpss
   case(-393)
     vdw_sr6=vdw_sr6_pbe38; vdw_s8=vdw_s8_pbe38
   case(-133116, -116133)
     vdw_sr6=vdw_sr6_pbesol; vdw_s8=vdw_s8_pbesol
   case(-312089, -89312)
     vdw_sr6=vdw_sr6_revssb; vdw_s8=vdw_s8_revssb
   case(-91089, -89091)
     vdw_sr6=vdw_sr6_ssb; vdw_s8=vdw_s8_ssb
   case(-401)
     vdw_sr6=vdw_sr6_b3pw91; vdw_s8=vdw_s8_b3pw91
   case(-280279, -279280)
     vdw_sr6=vdw_sr6_bmk; vdw_s8=vdw_s8_bmk
   case(-433)
     vdw_sr6=vdw_sr6_camb3lyp; vdw_s8=vdw_s8_camb3lyp
   case(-478)
     vdw_sr6=vdw_sr6_lcwpbe; vdw_s8=vdw_s8_lcwpbe
   case(-439238, -238439)
     vdw_sr6=vdw_sr6_m052x; vdw_s8=vdw_s8_m052x
   case(-438237, -237438)
     vdw_sr6=vdw_sr6_m05; vdw_s8=vdw_s8_m05
   case(-450236, -236450)
     vdw_sr6=vdw_sr6_m062x; vdw_s8=vdw_s8_m062x
   case(-444234, -234444)
     vdw_sr6=vdw_sr6_m06hf; vdw_s8=vdw_s8_m06hf
   case(-203233, -233203)
     vdw_sr6=vdw_sr6_m06l; vdw_s8=vdw_s8_m06l
   case(-449235, -235449)
     vdw_sr6=vdw_sr6_m06; vdw_s8=vdw_s8_m06
   case(-162)
     vdw_sr6=vdw_sr6_hcth120; vdw_s8=vdw_s8_hcth120
   case(-456)
     vdw_sr6=vdw_sr6_revpbe0; vdw_s8=vdw_s8_revpbe0
   case(-30108, -108030)
     vdw_sr6=vdw_sr6_rpw86pbe; vdw_s8=vdw_s8_rpw86pbe     
   case default
     write(msg,'(a,i8,a)')'  Van der Waals DFT-D3 correction not compatible with ixc=',ixc,' !'
     ABI_ERROR(msg)
   end select
! Case DFT-D3(BJ)
 elseif (vdw_xc == 7) then
   select case (ixc)
   case(11, -130101, -101130)
     vdw_s8=vdwbj_s8_pbe; vdw_a1=vdwbj_a1_pbe; vdw_a2=vdwbj_a2_pbe
   case(14, -130102, -102130)
     vdw_s8=vdwbj_s8_revpbe; vdw_a1=vdwbj_a1_revpbe; vdw_a2=vdwbj_a2_revpbe
   case(18, -131106, -106131)
     vdw_s8=vdwbj_s8_blyp; vdw_a1=vdwbj_a1_blyp; vdw_a2=vdwbj_a2_blyp
   case(19, -132106, -106132)
     vdw_s8=vdwbj_s8_bp86; vdw_a1=vdwbj_a1_bp86; vdw_a2=vdwbj_a2_bp86
   case(41, -406)
     vdw_s8=vdwbj_s8_pbe0; vdw_a1=vdwbj_a1_pbe0; vdw_a2=vdwbj_a2_pbe0
   case(-440)
     vdw_s8=vdwbj_s8_b1b95; vdw_a1=vdwbj_a1_b1b95; vdw_a2=vdwbj_a2_b1b95
   case(-401)
     vdw_s8=vdwbj_s8_b3pw91; vdw_a1=vdwbj_a1_b3pw91; vdw_a2=vdwbj_a2_b3pw91
   case(-435)
     vdw_s8=vdwbj_s8_bhlyp; vdw_a1=vdwbj_a1_bhlyp; vdw_a2=vdwbj_a2_bhlyp
   case(-280279, -279280)
     vdw_s8=vdwbj_s8_bmk; vdw_a1=vdwbj_a1_bmk; vdw_a2=vdwbj_a2_bmk
   case(-636)
     vdw_s8=vdwbj_s8_bop; vdw_a1=vdwbj_a1_bop; vdw_a2=vdwbj_a2_bop
   case(-130106, -106130)
     vdw_s8=vdwbj_s8_bpbe; vdw_a1=vdwbj_a1_bpbe; vdw_a2=vdwbj_a2_bpbe
   case(-433)
     vdw_s8=vdwbj_s8_camb3lyp; vdw_a1=vdwbj_a1_camb3lyp; vdw_a2=vdwbj_a2_camb3lyp
   case(-478)
     vdw_s8=vdwbj_s8_lcwpbe; vdw_a1=vdwbj_a1_lcwpbe; vdw_a2=vdwbj_a2_lcwpbe
   case(-445)
     vdw_s8=vdwbj_s8_mpw1b95; vdw_a1=vdwbj_a1_mpw1b95; vdw_a2=vdwbj_a2_mpw1b95
   case(-446)
     vdw_s8=vdwbj_s8_mpwb1k; vdw_a1=vdwbj_a1_mpwb1k; vdw_a2=vdwbj_a2_mpwb1k
   case(-174)
     vdw_s8=vdwbj_s8_mpwlyp; vdw_a1=vdwbj_a1_mpwlyp; vdw_a2=vdwbj_a2_mpwlyp
   case(-67)
     vdw_s8=vdwbj_s8_olyp; vdw_a1=vdwbj_a1_olyp; vdw_a2=vdwbj_a2_olyp
   case(-65)
     vdw_s8=vdwbj_s8_opbe; vdw_a1=vdwbj_a1_opbe; vdw_a2=vdwbj_a2_opbe
   case(-64)
     vdw_s8=vdwbj_s8_otpss; vdw_a1=vdwbj_a1_otpss; vdw_a2=vdwbj_a2_otpss
   case(-393)
     vdw_s8=vdwbj_s8_pbe38; vdw_a1=vdwbj_a1_pbe38; vdw_a2=vdwbj_a2_pbe38
   case(-133116, -116133)
     vdw_s8=vdwbj_s8_pbesol; vdw_a1=vdwbj_a1_pbesol; vdw_a2=vdwbj_a2_pbesol
   case(-452)
     vdw_s8=vdwbj_s8_pwb6k; vdw_a1=vdwbj_a1_pwb6k; vdw_a2=vdwbj_a2_pwb6k
   case(-312089, -89312)
     vdw_s8=vdwbj_s8_revssb; vdw_a1=vdwbj_a1_revssb; vdw_a2=vdwbj_a2_revssb
   case(-91089, -89091)
     vdw_s8=vdwbj_s8_ssb; vdw_a1=vdwbj_a1_ssb; vdw_a2=vdwbj_a2_ssb
   case(-457)
     vdw_s8=vdwbj_s8_tpssh; vdw_a1=vdwbj_a1_tpssh; vdw_a2=vdwbj_a2_tpssh
   case(-162)
     vdw_s8=vdwbj_s8_hcth120; vdw_a1=vdwbj_a1_hcth120; vdw_a2=vdwbj_a2_hcth120
   case(-402)
     vdw_s8=vdwbj_s8_b3lyp; vdw_a1=vdwbj_a1_b3lyp; vdw_a2=vdwbj_a2_b3lyp
   case(-170)
     vdw_s8=vdwbj_s8_b97d; vdw_a1=vdwbj_a1_b97d; vdw_a2=vdwbj_a2_b97d
   case(-451)
     vdw_s8=vdwbj_s8_pw6b95; vdw_a1=vdwbj_a1_pw6b95; vdw_a2=vdwbj_a2_pw6b95
   case(-30108, -108030)
     vdw_s8=vdwbj_s8_rpw86pbe; vdw_a1=vdwbj_a1_rpw86pbe; vdw_a2=vdwbj_a2_rpw86pbe
   case(-396)
     vdw_s8=vdwbj_s8_tpss0; vdw_a1=vdwbj_a1_tpss0; vdw_a2=vdwbj_a2_tpss0
   case(-231202, -202231)
     vdw_s8=vdwbj_s8_tpss; vdw_a1=vdwbj_a1_tpss; vdw_a2=vdwbj_a2_tpss
   case default
     write(msg,'(a,i8,a)')'  Van der Waals DFT-D3(BJ) correction not compatible with ixc=',ixc,' !'
     ABI_ERROR(msg)
   end select
 end if

! --------------------------------------------------------------
! Retrieve the data for the referenced c6, cn and r0 coefficients
!---------------------------------------------------------------
 refmax = 5

 ABI_MALLOC(vdw_c6ref,(ntypat,ntypat,refmax,refmax))
 ABI_MALLOC(vdw_cnrefi,(ntypat,ntypat,refmax,refmax))
 ABI_MALLOC(vdw_cnrefj,(ntypat,ntypat,refmax,refmax))
 ABI_MALLOC(vdw_r0,(ntypat,ntypat))
 if (bol_3bt) then
   ABI_MALLOC(r0ijk,(ntypat,ntypat,ntypat))
 end if

 vdw_c6ref = zero
 vdw_cnrefi = 100 ; vdw_cnrefj = 100 ;

 do refi=1,refmax
   do refj=1,refi
     do itypat=1,ntypat
       do jtypat=1,ntypat
         indi = ivdw(itypat)+100*(refi-1)
         indj = ivdw(jtypat)+100*(refj-1)
         found = .false.
         do ia=1,size(index_c6)
           do ja=1,size(index_c6)
             if (index_c6(ia)==indi.and.index_c6(ja)==indj) then
               if (ia>=ja)then
                 nline = ia*(ia-1)/2 + ja
               else
                 nline = ja*(ja-1)/2 + ia
               endif
               vdw_c6ref(itypat,jtypat,refi,refj) = vdw_dftd3_c6(nline)
               vdw_c6ref(jtypat,itypat,refj,refi) = vdw_dftd3_c6(nline)
               found = .false.
               do la=1,size(index_cni)
                 if (index_cni(la)==nline) then
                   found=.true.
                   vdw_cnrefi(itypat,jtypat,refi,refj)= vdw_dftd3_cni(la)
                   vdw_cnrefj(jtypat,itypat,refj,refi)= vdw_dftd3_cni(la)
                 else
                   vdw_cnrefi(itypat,jtypat,refi,refj) = zero
                   vdw_cnrefj(jtypat,itypat,refj,refi) = zero
                 end if
                 if (found) exit
               end do
               found = .false.
               do la=1,size(index_cnj)
                 if (index_cnj(la)==nline) then
                   found=.true.
                   vdw_cnrefj(itypat,jtypat,refi,refj)= vdw_dftd3_cnj(la)
                   vdw_cnrefi(jtypat,itypat,refj,refi)= vdw_dftd3_cnj(la)
                 else
                   vdw_cnrefj(itypat,jtypat,refi,refj) = zero
                   vdw_cnrefi(jtypat,itypat,refj,refi) = zero
                 end if
                 if (found) exit
               end do
               found = .true.
             end if
             if (found) exit
           end do
           if (found) exit
         end do
         if (refi.eq.1.and.refj.eq.1) then
           nline = ia*(ia-1)/2 + ja
           vdw_r0(itypat,jtypat)=vdw_dftd3_r0(nline)/Bohr_Ang
           if (bol_3bt) then
             do ktypat=1,ntypat
               r0ijk(itypat,jtypat,ktypat)=one/(vdw_r0(itypat,jtypat)*vdw_r0(jtypat,ktypat)*vdw_r0(ktypat,itypat))**third
             end do ! ka atom
           end if    ! Only if 3bt required
         end if       ! Only for the first set of references
       end do          ! Loop on references j
     end do             ! Loop on references i
   end do                ! Loop on atom j
 end do                   ! Loop on atom i

 !if (vdw_d3_cov==1) then
 !   vdw_cnrefi(:,:,:,refmax) =vdw_cnrefi(:,:,:,refmax-1)
 !   vdw_cnrefi(:,:,refmax,:) = 14.0_dp
 !   vdw_cnrefj(:,:,refmax,:) =vdw_cnrefj(:,:,refmax-1,:)
 !   vdw_cnrefj(:,:,:,refmax) = 14.0_dp
 !end if
!Retrieve cell geometry data
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Map reduced coordinates into [0,1[
 ABI_MALLOC(xred01,(3,natom))
 do ia=1,natom
   xred01(:,ia)=xred(:,ia)-aint(xred(:,ia)) ! Map into ]-1,1[
   do alpha=1,3
     if (abs(xred01(alpha,ia)).ge.tol8) xred01(alpha,ia) = xred01(alpha,ia)+half-sign(half,xred(alpha,ia))
   end do
 end do

! -------------------------------------------------------------------
! Computation of the coordination number (CN) for the different atoms
! -------------------------------------------------------------------

 write(msg,'(3a)')&
& '  Begin the computation of the Coordination Numbers (CN)',ch10,&
& '  required for DFT-D3 energy corrections...'
 call wrtout(std_out,msg,'COLL')

! Allocation of the CN coefficients and derivatives
 ABI_MALLOC(cn,(natom))
 ABI_MALLOC(dcn,(3,natom,natom))
 ABI_MALLOC(dcn_cart,(3,natom,natom))
 ABI_MALLOC(str_dcn,(6,natom))
! ABI_MALLOC_OR_DIE(d2cn, (2,3,natom,3,natom,natom), ierr)
 ABI_MALLOC_OR_DIE(d2cn_iii, (2,3,3,natom), ierr)
 ABI_MALLOC_OR_DIE(d2cn_jji, (2,3,3,natom,natom), ierr)
 ABI_MALLOC_OR_DIE(d2cn_iji, (2,3,3,natom,natom), ierr)
 ABI_MALLOC_OR_DIE(d2cn_jii, (2,3,3,natom,natom), ierr)
 ABI_MALLOC(d2cn_tmp, (2))
 ABI_MALLOC(fdcn,(2,3,natom,natom))
 ABI_MALLOC(cfdcn,(2,3,natom,natom))
 ABI_MALLOC(elt_cn,(6+3*natom,6,natom))

! Initializing the computed quantities to zero
 nshell = 0
 cn = zero
! Initializing the derivative of the computed quantities to zero (if required)
 dcn = zero ; str_dcn = zero
 d2cn_iii = zero
 d2cn_jji = zero
 d2cn_iji = zero
 d2cn_jii = zero
 fdcn = zero; cfdcn = zero
 elt_cn = zero ; dcn_cart = zero

 re_arg = zero ; im_arg = zero
 if (need_hess) then
   mcart = zero
   do alpha=1,3
     mcart(alpha,alpha) = one
   end do
 end if
 rcutcn = 200**2 ! Bohr
!Loop over shells of cell replicas
 do
   newshell=.false.;nshell=nshell+1
!   Loop over cell replicas in the shell
   do is3=-nshell,nshell
     do is2=-nshell,nshell
       do is1=-nshell,nshell
         if (nshell==1.or. &
&         abs(is3)==nshell.or.abs(is2)==nshell.or.abs(is1)==nshell) then
           is(3) = is3 ; is(2) = is2 ; is(1) = is1
!               Computation of phase factor for discrete Fourier transform
!               if phonon at qphon is required
           if (need_dynmat) then
             arg=two_pi*dot_product(qphon,is)
             re_arg=cos(arg) ; im_arg=sin(arg)
           end if
!               Loop over atoms ia and ja
           do ia=1,natom
             itypat=typat(ia)
             do ja=1,natom
               jtypat=typat(ja)
               r(:)=xred01(:,ia)-xred01(:,ja)-dble(is(:))
               rsq=dot_product(r,matmul(rmet,r))
!                     atom i =/= j
               if (rsq.ge.tol16.and.rsq<rcutcn) then
                 newshell=.true.
                 rr = sqrt(rsq)
                 rcovij = rcov(ivdw(itypat))+rcov(ivdw(jtypat))

!                         Computation of partial contribution to cn coefficients
                 exp_cn = exp(-k1*(rcovij/rr-one))
                 frac_cn= one/(one+exp_cn)
!                         Introduction of a damping function for the coordination
!                         number because of the divergence with increasing
!                         number of cells of this quantity in periodic systems
!                         See Reckien et al., J. Comp. Chem. 33, 2023 (2012) [[cite:Reckien2012]]
                 dr = rr-k2*rcovij
                 cn_dmp = half*abi_derfc(dr)
                 cn(ia) = cn(ia)+frac_cn*cn_dmp

!                         If force, stress, IFC or Elastic constants are required,
!                         computation of the first derivative of CN
                 if (need_grad) then
                   rcart=matmul(rprimd,r)
                   dexp_cn= k1*rcovij*exp_cn/rsq
                   dcn_dmp = -one/sqrt(pi)*exp(-dr*dr)
                   grad=(-frac_cn*frac_cn*cn_dmp*dexp_cn+dcn_dmp*frac_cn)/rr
                   if (need_forces.and.ia/=ja) then
!                               Variation of CN(ia) w.r. to displacement of atom
!                               ja. If ia==ka then all the other atoms contribute
!                               to the derivative. Required for the computation of
!                               the forces applied on atom k
                     rred = matmul(transpose(rprimd),rcart)
                     dcn(:,ia,ia) = dcn(:,ia,ia)+grad*rred(:)
                     dcn(:,ia,ja) = dcn(:,ia,ja)-grad*rred(:)
                   elseif (need_stress.or.need_elast) then
!                               The following quantity (str_dcn) is used for the computation
!                               of the DFT-D3 contribution to stress and elastic constants
                     vec(1:3)=rcart(1:3)*rcart(1:3); vec(4)=rcart(2)*rcart(3)
                     vec(5)=rcart(1)*rcart(3); vec(6)=rcart(1)*rcart(2)
                     str_dcn(:,ia)=str_dcn(:,ia)+grad*vec(:)
                   end if
!                            If dynamical matrix or elastic constants are required, compute
!                            the second derivative
                   if (need_hess) then
                     d2cn_dmp = two*(rr-k2*rcovij)/sqrt(pi)*exp(-(rr-k2*rcovij)**two)
                     d2cn_exp = dexp_cn*(k1*rcovij/rsq-two/rr)
                     d2frac_cn =frac_cn**two*(two*frac_cn*dexp_cn**two-d2cn_exp)
                     hess = (d2frac_cn*cn_dmp+d2cn_dmp*frac_cn-&
&                     two*dcn_dmp*frac_cn**two*dexp_cn-grad)/rsq
                     if (need_dynmat) then
!                                  Discrete Fourier Transform of dCN/drk in cartesian
!                                  coordinates. Note that the phase factor is different
!                                  if (ka=ia or ka/=ia).
!                                  This Fourier transform is summed over cell replica
!                                  See NOTE: add reference for more information
                       fdcn(1,:,ia,ia) = fdcn(1,:,ia,ia)+grad*rcart(:)
                       fdcn(1,:,ia,ja) = fdcn(1,:,ia,ja)-grad*rcart(:)*re_arg
                       fdcn(2,:,ia,ja) = fdcn(2,:,ia,ja)-grad*rcart(:)*im_arg
!                      Conjugate of fdcn
                       cfdcn(1,:,ia,ia) = cfdcn(1,:,ia,ia)+grad*rcart(:)
                       cfdcn(1,:,ia,ja) = cfdcn(1,:,ia,ja)-grad*rcart(:)*re_arg
                       cfdcn(2,:,ia,ja) = cfdcn(2,:,ia,ja)+grad*rcart(:)*im_arg
!                                  Computation of second derivative of CN required for the
!                                  interatomic force constants in reciprocal space
                       do alpha=1,3
                         rcart2(alpha,:) = rcart(alpha)*rcart(:)
                       end do
!                                  Computation of second derivative of CN required for the
!                                  interatomic force constants in reciprocal space
!                                  This Fourier transform is summed over cell replica
!                                  as it appears in the theory
                       do alpha=1,3
                         if (ia/=ja) then
                                         ! ka = ia ; la = ia
                           d2cn_iii(1,alpha,:,ia) = d2cn_iii(1,alpha,:,ia)+&
&                           (hess*rcart2(alpha,:)+grad*mcart(alpha,:))
                                         ! ka = ja ; la = ja
                           d2cn_jji(1,alpha,:,ja,ia) = d2cn_jji(1,alpha,:,ja,ia)+&
&                           (hess*rcart2(alpha,:)+grad*mcart(alpha,:))
                           if (abs(re_arg)>tol12) then
                                            ! ka = ia ; la = ja
                             d2cn_iji(1,alpha,:,ja,ia) = d2cn_iji(1,alpha,:,ja,ia)-&
&                             (hess*rcart2(alpha,:)+grad*mcart(alpha,:))*re_arg
                                            ! ka = ja ; la = ia
                             d2cn_jii(1,alpha,:,ja,ia) = d2cn_jii(1,alpha,:,ja,ia)-&
&                             (hess*rcart2(alpha,:)+grad*mcart(alpha,:))*re_arg
                           end if
                           if (abs(im_arg)>tol12) then
                                            ! ka = ia ; la = ja
                             d2cn_iji(2,alpha,:,ja,ia) = d2cn_iji(2,alpha,:,ja,ia)-&
&                             (hess*rcart2(alpha,:)+grad*mcart(alpha,:))*im_arg
                                            ! ka = ja ; la = ia
                             d2cn_jii(2,alpha,:,ja,ia) = d2cn_jii(2,alpha,:,ja,ia)+&
&                             (hess*rcart2(alpha,:)+grad*mcart(alpha,:))*im_arg
                           end if
                         else
                           if (abs(re_arg-one)>tol12) then
                             d2cn_iji(1,alpha,:,ja,ia) = d2cn_iji(1,alpha,:,ja,ia)+&
&                             two*(hess*rcart2(alpha,:)+grad*mcart(alpha,:))*(one-re_arg)
                           end if
                         end if
                       end do
                     end if ! Boolean Need_dynmat
                     if (need_elast) then
!                                  Derivative of str_dcn w.r. to strain for elastic tensor
                       vec(1:3)=rcart(1:3)*rcart(1:3); vec(4)=rcart(2)*rcart(3)
                       vec(5)=rcart(1)*rcart(3); vec(6)=rcart(1)*rcart(2)
                       do alpha=1,6
                         ii = voigt1(alpha) ; jj=voigt2(alpha)
                         do beta=1,6
                           kk = voigt1(beta) ; ll=voigt2(beta)
                           elt_cn(alpha,beta,ia) = elt_cn(alpha,beta,ia)+hess*rcart(ii)*rcart(jj)*rcart(kk)*rcart(ll)
                           if (ii==kk) elt_cn(alpha,beta,ia) = elt_cn(alpha,beta,ia)+half*grad*rcart(jj)*rcart(ll)
                           if (jj==kk) elt_cn(alpha,beta,ia) = elt_cn(alpha,beta,ia)+half*grad*rcart(ii)*rcart(ll)
                           if (ii==ll) elt_cn(alpha,beta,ia) = elt_cn(alpha,beta,ia)+half*grad*rcart(jj)*rcart(kk)
                           if (jj==ll) elt_cn(alpha,beta,ia) = elt_cn(alpha,beta,ia)+half*grad*rcart(ii)*rcart(kk)
                         end do
                       end do
                                   ! Derivative of str_dcn w.r. to atomic displacement
                                   ! for internal strains
                       dcn_cart(:,ia,ia) = dcn_cart(:,ia,ia)+grad*rcart(:)
                       dcn_cart(:,ia,ja) = dcn_cart(:,ia,ja)-grad*rcart(:)
                       if (ia/=ja) then
                         index_ia = 6+3*(ia-1)
                         index_ja = 6+3*(ja-1)
                         do alpha=1,6
                           do beta=1,3
                             ii = voigt1(alpha) ; jj=voigt2(alpha)
                             elt_cn(index_ia+beta,alpha,ia)=elt_cn(index_ia+beta,alpha,ia)&
&                             +hess*vec(alpha)*rcart(beta)
                             elt_cn(index_ja+beta,alpha,ia)=elt_cn(index_ja+beta,alpha,ia)&
&                             -hess*vec(alpha)*rcart(beta)
                             if (ii==beta) then
                               elt_cn(index_ia+beta,alpha,ia)=elt_cn(index_ia+beta,alpha,ia)+grad*rcart(jj)
                               elt_cn(index_ja+beta,alpha,ia)=elt_cn(index_ja+beta,alpha,ia)-grad*rcart(jj)
                             end if
                             if (jj==beta) then
                               elt_cn(index_ia+beta,alpha,ia)=elt_cn(index_ia+beta,alpha,ia)+grad*rcart(ii)
                               elt_cn(index_ja+beta,alpha,ia)=elt_cn(index_ja+beta,alpha,ia)-grad*rcart(ii)
                             end if
                           end do
                         end do
                       end if ! ia/=ja
                     end if ! Need strain derivative
                   end if ! Boolean second derivative
                 end if !  Boolean first derivative
               end if ! Tolerence
             end do ! Loop over ia atom
           end do ! Loop over ja atom
         end if ! Bondary Condition
       end do ! Loop over is1
     end do ! Loop over is2
   end do ! Loop over is3
   if(.not.newshell) exit ! Check if a new shell must be considered
 end do ! Loop over shell
 write(msg,'(3a,f8.5,1a,i3,1a,f8.5,1a,i3,1a)')&
& '                                            ... Done.',ch10,&
& '  max(CN) =', maxval(cn), ' (atom ',maxloc(cn),') ;  min(CN) =', minval(cn), ' (atom ', minloc(cn),')'
 call wrtout(std_out,msg,'COLL')

!----------------------------------------------------------------
! Computation of the C6 coefficient
! ---------------------------------------------------------------

 write(msg,'(3a)')&
& '  Begin the computation of the C6(CN)',ch10,&
& '  required for DFT-D3 energy corrections...'
 call wrtout(std_out,msg,'COLL')
 ! Allocation
 ABI_MALLOC(vdw_c6,(natom,natom))
 ABI_MALLOC(vdw_c8,(natom,natom))
 ABI_MALLOC(dc6ri,(natom,natom))
 ABI_MALLOC(dc6rj,(natom,natom))
 ABI_MALLOC(d2c6ri,(natom,natom))
 ABI_MALLOC(d2c6rj,(natom,natom))
 ABI_MALLOC(d2c6rirj,(natom,natom))
 if (bol_3bt) then
   ABI_MALLOC(vdw_c9,(natom,natom,natom))
   ABI_MALLOC(dc9ijri,(natom,natom))
   ABI_MALLOC(dc9ijrj,(natom,natom))
 end if
! Set accumulating quantities to zero
 vdw_c6 = zero ; vdw_c8 = zero
 dc6ri = zero ; dc6rj = zero
 d2c6ri = zero ; d2c6rj = zero
 d2c6rirj = zero
 if (bol_3bt) then
   dc9ijri = zero; dc9ijrj = zero
 end if
! C6 coefficients are interpolated from tabulated
! ab initio C6 values (following loop).
! C8 coefficients are obtained by:
! C8 = vdw_dftd3_q(itypat)*vdw_dftd3_q(jtypat)*C6

 do ia=1,natom
   itypat=typat(ia)
   do ja=1,natom
     jtypat=typat(ja)
!      Set accumulating quantities to zero
     ltot=zero
     sum_dlri = zero ; sum_dlc6ri= zero
     sum_dlrj = zero ; sum_dlc6rj= zero
     sum_d2lri = zero ; sum_d2lc6ri = zero
     sum_d2lrj = zero ; sum_d2lc6rj = zero
     sum_d2lrirj = zero ; sum_d2lc6rirj = zero
     min_dsys = 10000
     max_vdw_c6 = zero
!      Loop over references
     do refi=1,refmax
       do refj=1,refmax
         dsysref_a = cn(ia)-vdw_cnrefi(itypat,jtypat,refi,refj)
         dsysref_b = cn(ja)-vdw_cnrefj(itypat,jtypat,refi,refj)
         dsysref=(dsysref_a)**two+(dsysref_b)**two
         if (dsysref<min_dsys) then
!               Keep in memory the smallest value of dsysref
!               And the associated tabulated C6 value
           min_dsys = dsysref
           max_vdw_c6 = vdw_c6ref(itypat,jtypat,refi,refj)
         end if
         l = dexp(-k3*dsysref)
         ltot = ltot+l
         vdw_c6(ia,ja)=vdw_c6(ia,ja)+vdw_c6ref(itypat,jtypat,refi,refj)*l

         if (need_grad) then
!               Derivative of l(ia,ja) with respect to the displacement
!               of atom ka in reduced coordinates.
!               This factor is identical in the case of stress.
!               In purpose of speed up this routine, the prefactor of
!               dCNi/drk and dCNj/drk are separated
!               See NOTE: article to be added
           dlri=-k3*l*two*dsysref_a ;dlrj=-k3*l*two*dsysref_b
           sum_dlri=sum_dlri+dlri ; sum_dlrj=sum_dlrj+dlrj
           sum_dlc6ri=sum_dlc6ri+dlri*vdw_c6ref(itypat,jtypat,refi,refj)
           sum_dlc6rj=sum_dlc6rj+dlrj*vdw_c6ref(itypat,jtypat,refi,refj)
           if (need_hess) then
!                  Second derivative of l(ia,ja). Once again, it is separately in
!                  different contributions:
!                  d2lri: prefactor of dCNi/drk*dCNi/drl
!                  d2lrj: prefactor of dCNj/drk*dCNj/drl
!                  d2lrirj: prefacto of dCNi/drk*dCNj/drl
!                  The prefactor for d2CNi/drkdrl is dlri; for d2CNj/drkdrl is dlrj
             d2lri = -two*k3*l*(one-two*k3*dsysref_a**two)
             d2lrj = -two*k3*l*(one-two*k3*dsysref_b**two)
             d2lrirj = four*k3*k3*l*(dsysref_a*dsysref_b)
             sum_d2lri=sum_d2lri+d2lri ; sum_d2lrj=sum_d2lrj+d2lrj
             sum_d2lrirj = sum_d2lrirj+d2lrirj
             sum_d2lc6ri=sum_d2lc6ri+d2lri*vdw_c6ref(itypat,jtypat,refi,refj)
             sum_d2lc6rj=sum_d2lc6rj+d2lrj*vdw_c6ref(itypat,jtypat,refi,refj)
             sum_d2lc6rirj = sum_d2lc6rirj + d2lrirj*vdw_c6ref(itypat,jtypat,refi,refj)
           end if ! Boolean second derivative
         end if ! Boolean gradient
       end do ! Loop over references
     end do ! Loop over references
!      In some specific case (really covalently bound compounds) ltot -> 0
!      which may cause numerical problems for all quantities related to dispersion coefficient.
!      To be consistent with VASP implementation, the c6 value is taken as the last
!      referenced value of the dispersion coefficient.
     if (ltot>tol12) then
       vdw_c6(ia,ja)=vdw_c6(ia,ja)/ltot
       vdw_c8(ia,ja)=three*vdw_q_dftd3(ivdw(itypat))*vdw_q_dftd3(ivdw(jtypat))*vdw_c6(ia,ja)
!         If Force of Stress is required
       if (need_grad) then
!            Computation of the derivative of C6 w.r.to the displacement
!            of atom ka, in reduced coordinates (separated for dCNi/drk and dCNj/drk)
!            This is the crucial step to reduce the scaling from O(N^3) to O(N^2) for
!            the gradients
         dc6ri(ia,ja)=(sum_dlc6ri-vdw_c6(ia,ja)*sum_dlri)/ltot
         dc6rj(ia,ja)=(sum_dlc6rj-vdw_c6(ia,ja)*sum_dlrj)/ltot
         if (need_hess) then
!               Computation of the second derivative of C6 w.r.to the displacement of atom ka
!               and atom la
           d2c6ri(ia,ja)=(sum_d2lc6ri-vdw_c6(ia,ja)*sum_d2lri-two*dc6ri(ia,ja)*sum_dlri)/ltot
           d2c6rj(ia,ja)=(sum_d2lc6rj-vdw_c6(ia,ja)*sum_d2lrj-two*dc6rj(ia,ja)*sum_dlrj)/ltot
           d2c6rirj(ia,ja) = (sum_d2lc6rirj-vdw_c6(ia,ja)*sum_d2lrirj-dc6ri(ia,ja)*&
&           sum_dlrj-dc6rj(ia,ja)*sum_dlri)/ltot
         end if ! Boolean second derivative
       end if ! Boolean gradient
     else
       vdw_c6(ia,ja)= max_vdw_c6
       vdw_c8(ia,ja)=three*vdw_q_dftd3(ivdw(itypat))*vdw_q_dftd3(ivdw(jtypat))*vdw_c6(ia,ja)
     end if
   end do
 end do
! Computation of the three-body term dispersion coefficient
 if (bol_3bt) then
   do ia=1,natom
     do ja=1,natom
       do ka=1,natom
         vdw_c9(ia,ja,ka) =-sqrt(vdw_c6(ia,ja)*vdw_c6(ja,ka)*vdw_c6(ka,ia))
       end do
       if (need_grad) then
         dc9ijri(ia,ja) = half/vdw_c6(ia,ja)*dc6ri(ia,ja)
         dc9ijrj(ia,ja) = half/vdw_c6(ia,ja)*dc6rj(ia,ja)
       end if
     end do
   end do
 end if

 write(msg,'(3a,f10.5,1a,f10.5)')&
& '                                            ... Done.',ch10,&
& '  max(C6) =', maxval(vdw_c6),' ;  min(C6) =', minval(vdw_c6)
 call wrtout(std_out,msg,'COLL')

! Deallocation of used variables not needed anymore
 ABI_FREE(vdw_c6ref)
 ABI_FREE(vdw_cnrefi)
 ABI_FREE(vdw_cnrefj)

!----------------------------------------------------
! Computation of cut-off radii according to tolerance
!----------------------------------------------------

 ! Cut-off radius for pair-wise term
 if (vdw_tol<zero) then
   rcut=max((vdw_s6/vdw_tol_default*maxval(vdw_c6))**sixth, &
&   (vdw_s8/vdw_tol_default*maxval(vdw_c8))**(one/eight))
 else
   rcut=max((vdw_s6/vdw_tol*maxval(vdw_c6))**sixth,&
&   (vdw_s8/vdw_tol*maxval(vdw_c8))**(one/eight))
 end if
 ! Cut-off radius for three-body term
 rcut9 = zero
 if (bol_3bt) then
   rcut9=(128.0_dp*vdw_s6/(vdw_tol_3bt)*maxval(vdw_c6)**(3.0/2.0))**(1.0/9.0)
 end if
 rcut2=rcut*rcut

!--------------------------------------------------------------------
! Computation of the two bodies contribution to the dispersion energy
!--------------------------------------------------------------------

 write(msg,'(3a)')&
& '  Begin the computation of pair-wise term',ch10,&
& '  of DFT-D3 energy contribution...'
 call wrtout(std_out,msg,'COLL')
 nshell=0
 npairs=0
 ABI_MALLOC(e_alpha1,(natom))
 ABI_MALLOC(e_alpha2,(natom))
 ABI_MALLOC(e_alpha3,(natom))
 ABI_MALLOC(e_alpha4,(natom))
 ABI_MALLOC(e_no_c,(natom,natom))
 ABI_MALLOC(fe_no_c,(2,natom,natom))
 e_alpha1 =zero ; e_alpha2 = zero
 e_alpha3 =zero ; e_alpha4 = zero
 e_no_c=zero  ; fe_no_c = zero
 ABI_MALLOC(grad_no_cij,(3,natom,natom))
 ABI_MALLOC(fgrad_no_c,(2,3,natom,natom))
 ABI_MALLOC(cfgrad_no_c,(2,3,natom,natom))
 ABI_MALLOC(str_no_c,(6,natom,natom))
 ABI_MALLOC(str_alpha1,(6,natom))
 ABI_MALLOC(str_alpha2,(6,natom))
 grad_no_cij=zero ; str_no_c=zero
 fgrad_no_c = zero ;  cfgrad_no_c = zero
 str_alpha1 = zero ; str_alpha2 = zero

 re_arg = zero ; im_arg = zero
 dmp6 = zero ; dmp8 = zero
 e_no_c6 = zero ; e_no_c8 = zero
 fdmp6 = zero ; fdmp8 = zero
 grad6 = zero ; grad8 = zero
 grad6_no_c6 = zero ; grad8_no_c8 = zero
 hess6 = zero ; hess8 = zero
 do
   newshell=.false.;nshell=nshell+1
   do is3=-nshell,nshell
     do is2=-nshell,nshell
       do is1=-nshell,nshell
         if (nshell==1.or.abs(is3)==nshell.or.abs(is2)==nshell.or.abs(is1)==nshell) then
           is(1) = is1 ; is(2) = is2 ; is(3) = is3
!              Computation of phase factor for discrete Fourier transform
!              if phonon at qphon is required
           if (need_dynmat) then
             arg= two_pi*dot_product(qphon,is)
             re_arg=cos(arg)
             im_arg=sin(arg)
           end if
           do ia=1,natom
             itypat=typat(ia)
             do ja=1,natom
               jtypat=typat(ja)
               r(:)=xred01(:,ia)-xred01(:,ja)-dble(is(:))
               rsq=dot_product(r,matmul(rmet,r))
               if (rsq>=tol16.and.rsq<rcut2) then
                 npairs=npairs+1;newshell=.true.
                 sfact6 = half*vdw_s6 ; sfact8 = half*vdw_s8
                 rr=sqrt(rsq); r6 = rr**six ; r8 = rr**eight
                 c6=vdw_c6(ia,ja) ; c8=vdw_c8(ia,ja)
                 vdw_q = three*vdw_q_dftd3(ivdw(itypat))*vdw_q_dftd3(ivdw(jtypat))
                 r0=vdw_r0(itypat,jtypat)
!                        Computation of e_vdw_dftd3 (case DFT+D3)
                 if (vdw_xc == 6) then
                   dmp6=six*(rr/(vdw_sr6*r0))**(-alpha6)
                   fdmp6=one/(one+dmp6)
                   dmp8=six*(rr/(vdw_sr8*r0))**(-alpha8)
                   fdmp8=one/(one+dmp8)
!                           Contribution to energy
                   e_no_c6 = -sfact6*fdmp6/r6 ; e_no_c8 = -sfact8*fdmp8/r8
                   e_vdw_dftd3=e_vdw_dftd3+e_no_c6*c6 +e_no_c8*c8
!                        Computation of e_vdw_dftd3 (case DFT+D3-BJ)
                 elseif (vdw_xc == 7) then
                   dmp = (vdw_a1*sqrt(vdw_q)+vdw_a2)
                   fdmp6 = one/(dmp**six+rr**six)
                   fdmp8 = one/(dmp**eight+rr**eight)
                   e_no_c6 = -sfact6*fdmp6 ; e_no_c8 = -sfact8*fdmp8
                   e_vdw_dftd3=e_vdw_dftd3-sfact6*c6*fdmp6-sfact8*c8*fdmp8
                 end if
!                        Computation of the gradients (if required)
                 if (need_grad) then
                   if (vdw_xc == 6) then
                     gr6 = alpha6*dmp6*fdmp6**two
                     grad6_no_c6 = sfact6*(gr6-six*fdmp6)/r8
                     grad6 = grad6_no_c6*c6
                     gr8 = alpha8*dmp8*fdmp8**two
                     grad8_no_c8 = sfact8*(gr8-eight*fdmp8)/r8/rsq
                     grad8 = grad8_no_c8*c8
                   elseif (vdw_xc == 7) then
                     grad6_no_c6 = -sfact6*six*(fdmp6*rsq)**two
                     grad6 = grad6_no_c6*c6
                     grad8_no_c8 = -sfact8*eight*(fdmp8)**two*rsq**three
                     grad8 = grad8_no_c8*c8
                   end if
                   grad =grad6+grad8
                   grad_no_c = grad6_no_c6+grad8_no_c8*vdw_q
                   rcart=matmul(rprimd,r)
                   rred= matmul(transpose(rprimd),rcart)
!                           Additional contribution due to c6(cn(r))
!                           Not yet multiply by dCN/drk and summed to reduce
!                           computational time
                   e_no_c(ia,ja) = e_no_c(ia,ja)+(e_no_c6+vdw_q*e_no_c8)
!                           Part related to alpha1ij/alpha2ij
                   e_alpha1(ia) = e_alpha1(ia)+(e_no_c6+vdw_q*e_no_c8)*dc6ri(ia,ja)
                   e_alpha2(ja) = e_alpha2(ja)+(e_no_c6+vdw_q*e_no_c8)*dc6rj(ia,ja)
!                           Contribution to gradients wr to atomic displacement
!                           (forces)
                   if (need_forces.and.ia/=ja) then
                     gred(:)=grad*rred(:)
                     do alpha=1,3
                       gred_vdw_dftd3(alpha,ia)=gred_vdw_dftd3(alpha,ia)-gred(alpha)
                       gred_vdw_dftd3(alpha,ja)=gred_vdw_dftd3(alpha,ja)+gred(alpha)
                     end do
                   elseif (need_stress) then
!                              Computation of the DFT-D3 contribution to stress
                     vec(1:3)=rcart(1:3)*rcart(1:3); vec(4)=rcart(2)*rcart(3)
                     vec(5)=rcart(1)*rcart(3); vec(6)=rcart(1)*rcart(2)
                     do alpha=1,6
                       str_vdw_dftd3(alpha)=str_vdw_dftd3(alpha)-grad*vec(alpha)
                     end do
                   end if
!                           Second derivative (if required)
                   if (need_hess) then
                     if (vdw_xc==6) then
                       hess6 = (grad6*(alpha6*fdmp6*dmp6-8.0_dp)+&
&                       sfact6*c6/r6*dmp6*((alpha6*fdmp6)**two)*&
&                       (fdmp6*dmp6-one)/rsq)/rsq
                       hess8 = (grad8*(alpha8*fdmp8*dmp8-10.0_dp)+&
&                       sfact8*c8/r8*dmp8*((alpha8*fdmp8)**two)*&
&                       (fdmp8*dmp8-one)/rsq)/rsq
                     elseif (vdw_xc==7) then
                       hess6 = -four*grad6*(three*rsq**two*fdmp6-one/rsq)
                       hess8 = -two*grad8*(eight*rsq**three*fdmp8-three/rsq)
                     end if
!                              Contribution of d2C6 to the interatomic force constants
!                              Not yet multiply by CN derivative and summed to reduce the scaling from O(N^3) to O(N^2)
                     hessij = hess6+hess8
!                              Contribution of cross-derivative dC6 and grad
                     do alpha=1,3
                       grad_no_cij(alpha,ia,ja) = grad_no_cij(alpha,ia,ja) - grad_no_c*rcart(alpha)
                     end do
                     e_alpha3(ia) = e_alpha3(ia)+(e_no_c6+vdw_q*e_no_c8)*d2c6ri(ia,ja)
                     e_alpha4(ja) = e_alpha4(ja)+(e_no_c6+vdw_q*e_no_c8)*d2c6rj(ia,ja)
                     if (need_dynmat) then
!                                 Fourier transform of the partial contribution to the dispersion potential
                       fe_no_c(1,ia,ja) = fe_no_c(1,ia,ja)+(e_no_c6+vdw_q*e_no_c8)*re_arg
                       fe_no_c(2,ia,ja) = fe_no_c(2,ia,ja)+(e_no_c6+vdw_q*e_no_c8)*im_arg
                       do alpha=1,3
!                                    Fourier transform of the gradient (required for the IFCs)
                         fgrad_no_c(1,alpha,ia,ja) = fgrad_no_c(1,alpha,ia,ja)-grad_no_c*rcart(alpha)*re_arg
                         fgrad_no_c(2,alpha,ia,ja) = fgrad_no_c(2,alpha,ia,ja)-grad_no_c*rcart(alpha)*im_arg
!                                    Complex conjugated of the Fourier transform of the gradient
                         cfgrad_no_c(1,alpha,ia,ja) = cfgrad_no_c(1,alpha,ia,ja)-grad_no_c*rcart(alpha)*re_arg
                         cfgrad_no_c(2,alpha,ia,ja) = cfgrad_no_c(2,alpha,ia,ja)+grad_no_c*rcart(alpha)*im_arg
                       end do
!                                 Contribution to the IFCs (reciprocal space) of the 2nd derivative of e_no_c part
                       do alpha=1,3
                         do beta=1,3
                           rcart2(alpha,beta) = rcart(alpha)*rcart(beta)
                         end do
                       end do
                       if (ia/=ja) then
                         do alpha=1,3
                           dyn_vdw_dftd3(1,alpha,ja,:,ja) = dyn_vdw_dftd3(1,alpha,ja,:,ja) -&
&                           (hessij*rcart2(alpha,:)+grad*mcart(alpha,:))
                           dyn_vdw_dftd3(1,alpha,ia,:,ia) = dyn_vdw_dftd3(1,alpha,ia,:,ia) -&
&                           (hessij*rcart2(alpha,:)+grad*mcart(alpha,:))
                           if (abs(re_arg)>tol12) then
                             dyn_vdw_dftd3(1,alpha,ia,:,ja) = dyn_vdw_dftd3(1,alpha,ia,:,ja) +&
&                             (hessij*rcart2(alpha,:)+grad*mcart(alpha,:))*re_arg
                             dyn_vdw_dftd3(1,alpha,ja,:,ia) = dyn_vdw_dftd3(1,alpha,ja,:,ia) +&
&                             (hessij*rcart2(alpha,:)+grad*mcart(alpha,:))*re_arg
                           end if
                           if (abs(im_arg)>tol12) then
                             dyn_vdw_dftd3(2,alpha,ia,:,ja) = dyn_vdw_dftd3(2,alpha,ia,:,ja) +&
&                             (hessij*rcart2(alpha,:)+grad*mcart(alpha,:))*im_arg
                             dyn_vdw_dftd3(2,alpha,ja,:,ia) = dyn_vdw_dftd3(2,alpha,ja,:,ia) -&
&                             (hessij*rcart2(alpha,:)+grad*mcart(alpha,:))*im_arg
                           end if
                         end do
                       else ! ia==ja
                         do alpha=1,3
                           if (abs(re_arg-one)>tol12) then
                             dyn_vdw_dftd3(1,alpha,ia,:,ia) = dyn_vdw_dftd3(1,alpha,ia,:,ia) -&
&                             two*(hessij*rcart2(alpha,:)+grad*mcart(alpha,:))*(one-re_arg)
                           end if
                         end do
                       end if
                     end if
!                              Now compute the contribution to the elastic constants !!! Still under development
                     if (need_elast) then
                       vec(1:3)=rcart(1:3)*rcart(1:3); vec(4)=rcart(2)*rcart(3)
                       vec(5)=rcart(1)*rcart(3); vec(6)=rcart(1)*rcart(2)
                       str_no_c(:,ia,ja)=str_no_c(:,ia,ja)-grad_no_c*vec(:)
                       str_alpha1(:,ia)=str_alpha1(:,ia)-dc6ri(ia,ja)*grad_no_c*vec(:)
                       str_alpha2(:,ja)=str_alpha2(:,ja)-dc6rj(ia,ja)*grad_no_c*vec(:)
!                                 Contribution to elastic constants of DFT-D3 dispersion potential (no C6 derivative)
                       do alpha=1,6
                         ii = voigt1(alpha) ; jj=voigt2(alpha)
                         do beta=1,6
                           kk = voigt1(beta) ; ll=voigt2(beta)
                           elt_vdw_dftd3(alpha,beta) = elt_vdw_dftd3(alpha,beta)-hessij*vec(alpha)*vec(beta)
                           if (ii==kk) elt_vdw_dftd3(alpha,beta) = elt_vdw_dftd3(alpha,beta)-half*grad*rcart(jj)*rcart(ll)
                           if (jj==kk) elt_vdw_dftd3(alpha,beta) = elt_vdw_dftd3(alpha,beta)-half*grad*rcart(ii)*rcart(ll)
                           if (ii==ll) elt_vdw_dftd3(alpha,beta) = elt_vdw_dftd3(alpha,beta)-half*grad*rcart(jj)*rcart(kk)
                           if (jj==ll) elt_vdw_dftd3(alpha,beta) = elt_vdw_dftd3(alpha,beta)-half*grad*rcart(ii)*rcart(kk)
                         end do
                       end do
!                                 Contribution to internal strain of DFT-D3 dispersion potential (no C6 derivative)
                       if (ia/=ja) then
                         index_ia = 6+3*(ia-1)
                         index_ja = 6+3*(ja-1)
                         do alpha=1,6
                           do beta=1,3
                             ii = voigt1(alpha) ; jj=voigt2(alpha)
                             elt_vdw_dftd3(index_ia+beta,alpha)=elt_vdw_dftd3(index_ia+beta,alpha)-&
&                             hessij*vec(alpha)*rcart(beta)
                             elt_vdw_dftd3(index_ja+beta,alpha)=elt_vdw_dftd3(index_ja+beta,alpha)+&
&                             hessij*vec(alpha)*rcart(beta)
                             if (ii==beta) then
                               elt_vdw_dftd3(index_ia+beta,alpha)=elt_vdw_dftd3(index_ia+beta,alpha)-grad*rcart(jj)
                               elt_vdw_dftd3(index_ja+beta,alpha)=elt_vdw_dftd3(index_ja+beta,alpha)+grad*rcart(jj)
                             end if
                             if (jj==beta) then
                               elt_vdw_dftd3(index_ia+beta,alpha)=elt_vdw_dftd3(index_ia+beta,alpha)-grad*rcart(ii)
                               elt_vdw_dftd3(index_ja+beta,alpha)=elt_vdw_dftd3(index_ja+beta,alpha)+grad*rcart(ii)
                             end if
                           end do ! Direction beta
                         end do    ! Strain alpha
                       end if       ! ia/=ja
                     end if          ! Need elastic constant
                   end if             ! Need hessian
                 end if                ! Need gradient
               end if                   ! Tolerance
             end do                      ! Loop over atom j
           end do                         ! Loop over atom i
         end if                            ! Triple loop over cell replicas in shell
       end do                               ! Is1
     end do                                 ! Is2
   end do                                    ! Is3
   if(.not.newshell) exit ! Check if new shell must be calculated
 end do ! Loop over shell
 ABI_MALLOC(temp_prod,(2,natom))
 if (need_grad) then
!   Additional contribution to force due dc6_drk
   if (need_forces) then
     do ka=1,natom
       do ia=1,natom
             !do ja=1,natom
                !gred_vdw_dftd3(:,ka) = gred_vdw_dftd3(:,ka)+e_no_c(ia,ja)*(&
!&               !dcn(:,ia,ka)*dc6ri(ia,ja)+dcn(:,ja,ka)*dc6rj(ia,ja))
             !end do
         gred_vdw_dftd3(:,ka) = gred_vdw_dftd3(:,ka)+e_alpha1(ia)*dcn(:,ia,ka)+&
&         e_alpha2(ia)*dcn(:,ia,ka)
       end do
     end do
   elseif (need_stress) then
     do ia=1,natom
          !do ja=1,natom
          !   str_vdw_dftd3(:) = str_vdw_dftd3(:)+e_no_c(ia,ja)*(str_dcn(:,ia)*&
!&         !   dc6ri(ia,ja)+str_dcn(:,ja)*dc6rj(ia,ja))
          !end do
       str_vdw_dftd3(:) = str_vdw_dftd3(:)+e_alpha1(ia)*str_dcn(:,ia)+&
&       e_alpha2(ia)*str_dcn(:,ia)
     end do
   end if ! Optimization
!   If dynmat is required, add all the terms related to dc6, d2c6, ...
   if (need_hess) then
     if (need_dynmat) then
       do ka=1,natom
         do la =1,natom
           do alpha=1,3
             do beta=1,3
               do ia=1,natom
!TODO: avoid stupid if clauses inside the loops
                 d2cn_tmp = zero
                 if (ia==la) then 
                   if (ia==ka) then    ! iii
                     d2cn_tmp(:) = d2cn_iii(:,alpha,beta,ia)
                   else                ! jii
                     d2cn_tmp(:) = d2cn_jii(:,alpha,beta,ka,ia)
                   end if
                 else if (ia==ka) then !iji
                   d2cn_tmp(:) = d2cn_iji(:,alpha,beta,la,ia)
                 else if (ka==la) then    ! jji
                   d2cn_tmp(:) = d2cn_jji(:,alpha,beta,ka,ia)
                 end if 
!                 Add the second derivative of C6 contribution to the dynamical matrix
!                 First, add the second derivative of CN-related term
                 dyn_vdw_dftd3(:,alpha,ka,beta,la)=dyn_vdw_dftd3(:,alpha,ka,beta,la)+&
&                 (e_alpha1(ia)+e_alpha2(ia))*d2cn_tmp(:)
!OLDVERSION                 dyn_vdw_dftd3(:,alpha,ka,beta,la)=dyn_vdw_dftd3(:,alpha,ka,beta,la)+&
!&                 (e_alpha1(ia)+e_alpha2(ia))*d2cn(:,alpha,ka,beta,la,ia)



!                 Then the term related to dCNi/dr*dCNi/dr and dCNj/dr*dCNj/dr
                 call comp_prod(cfdcn(:,alpha,ia,ka),fdcn(:,beta,ia,la),temp_comp)
                 dyn_vdw_dftd3(:,alpha,ka,beta,la)=dyn_vdw_dftd3(:,alpha,ka,beta,la)+&
&                 (e_alpha3(ia)+e_alpha4(ia))*temp_comp(:)
!                 Add the cross derivative of fdmp/rij**6 and C6 contribution to the dynamical matrix
!                 !!!! The products are kind of tricky...
!                 First, add the dCNk/drl gradik and dCNk/drl gradjk terms...
                 dyn_vdw_dftd3(:,alpha,ka,beta,la)=dyn_vdw_dftd3(:,alpha,ka,beta,la)+&
&                 cfdcn(:,alpha,la,ka)*grad_no_cij(beta,la,ia)*(dc6ri(la,ia)+dc6rj(ia,la))
!                Then the dCNk/drl gradjk and dCNk/drl gradik terms...
                 call comp_prod(cfdcn(:,alpha,ia,ka),cfgrad_no_c(:,beta,la,ia),temp_comp)
                 dyn_vdw_dftd3(:,alpha,ka,beta,la)=dyn_vdw_dftd3(:,alpha,ka,beta,la)+&
&                 temp_comp(:)*(dc6ri(ia,la)+dc6rj(la,ia))
!                Here the symmetrical term (for dCNl/drk) are added...
                 dyn_vdw_dftd3(:,alpha,ka,beta,la)=dyn_vdw_dftd3(:,alpha,ka,beta,la)+&
&                 fdcn(:,beta,ka,la)*grad_no_cij(alpha,ka,ia)*(dc6ri(ka,ia)+dc6rj(ia,ka))
                 call comp_prod(fdcn(:,beta,ia,la),fgrad_no_c(:,alpha,ka,ia),temp_comp)
                 dyn_vdw_dftd3(:,alpha,ka,beta,la)=dyn_vdw_dftd3(:,alpha,ka,beta,la)+&
&                 temp_comp(:)*(dc6ri(ia,ka)+dc6rj(ka,ia))
               end do ! ia
             end do ! alpha
           end do ! beta
         end do ! la
         do alpha=1,3
           temp_prod(:,:) = zero
           do ja=1,natom
             do ia=1,natom
!             Finally the cross derivative dCNi/dr dCNj/dr
               temp_comp2(:) = d2c6rirj(ia,ja)*fe_no_c(:,ia,ja)
               call comp_prod(cfdcn(:,alpha,ia,ka),temp_comp2,temp_comp)
               temp_prod(:,ja) = temp_prod(:,ja)+temp_comp(:)
             end do
             do la = 1,natom
               do beta=1,3
                 call comp_prod(fdcn(:,beta,ja,la),temp_prod(:,ja),temp_comp2)
                 dyn_vdw_dftd3(:,alpha,ka,beta,la)=dyn_vdw_dftd3(:,alpha,ka,beta,la)+&
&                 two*temp_comp2
               end do ! beta
             end do ! la
           end do ! ja
         end do ! alpha
       end do !ka
!               Transformation from cartesian coordinates to reduced coordinates
       do ka=1,natom
         do la=1,natom
           do kk=1,2
             do alpha=1,3
               vec(1:3)=dyn_vdw_dftd3(kk,1:3,ka,alpha,la)
               call d3_cart2red(vec)
               dyn_vdw_dftd3(kk,1:3,ka,alpha,la)=vec(1:3)
             end do
             do alpha=1,3
               vec(1:3)=dyn_vdw_dftd3(kk,alpha,ka,1:3,la)
               call d3_cart2red(vec)
               dyn_vdw_dftd3(kk,alpha,ka,1:3,la)=vec(1:3)
             end do ! alpha
           end do ! real/im
         end do       ! Atom la
       end do          ! Atom ka
     end if             ! Boolean dynamical matrix
     if (need_elast) then
       do ia=1,natom
         index_ia = 6+3*(ia-1)
         do alpha=1,6
!          Add the second derivative of C6 contribution to the elastic tensor
!          First, the second derivative of CN with strain
           elt_vdw_dftd3(alpha,:) = elt_vdw_dftd3(alpha,:)+(&
&           e_alpha1(ia)+e_alpha2(ia))*elt_cn(alpha,:,ia)
!          Then, the derivative product dCNi/deta dCNi/deta
           elt_vdw_dftd3(alpha,:) = elt_vdw_dftd3(alpha,:)+(&
&           e_alpha3(ia)+e_alpha4(ia))*str_dcn(alpha,ia)*str_dcn(:,ia)
!          Then, the dCNi/deta dCNj/deta
           do ja=1,natom
             elt_vdw_dftd3(alpha,:)=elt_vdw_dftd3(alpha,:)+two*e_no_c(ia,ja)*&
&             d2c6rirj(ia,ja)*str_dcn(alpha,ia)*str_dcn(:,ja)
           end do
!          Add the cross derivative of fij and C6 contribution to the elastic tensor
           elt_vdw_dftd3(alpha,:)=elt_vdw_dftd3(alpha,:)+str_alpha1(:,ia)*str_dcn(alpha,ia)+str_alpha1(alpha,ia)*str_dcn(:,ia)&
&           +str_dcn(alpha,ia)*str_alpha2(:,ia)+str_dcn(:,ia)*str_alpha2(alpha,ia)
         end do
         do alpha=1,6
           ii = voigt1(alpha) ; jj=voigt2(alpha)
           do ka=1,natom
             index_ka = 6+3*(ka-1)
             do beta=1,3
               ! Add the second derivative of C6 contribution to the internal strains
               ii = voigt1(alpha) ; jj=voigt2(alpha)
               ! Second derivative of CN
               elt_vdw_dftd3(index_ka+beta,alpha)=elt_vdw_dftd3(index_ka+beta,alpha)+(e_alpha1(ia)+e_alpha2(ia))*&
&               elt_cn(index_ka+beta,alpha,ia)
               ! Cross-derivatives of CN
               elt_vdw_dftd3(index_ka+beta,alpha)=elt_vdw_dftd3(index_ka+beta,alpha)+(e_alpha3(ia)+e_alpha4(ia))*&
&               dcn_cart(beta,ia,ka)*str_dcn(alpha,ia) !OK
               do ja=1,natom
                 index_ja = 6+3*(ja-1)
                 elt_vdw_dftd3(index_ka+beta,alpha)=elt_vdw_dftd3(index_ka+beta,alpha)+two*d2c6rirj(ia,ja)*e_no_c(ia,ja)*&
&                 dcn_cart(beta,ia,ka)*str_dcn(alpha,ja)
               end do
               elt_vdw_dftd3(index_ka+beta,alpha)=elt_vdw_dftd3(index_ka+beta,alpha)+(str_alpha1(alpha,ia)+str_alpha2(alpha,ia))*&
&               dcn_cart(beta,ia,ka)
               elt_vdw_dftd3(index_ka+beta,alpha)=elt_vdw_dftd3(index_ka+beta,alpha)+grad_no_cij(beta,ka,ia)*&
&               (dc6ri(ia,ka)*str_dcn(alpha,ia)+dc6rj(ia,ka)*str_dcn(alpha,ka))-grad_no_cij(beta,ia,ka)*&
&               (dc6ri(ka,ia)*str_dcn(alpha,ka)+dc6rj(ka,ia)*str_dcn(alpha,ia))
             end do ! beta
           end do ! ka
         end do ! alpha
       end do !  ia
       do alpha=1,6
         index_ia=6
         do ia=1,natom
           elt_vdw_dftd3(index_ia+1:index_ia+3,alpha)=matmul(transpose(rprimd),elt_vdw_dftd3(index_ia+1:index_ia+3,alpha))
           index_ia=index_ia+3
         end do       ! Atom ia
       end do          ! Strain alpha
     end if             ! Boolean elastic tensor
   end if                ! Boolean hessian
 end if                   ! Boolean need_gradient
 ABI_FREE(temp_prod)
 write(msg,'(3a)')&
& '                                  ...Done.'
 call wrtout(std_out,msg,'COLL')

!print *, 'Evdw', e_vdw_dftd3
!if (need_forces) print *, 'fvdw', gred_vdw_dftd3(3,:)
!if (need_stress) print *, 'strvdw', str_vdw_dftd3(6)
!if (need_elast) print *, 'Elast(3,3,3,3)', elt_vdw_dftd3(3,3)
!if (need_elast) print *, 'Elast(3,3,3,3)', elt_vdw_dftd3(6,6)
!if (need_elast) print *, 'Elast(3,3,3,3)', elt_vdw_dftd3(3,6)
!if (need_elast) print *, 'Internal(1,3,3)', elt_vdw_dftd3(9,3)
!if (need_elast) print *, 'Internal(1,3,6)', elt_vdw_dftd3(9,6)
!if (need_dynmat) print *, 'Dynmat(1,3,:,3,:)', dyn_vdw_dftd3(1,3,:,3,:)
!if (need_dynmat) print *, 'Dynmat(1,3,:,3,:)', dyn_vdw_dftd3(1,2,:,1,:)

!---------------------------------------------
! Computation of the 3 body term (if required)
!---------------------------------------------

 e_3bt=zero
 if (bol_3bt) then
   if (need_grad) then
     ABI_MALLOC(e3bt_ij,(natom,natom))
     ABI_MALLOC(e3bt_jk,(natom,natom))
     ABI_MALLOC(e3bt_ki,(natom,natom))
     ABI_MALLOC(gred_vdw_3bt,(3,natom))
     e3bt_ij=zero; e3bt_jk=zero; e3bt_ki=zero
     if (need_forces) then
       gred_vdw_3bt = zero
     elseif (need_stress) then
       str_3bt=zero
     end if
   end if
   nshell_3bt(1) =  int(0.5+rcut9/sqrt(rmet(1,1)+rmet(2,1)+rmet(3,1)))
   nshell_3bt(2) =  int(0.5+rcut9/sqrt(rmet(1,2)+rmet(2,2)+rmet(3,2)))
   nshell_3bt(3) =  int(0.5+rcut9/sqrt(rmet(1,3)+rmet(2,3)+rmet(3,3)))

   do is3 = -nshell_3bt(3),nshell_3bt(3)
     do is2 = -nshell_3bt(2),nshell_3bt(2)
       do is1 = -nshell_3bt(1),nshell_3bt(1)
         is(1) = is1 ; is(2)=is2 ; is(3) = is3
         do alpha=1,3
           jmin(alpha) = max(-nshell_3bt(alpha), -nshell_3bt(alpha)+is(alpha))
           jmax(alpha) = min(nshell_3bt(alpha), nshell_3bt(alpha)+is(alpha))
         end do
         do js3=jmin(3),jmax(3)
           do js2=jmin(2),jmax(2)
             do js1=jmin(1),jmax(1)
               js(1) = js1 ; js(2)=js2 ; js(3) = js3
               do ia=1,natom
                 itypat=typat(ia)
                 do ja=1,ia
                   jtypat=typat(ja)
                   do ka=1,ja
                     ktypat=typat(ka)
                     rij(:) = xred01(:,ia)-xred01(:,ja)-dble(is(:))
                     rsqij = dot_product(rij(:),matmul(rmet,rij(:)))
                     rrij = dsqrt(rsqij)
                     rjk(:) = xred01(:,ja)-xred01(:,ka)+dble(is(:))-dble(js(:))
                     rsqjk = dot_product(rjk(:),matmul(rmet,rjk(:)))
                     rrjk = dsqrt(rsqjk)
                     rki(:) = xred01(:,ka)-xred01(:,ia)+dble(js(:))
                     rsqki = dot_product(rki(:),matmul(rmet,rki(:)))
                     rrki = dsqrt(rsqki)
                     if (rsqij>=tol16.and.rsqjk>=tol16.and.rsqki>=tol16) then
                       rmean = (rrij*rrjk*rrki)**third
                       if (rrij>rcut9.or.rrjk>rcut9.or.rrki>rcut9) cycle
                       sfact9=vdw_s6
                       if (ia==ja.and.ja==ka) sfact9 = sixth*sfact9
                       if (ia==ja.and.ja/=ka) sfact9 = half*sfact9
                       if (ia/=ja.and.ja==ka) sfact9 = half*sfact9
                       rijk = one/(rrij*rrjk*rrki)
                       dmp9 = six*(rmean*vdw_sr9*r0ijk(itypat,jtypat,ktypat))**(-alpha8)
                       fdmp9 = one/(one+dmp9)
                       cosa = half*rrjk*(rsqij+rsqki-rsqjk)*rijk
                       cosb = half*rrki*(rsqij+rsqjk-rsqki)*rijk
                       cosc = half*rrij*(rsqjk+rsqki-rsqij)*rijk
                       ang =  one+three*cosa*cosb*cosc
                       temp = sfact9*rijk*rijk*rijk
                       temp2 = temp*fdmp9*ang
!                                     Contribution to energy
                       e_3bt = e_3bt-temp2*vdw_c9(ia,ja,ka) !*temp2
                       e3bt_ij(ia,ja) = e3bt_ij(ia,ja)-temp2*vdw_c9(ia,ja,ka)
                       e3bt_jk(ja,ka) = e3bt_jk(ja,ka)-temp2*vdw_c9(ia,ja,ka)
                       e3bt_ki(ka,ia) = e3bt_ki(ka,ia)-temp2*vdw_c9(ia,ja,ka)
                       if (need_grad) then
                         dfdmp = third*alpha8*fdmp9*fdmp9*dmp9
                         if (ia/=ja.or.need_stress) then
                           d1_r3drij = -three*rrki*rrjk
                           dcosa_r3drij = (rrij-two*cosa*rrki)*rrjk
                           dcosb_r3drij = (rrij-two*cosb*rrjk)*rrki
                           dcosc_r3drij =-(rsqij+cosc*rrjk*rrki)
                           dfdmp_drij = dfdmp*rrki*rrjk
                           d_drij = vdw_c9(ia,ja,ka)*temp*rijk*((d1_r3drij+three*&
&                           (cosb*cosc*dcosa_r3drij+cosa*cosc*dcosb_r3drij+cosa*cosb*dcosc_r3drij))*&
&                           fdmp9+dfdmp_drij*ang)*rijk*rrjk*rrki
                           rcartij=matmul(rprimd,rij)
                         end if
                         if (ja/=ka.or.need_stress) then
                           d1_r3drjk = -three*rrij*rrki
                           dcosa_r3drjk =-(rsqjk+cosa*rrij*rrki)
                           dcosb_r3drjk = (rrjk-two*cosb*rrij)*rrki
                           dcosc_r3drjk = (rrjk-two*cosc*rrki)*rrij
                           dfdmp_drjk = dfdmp*rrij*rrki
                           d_drjk = vdw_c9(ia,ja,ka)*temp*rijk*((d1_r3drjk+three*&
&                           (cosb*cosc*dcosa_r3drjk+cosa*cosc*dcosb_r3drjk+cosa*cosb*dcosc_r3drjk))*&
&                           fdmp9+dfdmp_drjk*ang)*rijk*rrij*rrki
                           rcartjk=matmul(rprimd,rjk)
                         end if
                         if (ka/=ia.or.need_stress) then
                           d1_r3drki = -three*rrjk*rrij
                           dcosa_r3drki = (rrki-two*cosa*rrij)*rrjk
                           dcosb_r3drki =-(rsqki+cosb*rrij*rrjk)
                           dcosc_r3drki = (rrki-two*cosc*rrjk)*rrij
                           dfdmp_drki = dfdmp*rrij*rrjk
                           d_drki = vdw_c9(ia,ja,ka)*temp*rijk*((d1_r3drki+three*&
&                           (cosb*cosc*dcosa_r3drki+cosa*cosc*dcosb_r3drki+cosa*cosb*dcosc_r3drki))*&
&                           fdmp9+dfdmp_drki*ang)*rijk*rrij*rrjk
                           rcartki=matmul(rprimd,rki)
                         end if
!                                        Contribution to gradients wr to atomic displacement
!                                        (forces)
                         if (need_forces) then
                           if (ia/=ja) gredij=d_drij*matmul(transpose(rprimd),rcartij)
                           if (ja/=ka) gredjk=d_drjk*matmul(transpose(rprimd),rcartjk)
                           if (ka/=ia) gredki=d_drki*matmul(transpose(rprimd),rcartki)
                           if (ia/=ja.and.ka/=ia) then
                             gred_vdw_3bt(:,ia)=gred_vdw_3bt(:,ia)-gredij(:)+gredki(:)
                             gred_vdw_3bt(:,ja)=gred_vdw_3bt(:,ja)+gredij(:)-gredjk(:)
                             gred_vdw_3bt(:,ka)=gred_vdw_3bt(:,ka)-gredki(:)+gredjk(:)
                           else if (ia==ja.and.ia/=ka) then
                             gred_vdw_3bt(:,ia)=gred_vdw_3bt(:,ia)+gredki(:)
                             gred_vdw_3bt(:,ja)=gred_vdw_3bt(:,ja)-gredjk(:)
                             gred_vdw_3bt(:,ka)=gred_vdw_3bt(:,ka)-gredki(:)+gredjk(:)
                           elseif (ia==ka.and.ia/=ja) then
                             gred_vdw_3bt(:,ia)=gred_vdw_3bt(:,ia)-gredij(:)
                             gred_vdw_3bt(:,ja)=gred_vdw_3bt(:,ja)+gredij(:)-gredjk(:)
                             gred_vdw_3bt(:,ka)=gred_vdw_3bt(:,ka)+gredjk(:)
                           elseif (ja==ka.and.ia/=ja) then
                             gred_vdw_3bt(:,ia)=gred_vdw_3bt(:,ia)-gredij(:)+gredki(:)
                             gred_vdw_3bt(:,ja)=gred_vdw_3bt(:,ja)+gredij(:)
                             gred_vdw_3bt(:,ka)=gred_vdw_3bt(:,ka)-gredki(:)
                           end if
                         end if
!                                        Contribution to stress tensor
                         if (need_stress) then
                           vecij(1:3)=rcartij(1:3)*rcartij(1:3); vecij(4)=rcartij(2)*rcartij(3)
                           vecij(5)=rcartij(1)*rcartij(3); vecij(6)=rcartij(1)*rcartij(2)
                           vecjk(1:3)=rcartjk(1:3)*rcartjk(1:3); vecjk(4)=rcartjk(2)*rcartjk(3)
                           vecjk(5)=rcartjk(1)*rcartjk(3); vecjk(6)=rcartjk(1)*rcartjk(2)
                           vecki(1:3)=rcartki(1:3)*rcartki(1:3); vecki(4)=rcartki(2)*rcartki(3)
                           vecki(5)=rcartki(1)*rcartki(3); vecki(6)=rcartki(1)*rcartki(2)
                           str_3bt(:)=str_3bt(:)-d_drij*vecij(:)-d_drjk*vecjk(:)-d_drki*vecki(:) !-str_3bt_dcn9(:)
                         end if
                       end if ! Optimization
                     end if   ! Tolerance
                   end do  ! Loop over atom k
                 end do     ! Loop over atom j
               end do        ! Loop over atom i
             end do ! j3
           end do ! j2
         end do ! j1
       end do
     end do
   end do
   if (need_forces) then
     do ia=1,natom
       do ja=1,natom
         do la=1,natom
           gred_vdw_3bt(:,la) = gred_vdw_3bt(:,la)+e3bt_ij(ia,ja)*(dc9ijri(ia,ja)*dcn(:,ia,la)+dc9ijrj(ia,ja)*dcn(:,ja,la))
           gred_vdw_3bt(:,la) = gred_vdw_3bt(:,la)+e3bt_jk(ia,ja)*(dc9ijri(ia,ja)*dcn(:,ia,la)+dc9ijrj(ia,ja)*dcn(:,ja,la))
           gred_vdw_3bt(:,la) = gred_vdw_3bt(:,la)+e3bt_ki(ia,ja)*(dc9ijri(ia,ja)*dcn(:,ia,la)+dc9ijrj(ia,ja)*dcn(:,ja,la))
         end do
       end do
     end do
   elseif (need_stress) then
     do ia=1,natom
       do ja=1,natom
         str_3bt(:) = str_3bt(:)+e3bt_ij(ia,ja)*(dc9ijri(ia,ja)*str_dcn(:,ia)+dc9ijrj(ia,ja)*str_dcn(:,ja))
         str_3bt(:) = str_3bt(:)+e3bt_jk(ia,ja)*(dc9ijri(ia,ja)*str_dcn(:,ia)+dc9ijrj(ia,ja)*str_dcn(:,ja))
         str_3bt(:) = str_3bt(:)+e3bt_ki(ia,ja)*(dc9ijri(ia,ja)*str_dcn(:,ia)+dc9ijrj(ia,ja)*str_dcn(:,ja))
       end do
     end do
   end if
   e_vdw_dftd3 = e_vdw_dftd3+e_3bt
   if (need_forces) gred_vdw_dftd3= gred_vdw_dftd3+gred_vdw_3bt
   if (need_stress) str_vdw_dftd3 = str_vdw_dftd3+str_3bt
   ABI_FREE(dc9ijri)
   ABI_FREE(dc9ijrj)
   ABI_FREE(e3bt_ij)
   ABI_FREE(e3bt_jk)
   ABI_FREE(e3bt_ki)
   ABI_FREE(vdw_c9)
   ABI_FREE(r0ijk)
   ABI_FREE(gred_vdw_3bt)
 end if
 if (need_stress) str_vdw_dftd3=str_vdw_dftd3/ucvol

!Printing
 if (prtvol>0) then
   write(msg,'(2a)') ch10,&
&   '  --------------------------------------------------------------'
   call wrtout(std_out,msg,'COLL')
   if (vdw_xc==6) then
     write(msg,'(3a)') &
&     '   Van der Waals DFT-D3 semi-empirical dispersion potential as',ch10,&
&     '   proposed by Grimme et al., J. Chem. Phys. 132, 154104 (2010)' ! [[cite:Grimme2010]]
     call wrtout(std_out,msg,'COLL')
   elseif (vdw_xc==7) then
     write(msg,'(5a)') &
&     '    Van der Waals DFT-D3 semi-empirical dispersion potential  ' ,ch10,&
&     '    with Becke-Jonhson (BJ) refined by Grimme et al. J.    ',ch10,&
&     '    Comput. Chem. 32, 1456 (2011) ' ! [[cite:Grimme2011]]
     call wrtout(std_out,msg,'COLL')
   end if
   if (natom<5) then
     write(msg,'(3a)')&
&     '         Pair       C6 (a.u.)       C8 (a.u.)       R0 (Ang)  ',ch10,&
&     '  ---------------------------------------------------------------'
     call wrtout(std_out,msg,'COLL')
     do ia=1,natom
       do ja=1,ia
         itypat = typat(ia) ; jtypat = typat(ja)
         call atomdata_from_znucl(atom1,znucl(itypat))
         call atomdata_from_znucl(atom2,znucl(jtypat))
         write(msg,'(4X,2a,i2,3a,i2,1a,1X,es12.4,4X,es12.4,4X,es12.4,1X)') &
         atom1%symbol,'(',ia,')-',atom2%symbol,'(',ja,')', &
         vdw_c6(ia,ja), vdw_c8(ia,ja),&
         vdw_r0(itypat,jtypat)
         call wrtout(std_out,msg,'COLL')
       end do
     end do
   end if
   write(msg, '(3a,f6.3,a,f6.3)') &
&   '  ---------------------------------------------------------------',ch10,&
&   '      Scaling factors:       s6 = ', vdw_s6,',    s8 = ',vdw_s8
   call wrtout(std_out,msg,'COLL')
   if (vdw_xc==6) then
     write(msg,'(a,f6.3,a,f6.3)') &
&     '      Damping parameters:   sr6 = ', vdw_sr6,',   sr8 = ',vdw_sr8
     call wrtout(std_out,msg,'COLL')
   elseif (vdw_xc==7) then
     write(msg,'(a,f6.3,a,f6.3)') &
&     '      Damping parameters:    a1 = ', vdw_a1, ',    a2 = ', vdw_a2
     call wrtout(std_out,msg,'COLL')
   end if
   write(msg,'(a,es12.5,3a,i14,2a,es12.5,1a)') &
&   '      Cut-off radius   = ',rcut,' Bohr',ch10,&
&   '      Number of pairs contributing = ',npairs,ch10,&
&   '      DFT-D3 (no 3-body) energy contribution = ',e_vdw_dftd3-e_3bt,' Ha'
   call wrtout(std_out,msg,'COLL')
   if (bol_3bt) then
     write(msg,'(6a,i5,2a,es20.11,3a,es20.11,1a)')ch10,&
&     '  ---------------------------------------------------------------',ch10,&
&     '      3-Body Term Contribution:', ch10,&
&     '      Number of shells considered    = ', nshell, ch10,&
&     '      Additional 3-body contribution = ', e_3bt, ' Ha',ch10,&
&     '      Total E (2-body and 3-body)    = ', e_vdw_dftd3, 'Ha'
     call wrtout(std_out,msg,'COLL')
   end if
   write(msg,'(2a)')&
&   '  ----------------------------------------------------------------',ch10
   call wrtout(std_out,msg,'COLL')
 end if
 ABI_FREE(ivdw)
 ABI_FREE(xred01)
 ABI_FREE(vdw_r0)
 ABI_FREE(fe_no_c)
 ABI_FREE(e_no_c)
 ABI_FREE(e_alpha1)
 ABI_FREE(e_alpha2)
 ABI_FREE(e_alpha3)
 ABI_FREE(e_alpha4)
 ABI_FREE(grad_no_cij)
 ABI_FREE(fgrad_no_c)
 ABI_FREE(cfgrad_no_c)
 ABI_FREE(vdw_c6)
 ABI_FREE(vdw_c8)
 ABI_FREE(dc6ri)
 ABI_FREE(dc6rj)
 ABI_FREE(d2c6ri)
 ABI_FREE(d2c6rj)
 ABI_FREE(d2c6rirj)
 ABI_FREE(cn)
 ABI_FREE(d2cn_iii)
 ABI_FREE(d2cn_jji)
 ABI_FREE(d2cn_iji)
 ABI_FREE(d2cn_jii)
 ABI_FREE(d2cn_tmp)
 ABI_FREE(dcn)
 ABI_FREE(fdcn)
 ABI_FREE(cfdcn)
 ABI_FREE(str_dcn)
 ABI_FREE(elt_cn)
 ABI_FREE(str_no_c)
 ABI_FREE(str_alpha1)
 ABI_FREE(str_alpha2)
 ABI_FREE(dcn_cart)
 DBG_EXIT("COLL")

 contains
!! ***

!!****f*vdw_dftd3/comp_prod
!!
!! NAME
!! comp_prod
!!
!! FUNCTION
!! Return the product of two complex numbers stored in rank 1 array
!!
!! SOURCE

   subroutine comp_prod(a,b,c)

   implicit none
 !Arguments ----------------------
   real(dp),intent(in) :: a(2),b(2)
   real(dp),intent(out) :: c(2)

! *********************************************************************

   c(1) = a(1)*b(1)-a(2)*b(2)
   c(2) = a(1)*b(2)+a(2)*b(1)

 end subroutine comp_prod
!!***

!!****f*vdw_dftd3/d3_cart2red
!!
!! NAME
!! d3_cart2red
!!
!! FUNCTION
!! Convert gradients from cartesian to reduced coordinates
!!
!! SOURCE

subroutine d3_cart2red(grad)

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

 end subroutine d3_cart2red
!!***

end subroutine vdw_dftd3
!!***

end module m_vdw_dftd3
!!***
