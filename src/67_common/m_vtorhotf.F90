!!****m* ABINIT/m_vtorhotf
!! NAME
!!  m_vtorhotf
!!
!! FUNCTION
!! Computes the new density from a fixed potential (vtrial) using the Thomas-Fermi functional
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, MF, AR, MM)
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

module m_vtorhotf

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_dtset

 use defs_abitypes, only : MPI_type
 use m_time,     only : timab
 use m_spacepar,  only : symrhg

 implicit none

 private
!!***

 public :: vtorhotf
!!***

contains
!!***

!!****f* ABINIT/vtorhotf
!! NAME
!! vtorhotf
!!
!! FUNCTION
!! This routine computes the new density from a fixed potential (vtrial)
!! using the Thomas-Fermi functional
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=number of fft grid points
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetry elements in space group
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  gprimd(3,3)=dimensional real space primitive translations
!!  ucvol=unit cell volume in bohr**3.
!!  vtrial(nfft,nspden)=INPUT Vtrial(r).
!!
!! OUTPUT
!!  ek=kinetic energy part of total energy.
!!  enlx=nonlocal psp + potential Fock ACE part of total energy.
!!  entropy=entropy due to the occupation number smearing (if metal)
!!  fermie=fermi energy (Hartree)
!!  grnl(3*natom)=stores grads of nonlocal energy wrt length scales
!!   (3x3 tensor) and grads wrt atomic coordinates (3*natom)
!!
!! SIDE EFFECTS
!!  rhog(2,nfft)=array for Fourier transform of electron density
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!
!! PARENTS
!!      m_scfcv_core
!!
!! CHILDREN
!!
!! SOURCE

subroutine vtorhotf(dtset,ek,enlx,entropy,fermie,gprimd,grnl,&
&  irrzon,mpi_enreg,natom,nfft,nspden,nsppol,nsym,phnons,rhog,rhor,rprimd,ucvol,vtrial)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,nsppol,nsym
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: ek,enlx,entropy,fermie
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: irrzon((dtset%ngfft(1)*dtset%ngfft(1)*dtset%ngfft(1))**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: phnons(2,(dtset%ngfft(1)*dtset%ngfft(1)*dtset%ngfft(1))**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))
 real(dp),intent(in) :: rprimd(3,3),vtrial(nfft,nspden)
 real(dp),intent(inout) :: rhog(2,nfft),rhor(nfft,nspden)
 real(dp),intent(out) :: grnl(3*natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: jdichomax=20,level=111
 integer :: i1,i2,i3,ierr,ifft,ii,ir,iscf,jdicho
 integer :: me_fft,n1,n2,n3,nfftot,nproc_fft,prtvol
 real(dp),save :: cktf,fermie_tol,nelect_mid
 real(dp) :: dnelect_mid_dx,dxrtnewt,eektemp,eektf,feektemp,feektf
 real(dp) :: rtnewt,sum_rhor_mid,sum_rhor_middx
 logical,save :: lfirst_time_tf=.true.
 logical :: lnewtonraphson
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: betamumoinsV(:),rhor_mid(:),rhor_middx(:)

! *************************************************************************

!Keep track of total time spent in vtorho
 call timab(21,1,tsec)

!Structured debugging if prtvol==-level
 prtvol=dtset%prtvol
 if(prtvol==-level)then
   write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' vtorho : enter '
   call wrtout(std_out,message,'COLL')
 end if

 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 me_fft=dtset%ngfft(11) ; nproc_fft=dtset%ngfft(10)
 iscf=dtset%iscf
!Debugging : print vtrial and rhor
 if(prtvol==-level)then
   write(message,'(a)') '   ir              vtrial(ir)     rhor(ir) '
   call wrtout(std_out,message,'COLL')
   do ir=1,nfft
!    if(ir<=11 .or. mod(ir,301)==0 )then
     i3=(ir-1)/n1/(n2/nproc_fft)
     i2=(ir-1-i3*n1*n2/nproc_fft)/n1
     i1=ir-1-i3*n1*n2/nproc_fft-(i2-me_fft)*n1
     write(message,'(i5,3i3,a,2d13.6)')ir,i1,i2,i3,' ',vtrial(ir,1),rhor(ir,1)
     call wrtout(std_out,message,'COLL')
     if(nspden>=2)then
       write(message,'(a,2d13.6)')'               ',vtrial(ir,2),rhor(ir,2)
       call wrtout(std_out,message,'COLL')
     end if
!    end if
   end do
 end if

 ek=zero
 enlx=zero
 grnl(:)=zero

!Initialize rhor if needed
 if(iscf>0) rhor(:,:)=zero

!call Thomas-Fermi for the density
 call tf
!Compute energy terms
 call tfek

 call timab(21,2,tsec)
!End thomas fermi

 contains
!!***

!!****f* vtorhotf/tf
!! NAME
!! tf
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_vtorhotf
!!
!! CHILDREN
!!
!! SOURCE
  subroutine tf()

! *************************************************************************

   ABI_MALLOC(rhor_mid,(nfft))
   ABI_MALLOC(rhor_middx,(nfft))
   fermie_tol=1.e-10_dp
   cktf=one/two/pi**2*(two*dtset%tphysel)**1.5_dp

!  Should be made an input variable, if TF really needed for production
!  rtnewt=dtset%userra
   rtnewt=zero

!  Newton Raphson
   if (lfirst_time_tf) then
     lfirst_time_tf=.false.
   end if
   jdicho=0
   lnewtonraphson=.false.
   do while (.not.lnewtonraphson)
     jdicho=jdicho+1
!    do ifft=1,nfft
!    rhor_mid(ifft)=cktf*zfermi12((rtnewt-vtrial(ifft,1))/dtset%tphysel)
!    rhor_middx(ifft)=cktf*zfermim12((rtnewt-vtrial(ifft,1))/dtset%tphysel)
!    end do
     call fm12a1t(cktf,rtnewt,dtset%tphysel,vtrial(:,1),rhor_middx,rhor_mid,&
&     nfft)
     sum_rhor_mid=sum(rhor_mid(:))
     sum_rhor_middx=sum(rhor_middx(:))
     call xmpi_sum(sum_rhor_mid,mpi_enreg%comm_fft ,ierr)
     call xmpi_sum(sum_rhor_middx,mpi_enreg%comm_fft ,ierr)
     nelect_mid=sum_rhor_mid*ucvol/(nfft*nproc_fft)-dtset%nelect
     dnelect_mid_dx=sum_rhor_middx*ucvol/(nfft*nproc_fft)/dtset%tphysel/two
     dxrtnewt=nelect_mid/dnelect_mid_dx
     rtnewt=rtnewt-dxrtnewt
     if (abs(nelect_mid) < fermie_tol/2._dp) then
       lnewtonraphson=.true.
     end if
     if (jdicho > jdichomax) then
       ABI_ERROR('NEWTON RAPHSON NOT CONVERGED')
     end if
   end do
   fermie=rtnewt
   rhor(:,1)=rhor_mid(:)
   ABI_FREE(rhor_mid)
   ABI_FREE(rhor_middx)

!  DEBUG
!  write(std_out,*)'fmid,nmid,jdicho',fermie,nelect_mid,jdicho
!  ENDDEBUG

!  Compute rhog
   call timab(70,1,tsec)

   nfftot=dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)
   call symrhg(1,gprimd,irrzon,mpi_enreg,nfft,nfftot,dtset%ngfft,nspden,nsppol,nsym,phnons,&
&   rhog,rhor,rprimd,dtset%symafm,dtset%symrel,dtset%tnons)

!  We now have both rho(r) and rho(G), symmetrized, and if nsppol=2
!  we also have the spin-up density, symmetrized, in rhor(:,2).

   call timab(70,2,tsec)

end subroutine tf
!!***

!!****f* vtorhotf/tfek
!! NAME
!! tfek
!!
!! FUNCTION
!! This is the calculation of the kinetic energy for Thomas Fermi
!! Energy and free energy must be distinguished
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_vtorhotf
!!
!! CHILDREN
!!
!! SOURCE

  subroutine tfek()

! *************************************************************************

   ABI_MALLOC(betamumoinsV,(nfft))
   cktf=one/two/pi**2*(two*dtset%tphysel)**1.5_dp
   eektf=zero
   feektf=zero
   do ifft=1,nfft

!    betamumoinsv(ifft)=ifermi12(rhor(ifft,1)/cktf)
     betamumoinsv(ifft)=(rtnewt-vtrial(ifft,1))/dtset%tphysel
!    eektemp=zfermi32(betamumoinsV(ifft))/zfermi12(betamumoinsV(ifft))
     eektemp=fp32a1(betamumoinsV(ifft))/rhor(ifft,1)*cktf
     feektemp=betamumoinsV(ifft)-two/three*eektemp
     feektf=feektf+feektemp*rhor(ifft,1)
     eektf=eektf+eektemp*rhor(ifft,1)
   end do
!  Init mpi_comm
   call timab(48,1,tsec)
   call xmpi_sum(eektf,mpi_enreg%comm_fft ,ierr)
   call xmpi_sum(feektf,mpi_enreg%comm_fft ,ierr)
   call timab(48,2,tsec)
   eektf=eektf*dtset%tphysel
   eektf=eektf*ucvol/dble(nfft*nproc_fft)
   feektf=feektf*dtset%tphysel
   feektf=feektf*ucvol/dble(nfft*nproc_fft)
!  DEBUG
!  write(std_out,*)'eektf',eektf
!  stop ('vtorhotf')
!  ENDDEBUG
   ek=eektf
   entropy=(eektf-feektf)/dtset%tphysel
   ABI_FREE(betamumoinsV)
 end subroutine tfek
!!***

!!****f* ABINIT/zfermim12
!! NAME
!! zfermim12
!!
!! FUNCTION
!!..file contains fermi-dirac integral routines:
!!..
!!..function zfermim12 does a rational function fit for the order -1/2 integral
!!..function zfermi12 does a rational function fit for the order 1/2 integral
!!..function zfermi1 does a rational function fit for the order 1 integral
!!..function zfermi32 does a rational function fit for the order 3/2 integral
!!..function zfermi2 does a rational function fit for the order 2 integral
!!..function zfermi52 does a rational function fit for the order 5/2 integral
!!..function zfermi3 does a rational function fit for the order 3 integral
!!..
!!..function ifermim12 is a rational function fit for the inverse of order -1/2
!!..function ifermi12 is a rational function fit for the inverse of order 1/2
!!..function ifermi32 is a rational function fit for the inverse of order 3/2
!!..function ifermi52 is a rational function fit for the inverse of order 5/2
!!
!!..this routine applies a rational function expansion to get the fermi-dirac
!!..integral of order -1/2 evaluated at x. maximum error is 1.23d-12.
!!..reference: antia apjs 84,101 1993
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function zfermim12(xx)

!Arguments -------------------------------
 real(dp), intent(in) :: xx
 real(dp) :: zfermim12

!Local variables-------------------------------
 integer ::  ii,m1,k1,m2,k2
 real(dp) :: an,a1(12),b1(12),a2(12),b2(12),rn,den,xx1
!..load the coefficients of the expansion
 data  an,m1,k1,m2,k2 /-0.5e0_dp, 7, 7, 11, 11/
 data  (a1(ii),ii=1,8)/ 1.71446374704454e7_dp,    3.88148302324068e7_dp,&
& 3.16743385304962e7_dp,    1.14587609192151e7_dp,&
& 1.83696370756153E6_dp,    1.14980998186874e5_dp,&
& 1.98276889924768e3_dp,    1.0e0_dp/
 data  (b1(ii),ii=1,8)/ 9.67282587452899e6_dp,    2.87386436731785e7_dp,&
& 3.26070130734158e7_dp,    1.77657027846367e7_dp,&
& 4.81648022267831e6_dp,    6.13709569333207e5_dp,&
& 3.13595854332114e4_dp,    4.35061725080755e2_dp/
 data (a2(ii),ii=1,12)/-4.46620341924942e-15_dp, -1.58654991146236e-12_dp,&
& -4.44467627042232e-10_dp, -6.84738791621745e-8_dp,&
& -6.64932238528105e-6_dp,  -3.69976170193942e-4_dp,&
& -1.12295393687006e-2_dp,  -1.60926102124442e-1_dp,&
& -8.52408612877447e-1_dp,  -7.45519953763928e-1_dp,&
& 2.98435207466372e0_dp,    1.0e0_dp/
 data (b2(ii),ii=1,12)/-2.23310170962369e-15_dp, -7.94193282071464e-13_dp,&
& -2.22564376956228e-10_dp, -3.43299431079845e-8_dp,&
& -3.33919612678907e-6_dp,  -1.86432212187088e-4_dp,&
& -5.69764436880529e-3_dp,  -8.34904593067194e-2_dp,&
& -4.78770844009440e-1_dp,  -4.99759250374148e-1_dp,&
& 1.86795964993052e0_dp,    4.16485970495288e-1_dp/

! *************************************************************************

 if (xx .lt. 2.0e0_dp) then
   xx1 = exp(xx)
   rn = xx1 + a1(m1)
   do ii=m1-1,1,-1
     rn = rn*xx1 + a1(ii)
   end do
   den = b1(k1+1)
   do ii=k1,1,-1
     den = den*xx1 + b1(ii)
   end do
   zfermim12 = xx1 * rn/den
!  ..
 else
   xx1 = one/(xx*xx)
   rn = xx1 + a2(m2)
   do ii=m2-1,1,-1
     rn = rn*xx1 + a2(ii)
   end do
   den = b2(k2+1)
   do ii=k2,1,-1
     den = den*xx1 + b2(ii)
   end do
   zfermim12 = sqrt(xx)*rn/den
 end if

end function zfermim12
!!***

!!****f* ABINIT/zfermi12
!! NAME
!! zfermi12
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function zfermi12(xx)
!..
!..this routine applies a rational function expansion to get the fermi-dirac
!..integral of order 1/2 evaluated at x. maximum error is 5.47d-13.
!..reference: antia apjs 84,101 1993
!..
!..declare

!Arguments -------------------------------
 real(dp), intent(in) :: xx
 real(dp):: zfermi12

!Local variables-------------------------------
 integer ::         ii,m1,k1,m2,k2
 real(dp) :: an,a1(12),b1(12),a2(12),b2(12),rn,den,xx1

!..load the coefficients of the expansion
 data  an,m1,k1,m2,k2 /0.5e0_dp, 7, 7, 10, 11/
 data  (a1(ii),ii=1,8)/5.75834152995465e6_dp,   1.30964880355883e7_dp,&
& 1.07608632249013e7_dp,   3.93536421893014e6_dp,&
& 6.42493233715640e5_dp,   4.16031909245777e4_dp,&
& 7.77238678539648e2_dp,   1.0e0_dp/
 data  (b1(ii),ii=1,8)/6.49759261942269e6_dp,   1.70750501625775e7_dp,&
& 1.69288134856160e7_dp,   7.95192647756086e6_dp,&
& 1.83167424554505e6_dp,   1.95155948326832e5_dp,&
& 8.17922106644547e3_dp,   9.02129136642157e1_dp/
 data (a2(ii),ii=1,11)/4.85378381173415e-14_dp, 1.64429113030738e-11_dp,&
& 3.76794942277806e-9_dp,  4.69233883900644e-7_dp,&
& 3.40679845803144e-5_dp,  1.32212995937796e-3_dp,&
& 2.60768398973913e-2_dp,  2.48653216266227e-1_dp,&
& 1.08037861921488e0_dp,   1.91247528779676e0_dp,&
& 1.0e0_dp/
 data (b2(ii),ii=1,12)/7.28067571760518e-14_dp, 2.45745452167585e-11_dp,&
& 5.62152894375277e-9_dp,  6.96888634549649e-7_dp,&
& 5.02360015186394e-5_dp,  1.92040136756592e-3_dp,&
& 3.66887808002874e-2_dp,  3.24095226486468e-1_dp,&
& 1.16434871200131e0_dp,   1.34981244060549e0_dp,&
& 2.01311836975930e-1_dp, -2.14562434782759e-2_dp/

! *************************************************************************

 if (xx .lt. two) then
   xx1 = exp(xx)
   rn = xx1 + a1(m1)
   do ii=m1-1,1,-1
     rn = rn*xx1 + a1(ii)
   end do
   den = b1(k1+1)
   do ii=k1,1,-1
     den = den*xx1 + b1(ii)
   end do
   zfermi12 = xx1 * rn/den

 else
   xx1 = one/(xx*xx)
   rn = xx1 + a2(m2)
   do ii=m2-1,1,-1
     rn = rn*xx1 + a2(ii)
   end do
   den = b2(k2+1)
   do ii=k2,1,-1
     den = den*xx1 + b2(ii)
   end do
   zfermi12 = xx*sqrt(xx)*rn/den
 end if

end function zfermi12
!!***

!!****f* ABINIT/zfermi1
!! NAME
!! zfermi1
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function zfermi1(xx)
!..
!..this routine applies a rational function expansion to get the fermi-dirac
!..integral of order 1 evaluated at x. maximum error is 1.0e-8.
!..reference: antia  priv comm. 11sep94
!..
!..declare

!Arguments -------------------------------
 real(dp), intent(in) :: xx
 real(dp):: zfermi1
!Local variables-------------------------------
 integer ::  ii,m1,k1,m2,k2
 real(dp) :: an,a1(12),b1(12),a2(12),b2(12),rn,den,xx1

!..load the coefficients of the expansion
 data  an,m1,k1,m2,k2 /1.0_dp, 7, 4, 9, 5/
 data  (a1(ii),ii=1,8)/-7.606458638543e7_dp,  -1.143519707857e8_dp,&
& -5.167289383236e7_dp,  -7.304766495775e6_dp,&
& -1.630563622280e5_dp,   3.145920924780e3_dp,&
& -7.156354090495e1_dp,   1.0_dp/
 data  (b1(ii),ii=1,5)/-7.606458639561e7_dp,  -1.333681162517e8_dp,&
& -7.656332234147e7_dp,  -1.638081306504e7_dp,&
& -1.044683266663e6_dp/
 data (a2(ii),ii=1,10)/-3.493105157219e-7_dp, -5.628286279892e-5_dp,&
& -5.188757767899e-3_dp, -2.097205947730e-1_dp,&
& -3.353243201574_dp,    -1.682094530855e1_dp,&
& -2.042542575231e1_dp,   3.551366939795_dp,&
& -2.400826804233_dp,     1.0_dp/
 data  (b2(ii),ii=1,6)/-6.986210315105e-7_dp, -1.102673536040e-4_dp,&
& -1.001475250797e-2_dp, -3.864923270059e-1_dp,&
& -5.435619477378_dp,    -1.563274262745e1_dp/

! *************************************************************************

 if (xx .lt. 2.0_dp) then
   xx1 = exp(xx)
   rn = xx1 + a1(m1)
   do ii=m1-1,1,-1
     rn = rn*xx1 + a1(ii)
   end do
   den = b1(k1+1)
   do ii=k1,1,-1
     den = den*xx1 + b1(ii)
   end do
   zfermi1 = xx1 * rn/den

 else
   xx1 = 1.0_dp/(xx*xx)
   rn = xx1 + a2(m2)
   do ii=m2-1,1,-1
     rn = rn*xx1 + a2(ii)
   end do
   den = b2(k2+1)
   do ii=k2,1,-1
     den = den*xx1 + b2(ii)
   end do
   zfermi1 = xx*xx*rn/den
 end if

end function zfermi1
!!***

!!****f* ABINIT/zfermi32
!! NAME
!! zfermi32
!!
!! FUNCTION
!!  this routine applies a rational function expansion to get the fermi-dirac
!!  integral of order 3/2 evaluated at x. maximum error is 5.07d-13.
!!  reference: antia apjs 84,101 1993
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function zfermi32(xx)

!Arguments -------------------------------
 real(dp), intent(in) :: xx
 real(dp) :: zfermi32

!Local variables-------------------------------
 integer :: ii,m1,k1,m2,k2
 real(dp) :: an,a1(12),b1(12),a2(12),b2(12),rn,den,xx1

!..load the coefficients of the expansion
 data  an,m1,k1,m2,k2 /1.5e0_dp, 6, 7, 9, 10/
 data  (a1(ii),ii=1,7)/4.32326386604283e4_dp,   8.55472308218786e4_dp,&
& 5.95275291210962e4_dp,   1.77294861572005e4_dp,&
& 2.21876607796460e3_dp,   9.90562948053193e1_dp,&
& 1.0e0_dp/
 data  (b1(ii),ii=1,8)/3.25218725353467e4_dp,   7.01022511904373e4_dp,&
& 5.50859144223638e4_dp,   1.95942074576400e4_dp,&
& 3.20803912586318e3_dp,   2.20853967067789e2_dp,&
& 5.05580641737527e0_dp,   1.99507945223266e-2_dp/
 data (a2(ii),ii=1,10)/2.80452693148553e-13_dp, 8.60096863656367e-11_dp,&
& 1.62974620742993e-8_dp,  1.63598843752050e-6_dp,&
& 9.12915407846722e-5_dp,  2.62988766922117e-3_dp,&
& 3.85682997219346e-2_dp,  2.78383256609605e-1_dp,&
& 9.02250179334496e-1_dp,  1.0e0_dp/
 data (b2(ii),ii=1,11)/7.01131732871184e-13_dp, 2.10699282897576e-10_dp,&
& 3.94452010378723e-8_dp,  3.84703231868724e-6_dp,&
& 2.04569943213216e-4_dp,  5.31999109566385e-3_dp,&
& 6.39899717779153e-2_dp,  3.14236143831882e-1_dp,&
& 4.70252591891375e-1_dp, -2.15540156936373e-2_dp,&
& 2.34829436438087e-3_dp/

! *************************************************************************

 if (xx .lt. 2.0e0_dp) then
   xx1 = exp(xx)
   rn = xx1 + a1(m1)
   do ii=m1-1,1,-1
     rn = rn*xx1 + a1(ii)
   end do
   den = b1(k1+1)
   do ii=k1,1,-1
     den = den*xx1 + b1(ii)
   end do
   zfermi32 = xx1 * rn/den

 else
   xx1 = one/(xx*xx)
   rn = xx1 + a2(m2)
   do ii=m2-1,1,-1
     rn = rn*xx1 + a2(ii)
   end do
   den = b2(k2+1)
   do ii=k2,1,-1
     den = den*xx1 + b2(ii)
   end do
   zfermi32 = xx*xx*sqrt(xx)*rn/den
 end if

end function zfermi32
!!***

!!****f* ABINIT/zfermi2
!! NAME
!! zfermi2
!!
!! FUNCTION
!!
!!  this routine applies a rational function expansion to get the fermi-dirac
!!  integral of order 2 evaluated at x. maximum error is 1.0e-8.
!!  reference: antia  priv comm. 11sep94
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function zfermi2(xx)

!Arguments -------------------------------
 real(dp), intent(in) :: xx
 real(dp) :: zfermi2
!Local variables-------------------------------
 integer ::  ii,m1,k1,m2,k2
 real(dp) :: an,a1(12),b1(12),a2(12),b2(12),rn,den,xx1

!..load the coefficients of the expansion
 data  an,m1,k1,m2,k2 /2.0_dp, 7, 4, 5, 9/
 data  (a1(ii),ii=1,8)/-1.434885992395e8_dp,  -2.001711155617e8_dp,&
& -8.507067153428e7_dp,  -1.175118281976e7_dp,&
& -3.145120854293e5_dp,   4.275771034579e3_dp,&
& -8.069902926891e1_dp,   1.0e0_dp/
 data  (b1(ii),ii=1,5)/-7.174429962316e7_dp,  -1.090535948744e8_dp,&
& -5.350984486022e7_dp,  -9.646265123816e6_dp,&
& -5.113415562845e5_dp/
 data  (a2(ii),ii=1,6)/ 6.919705180051e-8_dp,  1.134026972699e-5_dp,&
& 7.967092675369e-4_dp,  2.432500578301e-2_dp,&
& 2.784751844942e-1_dp,  1.0e0_dp/
 data (b2(ii),ii=1,10)/ 2.075911553728e-7_dp,  3.197196691324e-5_dp,&
& 2.074576609543e-3_dp,  5.250009686722e-2_dp,&
& 3.171705130118e-1_dp, -1.147237720706e-1_dp,&
& 6.638430718056e-2_dp, -1.356814647640e-2_dp,&
& -3.648576227388e-2_dp,  3.621098757460e-2_dp/

! *************************************************************************

 if (xx .lt. 2.0e0_dp) then
   xx1 = exp(xx)
   rn = xx1 + a1(m1)
   do ii=m1-1,1,-1
     rn = rn*xx1 + a1(ii)
   end do
   den = b1(k1+1)
   do ii=k1,1,-1
     den = den*xx1 + b1(ii)
   end do
   zfermi2 = xx1 * rn/den

 else
   xx1 = one/(xx*xx)
   rn = xx1 + a2(m2)
   do ii=m2-1,1,-1
     rn = rn*xx1 + a2(ii)
   end do
   den = b2(k2+1)
   do ii=k2,1,-1
     den = den*xx1 + b2(ii)
   end do
   zfermi2 = xx*xx*xx*rn/den
 end if

end function zfermi2
!!***

!!****f* ABINIT/zfermi52
!! NAME
!! zfermi52
!!
!! FUNCTION
!!  this routine applies a rational function expansion to get the fermi-dirac
!!  integral of order 5/2 evaluated at x. maximum error is 2.47d-13.
!!  reference: antia apjs 84,101 1993
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function zfermi52(xx)

!Arguments -------------------------------
 real(dp), intent(in) :: xx
 real(dp) :: zfermi52

!Local variables-------------------------------
 integer :: ii,m1,k1,m2,k2
 real(dp) :: an,a1(12),b1(12),a2(12),b2(12),rn,den,xx1

!..load the coefficients of the expansion
 data  an,m1,k1,m2,k2 /2.5e0_dp, 6, 7, 10, 9/
 data  (a1(ii),ii=1,7)/6.61606300631656e4_dp,   1.20132462801652e5_dp,&
& 7.67255995316812e4_dp,   2.10427138842443e4_dp,&
& 2.44325236813275e3_dp,   1.02589947781696e2_dp,&
& 1.0e0_dp/
 data  (b1(ii),ii=1,8)/1.99078071053871e4_dp,   3.79076097261066e4_dp,&
& 2.60117136841197e4_dp,   7.97584657659364e3_dp,&
& 1.10886130159658e3_dp,   6.35483623268093e1_dp,&
& 1.16951072617142e0_dp,   3.31482978240026e-3_dp/
 data (a2(ii),ii=1,11)/8.42667076131315e-12_dp, 2.31618876821567e-9_dp,&
& 3.54323824923987e-7_dp,  2.77981736000034e-5_dp,&
& 1.14008027400645e-3_dp,  2.32779790773633e-2_dp,&
& 2.39564845938301e-1_dp,  1.24415366126179e0_dp,&
& 3.18831203950106e0_dp,   3.42040216997894e0_dp,&
& 1.0e0_dp/
 data (b2(ii),ii=1,10)/2.94933476646033e-11_dp, 7.68215783076936e-9_dp,&
& 1.12919616415947e-6_dp,  8.09451165406274e-5_dp,&
& 2.81111224925648e-3_dp,  3.99937801931919e-2_dp,&
& 2.27132567866839e-1_dp,  5.31886045222680e-1_dp,&
& 3.70866321410385e-1_dp,  2.27326643192516e-2_dp/

! *************************************************************************

 if (xx .lt. two) then
   xx1 = exp(xx)
   rn = xx1 + a1(m1)
   do ii=m1-1,1,-1
     rn = rn*xx1 + a1(ii)
   end do
   den = b1(k1+1)
   do ii=k1,1,-1
     den = den*xx1 + b1(ii)
   end do
   zfermi52 = xx1 * rn/den

 else
   xx1 = one/(xx*xx)
   rn = xx1 + a2(m2)
   do ii=m2-1,1,-1
     rn = rn*xx1 + a2(ii)
   end do
   den = b2(k2+1)
   do ii=k2,1,-1
     den = den*xx1 + b2(ii)
   end do
   zfermi52 = xx*xx*xx*sqrt(xx)*rn/den
 end if

end function zfermi52
!!***

!!****f* ABINIT/zfermi3
!! NAME
!! zfermi3
!!
!! FUNCTION
!!  this routine applies a rational function expansion to get the fermi-dirac
!!  integral of order 3 evaluated at x. maximum error is 1.0e-8.
!!  reference: antia  priv comm. 11sep94
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function zfermi3(xx)

!Arguments -------------------------------
 real(dp), intent(in) :: xx
 real(dp):: zfermi3

!Local variables-------------------------------
 integer :: ii,m1,k1,m2,k2
 real(dp) :: an,a1(12),b1(12),a2(12),b2(12),rn,den,xx1

!..load the coefficients of the expansion
 data  an,m1,k1,m2,k2 /3.0, 4, 6, 7, 7/
 data  (a1(ii),ii=1,5)/ 6.317036716422e2_dp,    7.514163924637e2_dp,&
& 2.711961035750e2_dp,    3.274540902317e1_dp,&
& 1.0_dp/
 data  (b1(ii),ii=1,7)/ 1.052839452797e2_dp,    1.318163114785e2_dp,&
& 5.213807524405e1_dp,    7.500064111991_dp,&
& 3.383020205492e-1_dp,   2.342176749453e-3_dp,&
& -8.445226098359e-6_dp/
 data  (a2(ii),ii=1,8)/ 1.360999428425e-8_dp,   1.651419468084e-6_dp,&
& 1.021455604288e-4_dp,   3.041270709839e-3_dp,&
& 4.584298418374e-2_dp,   3.440523212512e-1_dp,&
& 1.077505444383_dp,    1.0_dp/
 data  (b2(ii),ii=1,8)/ 5.443997714076e-8_dp,   5.531075760054e-6_dp,&
& 2.969285281294e-4_dp,   6.052488134435e-3_dp,&
& 5.041144894964e-2_dp,   1.048282487684e-1_dp,&
& 1.280969214096e-2_dp,  -2.851555446444e-3_dp/

! *************************************************************************

 if (xx .lt. two) then
   xx1 = exp(xx)
   rn = xx1 + a1(m1)
   do ii=m1-1,1,-1
     rn = rn*xx1 + a1(ii)
   end do
   den = b1(k1+1)
   do ii=k1,1,-1
     den = den*xx1 + b1(ii)
   end do
   zfermi3 = xx1 * rn/den

 else
   xx1 = one/(xx*xx)
   rn = xx1 + a2(m2)
   do ii=m2-1,1,-1
     rn = rn*xx1 + a2(ii)
   end do
   den = b2(k2+1)
   do ii=k2,1,-1
     den = den*xx1 + b2(ii)
   end do
   zfermi3 = xx*xx*xx*xx*rn/den
 end if

end function zfermi3
!!***

!!****f* ABINIT/ifermim12
!! NAME
!! ifermim12
!!
!! FUNCTION
!!  this routine applies a rational function expansion to get the inverse
!!  fermi-dirac integral of order -1/2 when it is equal to f.
!!  maximum error is 3.03d-9.   reference: antia apjs 84,101 1993
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function ifermim12(ff)

!Arguments -------------------------------
 real(dp), intent(in) :: ff
 real(dp) :: ifermim12
!Local variables-------------------------------
 integer :: ii,m1,k1,m2,k2
 real(dp) :: an,a1(12),b1(12),a2(12),b2(12),rn,den,ff1

!..load the coefficients of the expansion
 data  an,m1,k1,m2,k2 /-0.5e0_dp, 5, 6, 6, 6/
 data  (a1(ii),ii=1,6)/-1.570044577033e4_dp,   1.001958278442e4_dp,&
& -2.805343454951e3_dp,   4.121170498099e2_dp,&
& -3.174780572961e1_dp,   1.0e0_dp/
 data  (b1(ii),ii=1,7)/-2.782831558471e4_dp,   2.886114034012e4_dp,&
& -1.274243093149e4_dp,   3.063252215963e3_dp,&
& -4.225615045074e2_dp,   3.168918168284e1_dp,&
& -1.008561571363e0_dp/
 data  (a2(ii),ii=1,7)/ 2.206779160034e-8_dp,  -1.437701234283e-6_dp,&
& 6.103116850636e-5_dp,  -1.169411057416e-3_dp,&
& 1.814141021608e-2_dp,  -9.588603457639e-2_dp,&
& 1.0e0_dp/
 data  (b2(ii),ii=1,7)/ 8.827116613576e-8_dp,  -5.750804196059e-6_dp,&
& 2.429627688357e-4_dp,  -4.601959491394e-3_dp,&
& 6.932122275919e-2_dp,  -3.217372489776e-1_dp,&
& 3.124344749296e0_dp/

! *************************************************************************

 if (ff .lt. 4.0e0_dp) then
   rn = ff + a1(m1)
   do ii=m1-1,1,-1
     rn = rn*ff + a1(ii)
   end do
   den = b1(k1+1)
   do ii=k1,1,-1
     den = den*ff + b1(ii)
   end do
   ifermim12 = log(ff * rn/den)

 else
   ff1 = one/ff**(one/(one + an))
   rn = ff1 + a2(m2)
   do ii=m2-1,1,-1
     rn = rn*ff1 + a2(ii)
   end do
   den = b2(k2+1)
   do ii=k2,1,-1
     den = den*ff1 + b2(ii)
   end do
   ifermim12 = rn/(den*ff1)
 end if

end function ifermim12
!!***

!!****f* ABINIT/ifermi12
!! NAME
!! ifermi12
!!
!! FUNCTION
!!   this routine applies a rational function expansion to get the inverse
!!   fermi-dirac integral of order 1/2 when it is equal to f.
!!   maximum error is 4.19d-9.   reference: antia apjs 84,101 1993
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function ifermi12(ff)

!Arguments -------------------------------
 real(dp), intent(in) :: ff
 real(dp) :: ifermi12
!Local variables-------------------------------
 integer :: ii,m1,k1,m2,k2
 real(dp) :: an,a1(12),b1(12),a2(12),b2(12),rn,den,ff1

!..load the coefficients of the expansion
 data  an,m1,k1,m2,k2 /0.5e0_dp, 4, 3, 6, 5/
 data  (a1(ii),ii=1,5)/ 1.999266880833e4_dp,   5.702479099336e3_dp,&
& 6.610132843877e2_dp,   3.818838129486e1_dp,&
& 1.0e0_dp/
 data  (b1(ii),ii=1,4)/ 1.771804140488e4_dp,  -2.014785161019e3_dp,&
& 9.130355392717e1_dp,  -1.670718177489e0_dp/
 data  (a2(ii),ii=1,7)/-1.277060388085e-2_dp,  7.187946804945e-2_dp,&
& -4.262314235106e-1_dp,  4.997559426872e-1_dp,&
& -1.285579118012e0_dp,  -3.930805454272e-1_dp,&
& 1.0e0_dp/
 data  (b2(ii),ii=1,6)/-9.745794806288e-3_dp,  5.485432756838e-2_dp,&
& -3.299466243260e-1_dp,  4.077841975923e-1_dp,&
& -1.145531476975e0_dp,  -6.067091689181e-2_dp/

! *************************************************************************

 if (ff .lt. 4.0e0_dp) then
   rn = ff + a1(m1)
   do ii=m1-1,1,-1
     rn = rn*ff + a1(ii)
   end do
   den = b1(k1+1)
   do ii=k1,1,-1
     den = den*ff + b1(ii)
   end do
   ifermi12 = log(ff * rn/den)

 else
   ff1 = one/ff**(one/(one + an))
   rn = ff1 + a2(m2)
   do ii=m2-1,1,-1
     rn = rn*ff1 + a2(ii)
   end do
   den = b2(k2+1)
   do ii=k2,1,-1
     den = den*ff1 + b2(ii)
   end do
   ifermi12 = rn/(den*ff1)
 end if

end function ifermi12
!!***

!!****f* ABINIT/ifermi32
!! NAME
!! ifermi32
!!
!! FUNCTION
!!   this routine applies a rational function expansion to get the inverse
!!   fermi-dirac integral of order 3/2 when it is equal to f.
!!   maximum error is 2.26d-9.   reference: antia apjs 84,101 1993
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function ifermi32(ff)

!Arguments -------------------------------
 real(dp), intent(in) :: ff
 real(dp) :: ifermi32
!Local variables-------------------------------
 integer :: ii,m1,k1,m2,k2
 real(dp) :: an,a1(12),b1(12),a2(12),b2(12),rn,den,ff1

!..load the coefficients of the expansion
 data  an,m1,k1,m2,k2 /1.5e0_dp, 3, 4, 6, 5/
 data  (a1(ii),ii=1,4)/ 1.715627994191e2_dp,   1.125926232897e2_dp,&
& 2.056296753055e1_dp,   1.0e0_dp/
 data  (b1(ii),ii=1,5)/ 2.280653583157e2_dp,   1.193456203021e2_dp,&
& 1.167743113540e1_dp,  -3.226808804038e-1_dp,&
& 3.519268762788e-3_dp/
 data  (a2(ii),ii=1,7)/-6.321828169799e-3_dp, -2.183147266896e-2_dp,&
& -1.057562799320e-1_dp, -4.657944387545e-1_dp,&
& -5.951932864088e-1_dp,  3.684471177100e-1_dp,&
& 1.0e0_dp/
 data  (b2(ii),ii=1,6)/-4.381942605018e-3_dp, -1.513236504100e-2_dp,&
& -7.850001283886e-2_dp, -3.407561772612e-1_dp,&
& -5.074812565486e-1_dp, -1.387107009074e-1_dp/

! *************************************************************************

 if (ff .lt. 4.0e0_dp) then
   rn = ff + a1(m1)
   do ii=m1-1,1,-1
     rn = rn*ff + a1(ii)
   end do
   den = b1(k1+1)
   do ii=k1,1,-1
     den = den*ff + b1(ii)
   end do
   ifermi32 = log(ff * rn/den)

 else
   ff1 = one/ff**(one/(one + an))
   rn = ff1 + a2(m2)
   do ii=m2-1,1,-1
     rn = rn*ff1 + a2(ii)
   end do
   den = b2(k2+1)
   do ii=k2,1,-1
     den = den*ff1 + b2(ii)
   end do
   ifermi32 = rn/(den*ff1)
 end if

end function ifermi32
!!***

!!****f* ABINIT/ifermi52
!! NAME
!! ifermi52
!!
!! FUNCTION
!!   this routine applies a rational function expansion to get the inverse
!!   fermi-dirac integral of order 5/2 when it is equal to f.
!!   maximum error is 6.17d-9.   reference: antia apjs 84,101 1993
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function ifermi52(ff)

!Arguments -------------------------------
 real(dp), intent(in) :: ff
 real(dp) :: ifermi52

!Local variables-------------------------------
 integer :: ii,m1,k1,m2,k2
 real(dp) :: an,a1(12),b1(12),a2(12),b2(12),rn,den,ff1

!..load the coefficients of the expansion
 data  an,m1,k1,m2,k2 /2.5e0_dp, 2, 3, 6, 6/
 data  (a1(ii),ii=1,3)/ 2.138969250409e2_dp,   3.539903493971e1_dp,&
& 1.0e0_dp/
 data  (b1(ii),ii=1,4)/ 7.108545512710e2_dp,   9.873746988121e1_dp,&
& 1.067755522895e0_dp,  -1.182798726503e-2_dp/
 data  (a2(ii),ii=1,7)/-3.312041011227e-2_dp,  1.315763372315e-1_dp,&
& -4.820942898296e-1_dp,  5.099038074944e-1_dp,&
& 5.495613498630e-1_dp, -1.498867562255e0_dp,&
& 1.0e0_dp/
 data  (b2(ii),ii=1,7)/-2.315515517515e-2_dp,  9.198776585252e-2_dp,&
& -3.835879295548e-1_dp,  5.415026856351e-1_dp,&
& -3.847241692193e-1_dp,  3.739781456585e-2_dp,&
& -3.008504449098e-2_dp/

! *************************************************************************

 if (ff .lt. 4.0e0_dp) then
   rn = ff + a1(m1)
   do ii=m1-1,1,-1
     rn = rn*ff + a1(ii)
   end do
   den = b1(k1+1)
   do ii=k1,1,-1
     den = den*ff + b1(ii)
   end do
   ifermi52 = log(ff * rn/den)

 else
   ff1 = one/ff**(one/(one + an))
   rn = ff1 + a2(m2)
   do ii=m2-1,1,-1
     rn = rn*ff1 + a2(ii)
   end do
   den = b2(k2+1)
   do ii=k2,1,-1
     den = den*ff1 + b2(ii)
   end do
   ifermi52 = rn/(den*ff1)
 end if

end function ifermi52
!!***

!!****f* ABINIT/fp12a1
!! NAME
!! fp12a1
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function fp12a1 (x)

! Arguments -------------------------------
 real(dp),intent(in) :: x
 real(dp) :: fp12a1

 real(dp) :: y

!**********************************************************************
!*                                                                    *
!*               Integrale de Fermi d'ordre 1/2                       *
!*    Fp12(x) = somme de 0 a l'infini de (dt*t**1/2)/(1+exp(t-x))     *
!*                                                                    *
!**********************************************************************
!
!H. M. Antia, Astrophys. J. Suppl. 84, 101 (1993)
!Erreur relative maximum annoncee 5.54 e-5
!Erreur relative maximum constatee : -5.53e-5 pour eta = 2
!
 if (x.lt.2._dp) then
   y=exp(x)
   fp12a1=y*(21.8168_dp+y*(13.1693_dp+y))&
&   /(24.6180_dp+y*(23.5546_dp+y*(4.76290_dp+y*0.134481_dp)))
 else
   y=one/(x*x)
   fp12a1=x*sqrt(x)*(0.0473011_dp+y*(0.548433_dp+y))&
&   /(0.0709478_dp+y*(0.737041_dp+y*0.382065_dp))
 end if
!
!**********************************************************************
 end function fp12a1
!!***

!!****f* ABINIT/fp32a1
!! NAME
!! fp32a1
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function fp32a1 (x)

!Arguments -------------------------------
 real(dp),intent(in) :: x
 real(dp) :: fp32a1

 real(dp) :: y,x2

!
!**********************************************************************
!*                                                                    *
!*               Integrale de Fermi d'ordre 3/2                       *
!*    Fp32(x) = somme de 0 a l'infini de (dt*t**3/2)/(1+exp(t-x))     *
!*                                                                    *
!**********************************************************************
!
!H. M. Antia, Astrophys. J. Suppl. 84, 101 (1993)
!Erreur relative maximum annoncee 6.54 e-5
!Erreur relative maximum constatee : -5.84e-5 pour eta = -5
!
 if (x.lt.two) then
   y=exp(x)
   fp32a1=y*(135.863_dp+y*(49.2764_dp+y))/(102.210_dp+y*(55.0312_dp+y*4.23365_dp))
 else
   x2=x*x
   y=1._dp/x2
   fp32a1=x2*sqrt(x)*(0.154699_dp+y*(1.20037_dp+y))&
&   /(0.386765_dp+y*(0.608119_dp-y*0.165665_dp))
 end if
!
!**********************************************************************
 end function fp32a1
!!***

!!****f* ABINIT/xp12a1
!! NAME
!! xp12a1
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function xp12a1 (y)

!Arguments -------------------------------
 real(dp) :: xp12a1
 real(dp),intent(in) :: y

 real(dp),parameter :: deux=2._dp,deuxs3=deux/3._dp
 real(dp) :: top,den,z

!
!**********************************************************************
!*                                                                    *
!*              Calcul de eta tel que fp12 (eta) = y                  *
!*          ou fp12 est l'integrale de Fermi d'ordre +1/2             *
!*                                                                    *
!**********************************************************************
!
!H. M. Antia, Astrophys. J. Suppl. 84, 101 (1993)
!Erreur relative maximum annoncee sur exp(eta) : 3.02 e-5
!
 if (y.lt.4._dp) then
   top=44.593646_dp+y*(11.288764_dp+y)
   den=39.519346_dp+y*(-5.7517464_dp+y*0.26594291_dp)
   xp12a1=log(y*top/den)
 else
   z=y**(-deuxs3)
   top=34.873722_dp+z*(-26.922515_dp+z)
   den=26.612832_dp+z*(-20.452930_dp+z*11.808945_dp)
   xp12a1=top/(z*den)
 end if
!
!**********************************************************************
 end function xp12a1
!!***

!!****f* ABINIT/fm12a1
!! NAME
!! fm12a1
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

 function fm12a1 (x)

!Arguments -------------------------------
 real(dp),intent(in) :: x
 real(dp) :: fm12a1

 real(dp) :: y

!
!**********************************************************************
!*                                                                    *
!*               Integrale de Fermi d'ordre -1/2                      *
!*    Fm12(x) = somme de 0 a l'infini de (dt*t**-1/2)/(1+exp(t-x))    *
!*                                                                    *
!**********************************************************************
!
!H. M. Antia, Astrophys. J. Suppl. 84, 101 (1993)
!Erreur relative maximum annoncee 4.75 e-5
!
 if (x.lt.2._dp) then
   y=exp(x)
   fm12a1=y*(23.1456_dp+y*(13.7820_dp+y))&
&   /(13.0586_dp+y*(17.0048_dp+y*(5.07527_dp+y*0.236620_dp)))
 else
   y=1./(x*x)
   fm12a1=sqrt(x)*(0.0153602_dp+y*(0.146815_dp+y))&
&   /(0.00768015_dp+y*(0.0763700_dp+y*0.570485_dp))
 end if
!
!**********************************************************************
 end function fm12a1
!!***

!!****f* ABINIT/fm12a1t
!! NAME
!! fm12a1t
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_vtorhotf
!!
!! CHILDREN
!!
!! SOURCE

 subroutine fm12a1t (cktf,rtnewt,tphysel,vtrial,rhor_middx,rhor_mid,nfft)

 integer,intent(in) :: nfft
 real(dp),intent(in) :: tphysel,rtnewt,cktf
 real(dp),intent(in) :: vtrial(nfft)
 real(dp),intent(out) :: rhor_middx(nfft),rhor_mid(nfft)

 !intrinsic exp,sqrt
 integer :: ifft
 real(dp) :: x,y,sqrtx

!
!**********************************************************************
!*                                                                    *
!*               Integrale de Fermi d'ordre -1/2                      *
!*    Fm12(x) = somme de 0 a l'infini de (dt*t**-1/2)/(1+exp(t-x))    *
!*                      ....                                              *
!**********************************************************************
!
!H. M. Antia, Astrophys. J. Suppl. 84, 101 (1993)
!Erreur relative maximum annoncee 4.75 e-5
!
 do ifft=1,nfft
   x=(rtnewt-vtrial(ifft))/tphysel
   if (x.lt.2._dp) then
     y=exp(x)
     rhor_middx(ifft)=cktf*y*(23.1456e0_dp+y*(13.7820e0_dp+y))&
&     /(13.0586e0_dp+y*(17.0048e0_dp+y*(5.07527e0_dp+y*0.236620e0_dp)))
     rhor_mid(ifft)=cktf*y*(21.8168_dp+y*(13.1693_dp+y))&
&     /(24.6180+y*(23.5546_dp+y*(4.76290_dp+y*0.134481_dp)))
   else
     y=1._dp/(x*x)
     sqrtx=sqrt(x)
     rhor_middx(ifft)=cktf*sqrtx*(0.0153602e0_dp+y*(0.146815e0_dp+y))&
&     /(0.00768015e0_dp+y*(0.0763700e0_dp+y*0.570485e0_dp))
     rhor_mid(ifft)=cktf*x*sqrtx*(0.0473011_dp+y*(0.548433_dp+y))&
&     /(0.0709478_dp+y*(0.737041_dp+y*0.382065_dp))
   end if
 end do
!
!**********************************************************************
 end subroutine fm12a1t
!!***

end subroutine vtorhotf
!!***

end module m_vtorhotf
!!***
