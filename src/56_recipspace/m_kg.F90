!!****m* ABINIT/m_kg
!! NAME
!! m_kg
!!
!! FUNCTION
!!  Low-level functions to operate of G-vectors.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (DCA, XG, GMR, MT, DRH, AR)
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

MODULE m_kg

 use defs_basis
 use m_errors
 use m_abicore
 use m_errors
 use m_xmpi
 use m_dtset

 use defs_abitypes, only : MPI_type
 use m_fftcore,     only : kpgsph, bound
 use m_mpinfo,      only : proc_distrb_cycle

 implicit none

 private
!!***

 public :: getcut       ! Compute cutoff G^2
 public :: getmpw       ! Compute recommended npw from ecut, ucvol and gmet
 public :: mkkin        ! Compute elements of kinetic energy operator in reciprocal space at a given k point
 public :: kpgio        ! Do initialization of kg data.
 public :: ph1d3d       ! Compute the three-dimensional phase factor $e^{i 2 \pi (k+G) cdot xred}$
 public :: getph        ! Compute three factors of one-dimensional structure factor phase
 public :: kpgstr       ! Derivative of kinetic energy operator in reciprocal space.
 public :: mkkpg        ! Compute all (k+G) vectors (dp, in reduced coordinates) for given k point
 public :: mkpwind_k    ! Make plane wave index at k point for basis at second k point
 public :: mkkpgcart    ! Compute all (k+G) vectors (dp, in cartesian coordinates) for given k point
 public :: mkkin_metdqdq ! Compute the second q-gradient of the derivative of the kinetic energy operator w.r.t a metric

contains
!!***

!!****f* m_kg/getcut
!!
!! NAME
!! getcut
!!
!! FUNCTION
!! For input kpt, fft box dim ngfft(1:3), recip space metric gmet,
!! and kinetic energy cutoff ecut (hartree), COMPUTES:
!! if iboxcut==0:
!!   gsqcut: cut-off on G^2 for "large sphere" of radius double that
!!            of the basis sphere corresponding to ecut
!!   boxcut: where boxcut == gcut(box)/gcut(sphere).
!!                 boxcut >=2 for no aliasing.
!!                 boxcut < 1 is wrong and halts subroutine.
!! if iboxcut==1:
!!   gsqcut: cut-off on G^2 for "large sphere"
!!            containing the whole fft box
!!   boxcut: no meaning (zero)
!!
!! INPUTS
!! ecut=kinetic energy cutoff for planewave sphere (hartree)
!! gmet(3,3)=reciprocal space metric (bohr^-2)
!! iboxcut=0: compute gsqcut and boxcut with boxcut>=1
!!         1: compute gsqcut for boxcut=1 (sphere_cutoff=box_cutoff)
!! iout=unit number for output file
!! kpt(3)=input k vector (reduced coordinates--in terms of reciprocal lattice primitive translations)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!
!! OUTPUT
!! boxcut=defined above (dimensionless), ratio of basis sphere
!!  diameter to fft box length (smallest value)
!! gsqcut=Fourier cutoff on G^2 for "large sphere" of radius double
!!  that of the basis sphere--appropriate for charge density rho(G),
!!  Hartree potential, and pseudopotentials
!!
!! NOTES
!! 2*gcut arises from rho(g)=sum g prime (psi(g primt)*psi(g prime+g))
!!               where psi(g) is only nonzero for |g| <= gcut).
!! ecut (currently in hartree) is proportional to gcut(sphere)**2.
!!
!! PARENTS
!!      m_common,m_dfpt_looppert,m_dfpt_lw,m_dfpt_scfcv,m_nonlinear,m_pawfgr
!!      m_positron,m_respfn_driver,m_scfcv_core
!!
!! CHILDREN
!!
!! SOURCE

subroutine getcut(boxcut,ecut,gmet,gsqcut,iboxcut,iout,kpt,ngfft)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iboxcut,iout
 real(dp),intent(in) :: ecut
 real(dp),intent(out) :: boxcut,gsqcut
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3),kpt(3)

!Local variables-------------------------------
!scalars
 integer :: plane
 real(dp) :: boxsq,cutrad,ecut_pw,effcut,largesq,sphsq
 character(len=1000) :: message
!arrays
 integer :: gbound(3)

! *************************************************************************

!This is to treat the case where ecut has not been initialized,
!for wavelet computations. The default for ecut is -1.0 , allowed only for wavelets calculations
 ecut_pw=ecut
 if(ecut<-tol8)ecut_pw=ten

!gcut(box)**2=boxsq; gcut(sphere)**2=sphsq
!get min. d**2 to boundary of fft box:
!(gmet sets dimensions: bohr**-2)
!ecut(sphere)=0.5*(2 pi)**2 * sphsq:
 call bound(largesq,boxsq,gbound,gmet,kpt,ngfft,plane)
 effcut=0.5_dp * (two_pi)**2 * boxsq
 sphsq=2._dp*ecut_pw/two_pi**2

 if (iboxcut/=0) then
   boxcut=10._dp
   gsqcut=(largesq/sphsq)*(2.0_dp*ecut)/two_pi**2

   write(message, '(a,a,3f8.4,a,3i4,a,a,f11.3,a,a)' ) ch10,&
&   ' getcut: wavevector=',kpt,'  ngfft=',ngfft(1:3),ch10,&
&   '         ecut(hartree)=',ecut_pw+tol8,ch10,'=> whole FFT box selected'
   if(iout/=std_out) call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')

 else

!  Get G^2 cutoff for sphere of double radius of basis sphere
!  for selecting G s for rho(G), V_Hartree(G), and V_psp(G)--
!  cut off at fft box boundary or double basis sphere radius, whichever
!  is smaller.  If boxcut were 2, then relation would be
!$ecut_eff = (1/2) * (2 Pi Gsmall)^2 and gsqcut=4*Gsmall^2$.
   boxcut = sqrt(boxsq/sphsq)
   cutrad = min(2.0_dp,boxcut)
   gsqcut = (cutrad**2)*(2.0_dp*ecut_pw)/two_pi**2

   if(ecut>-tol8)then

     write(message, '(a,a,3f8.4,a,3i4,a,a,f11.3,3x,a,f10.5)' ) ch10,&
&     ' getcut: wavevector=',kpt,'  ngfft=',ngfft(1:3),ch10,&
&     '         ecut(hartree)=',ecut+tol8,'=> boxcut(ratio)=',boxcut+tol8
     if(iout/=std_out) call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')

     if (boxcut<1.0_dp) then
       write(message, '(9a,f12.6,6a)' )&
&       '  Choice of acell, ngfft, and ecut',ch10,&
&       '  ===> basis sphere extends BEYOND fft box !',ch10,&
&       '  Recall that boxcut=Gcut(box)/Gcut(sphere)  must be > 1.',ch10,&
&       '  Action: try larger ngfft or smaller ecut.',ch10,&
&       '  Note that ecut=effcut/boxcut**2 and effcut=',effcut+tol8,ch10,&
&       '  This situation might happen when optimizing the cell parameters.',ch10,&
&       '  Your starting geometry might be crazy.',ch10,&
&       '  See https://wiki.abinit.org/doku.php?id=howto:troubleshooting#incorrect_initial_geometry .'
       if(iout/=std_out) call wrtout(iout,message,'COLL')
       MSG_ERROR(message)
     end if

     if (boxcut>2.2_dp) then
       write(message, '(a,a,a,a,a,a,a,a,a,a,a,f12.6,a,a)' ) ch10,&
&       ' getcut : COMMENT -',ch10,&
&       '  Note that boxcut > 2.2 ; recall that',' boxcut=Gcut(box)/Gcut(sphere) = 2',ch10,&
&       '  is sufficient for exact treatment of convolution.',ch10,&
&       '  Such a large boxcut is a waste : you could raise ecut',ch10,&
&       '  e.g. ecut=',effcut*0.25_dp+tol8,' Hartrees makes boxcut=2',ch10
       if(iout/=std_out) call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if

     if (boxcut<1.5_dp) then
       write(message, '(15a)' ) ch10,&
&       ' getcut : WARNING -',ch10,&
&       '  Note that boxcut < 1.5; this usually means',ch10,&
&       '  that the forces are being fairly strongly affected by','  the smallness of the fft box.',ch10,&
&       '  Be sure to test with larger ngfft(1:3) values.',ch10,&
&       '  This situation might happen when optimizing the cell parameters.',ch10,&
&       '  Your starting geometry might be crazy.',ch10,&
&       '  See https://wiki.abinit.org/doku.php?id=howto:troubleshooting#incorrect_initial_geometry .'
       if(iout/=std_out) call wrtout(iout,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if

   end if

 end if  ! iboxcut

end subroutine getcut
!!***

!!****f* m_kg/getmpw
!! NAME
!! getmpw
!!
!! FUNCTION
!! From input ecut, combined with ucvol and gmet, compute recommended mpw
!! mpw is the maximum number of plane-waves in the wave-function basis
!!  for one processor of the WF group
!!
!! INPUTS
!! ecut=plane-wave cutoff energy in Hartrees
!! exchn2n3d=if n1, n2 and n3 are exchanged
!! gmet(3,3)=reciprocal space metric (bohr**-2).
!! istwfk(nkpt)=input option parameter that describes the storage of wfs
!! kptns(3,nkpt)=real(dp) array for k points (normalisation is already taken into account)
!! mpi_enreg=information about MPI parallelization
!! nkpt=integer number of k points in the calculation
!!
!! OUTPUT
!! mpw=maximal number of plane waves over all k points of the processor
!!  (for one processor of the WF group)
!!
!! PARENTS
!!      m_dfpt_looppert,m_dfpt_lw,m_memeval,m_mpi_setup,m_scfcv_core
!!
!! CHILDREN
!!
!! SOURCE

subroutine getmpw(ecut,exchn2n3d,gmet,istwfk,kptns,mpi_enreg,mpw,nkpt)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: exchn2n3d,nkpt
 integer,intent(out) :: mpw
 real(dp),intent(in) :: ecut
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: istwfk(nkpt)
 real(dp),intent(in) :: gmet(3,3),kptns(3,nkpt)

!Local variables-------------------------------
!scalars
 integer :: ikpt,istwf_k,npw
! integer :: npwwrk,pad=50
! real(dp) :: scale=1.3_dp
 character(len=500) :: message
!arrays
 integer,allocatable :: kg(:,:)
 real(dp) :: kpoint(3)

! *************************************************************************

!An upper bound for mpw, might be obtained as follows
!the average number of plane-waves in the cutoff sphere is
!npwave = (2*ecut)**(3/2) * ucvol / (6*pi**2)
!the upper bound is calculated as
!npwwrk = int(scale * npwave) + pad
!rescale so an upper bound
!npwave=nint(ucvol*(2.0_dp*ecut)**1.5_dp/(6.0_dp*pi**2))
!npwwrk=nint(dble(npwave)*scale)+pad

 ABI_ALLOCATE(kg,(3,100))

!set mpw to zero, as needed for only counting in kpgsph
 mpw = 0

!Might be parallelized over k points ? !
 do ikpt = 1,nkpt
!  Do computation of G sphere, returning npw
   kpoint(:)=kptns(:,ikpt)
   istwf_k=istwfk(ikpt)
   call kpgsph(ecut,exchn2n3d,gmet,0,ikpt,istwf_k,kg,kpoint,0,mpi_enreg,0,npw)
   mpw = max(npw,mpw)
 end do

 write(message,'(a,i0)') ' getmpw: optimal value of mpw= ',mpw
 call wrtout(std_out,message,'COLL')

 ABI_DEALLOCATE(kg)

end subroutine getmpw
!!***

!!****f* m_kg/mkkin
!! NAME
!! mkkin
!!
!! FUNCTION
!! compute elements of kinetic energy operator in reciprocal space at a given k point
!!
!! INPUTS
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!  effmass_free=effective mass for electrons (1. in common case)
!!  gmet(3,3)=reciprocal lattice metric tensor ($\textrm{Bohr}^{-2}$)
!!  idir1 = 1st direction of the derivative (if 1 <= idir1 <= 3, not used otherwise)
!!  idir2 = 2st direction of the derivative (if 1 <= idir1,idir2 <= 3, not used otherwise))
!!  kg(3,npw)=integer coordinates of planewaves in basis sphere.
!!  kpt(3)=reduced coordinates of k point
!!  npw=number of plane waves at kpt.
!!
!! OUTPUT
!!  kinpw(npw)=(modified) kinetic energy (or derivative) for each plane wave (Hartree)
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
!! This smearing factor is also used to derived a modified kinetic
!! contribution to stress, in another routine (forstrnps.f)
!!
!! Also, in order to break slightly the symmetry between axes, that causes
!! sometimes a degeneracy of eigenvalues and do not allow to obtain
!! the same results on different machines, there is a modification
!! by one part over 1.0e12 of the metric tensor elements (1,1) and (3,3)
!!
!! PARENTS
!!      m_d2frnl,m_dfpt_nstwf,m_dfpt_scfcv,m_dfptnl_pert,m_dft_energy,m_getgh1c
!!      m_io_kss,m_ksdiago,m_mkffnl,m_orbmag,m_vhxc_me,m_vkbr,m_vtorho
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkkin (ecut,ecutsm,effmass_free,gmet,kg,kinpw,kpt,npw,idir1,idir2)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw
 integer,intent(in) :: idir1,idir2
 real(dp),intent(in) :: ecut,ecutsm,effmass_free
!arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: gmet(3,3),kpt(3)
 real(dp),intent(out) :: kinpw(npw)

!Local variables-------------------------------
!scalars
 integer :: ig,order
 real(dp),parameter :: break_symm=1.0d-11
 real(dp) :: ecutsm_inv,fsm,gpk1,gpk2,gpk3,htpisq,kinetic,kpg2,dkpg2,xx
 real(dp) :: d1kpg2,d2kpg2,ddfsm, dfsm
!arrays
 real(dp) :: gmet_break(3,3)

! *************************************************************************
!
!htpisq is (1/2) (2 Pi) **2:
 htpisq=0.5_dp*(two_pi)**2

 ecutsm_inv=0.0_dp
 if(ecutsm>1.0d-20)ecutsm_inv=1/ecutsm

 gmet_break(:,:)=gmet(:,:)
 gmet_break(1,1)=(1.0_dp+break_symm)*gmet(1,1)
 gmet_break(3,3)=(1.0_dp-break_symm)*gmet(3,3)

 order=0 ! Compute the kinetic operator
 if (idir1>0.and.idir1<4) then
   order=1 ! Compute the 1st derivative of the kinetic operator
   if (idir2>0.and.idir2<4) then
     order=2 ! Compute the 2nd derivative of the kinetic operator
   end if
 end if

!$OMP PARALLEL DO PRIVATE(dkpg2,d1kpg2,d2kpg2,gpk1,gpk2,gpk3,ig,kinetic,kpg2,xx,fsm,dfsm,ddfsm) &
!$OMP SHARED(kinpw,ecut,ecutsm,ecutsm_inv) &
!$OMP SHARED(gmet_break,htpisq,idir1,idir2,kg,kpt,npw)
 do ig=1,npw
   gpk1=dble(kg(1,ig))+kpt(1)
   gpk2=dble(kg(2,ig))+kpt(2)
   gpk3=dble(kg(3,ig))+kpt(3)
   kpg2=htpisq*&
&   ( gmet_break(1,1)*gpk1**2+         &
&   gmet_break(2,2)*gpk2**2+         &
&   gmet_break(3,3)*gpk3**2          &
&   +2.0_dp*(gpk1*gmet_break(1,2)*gpk2+  &
&   gpk1*gmet_break(1,3)*gpk3+  &
&   gpk2*gmet_break(2,3)*gpk3 )  )
   select case (order)
   case(0)
     kinetic=kpg2
   case(1)
     dkpg2=htpisq*2.0_dp*&
&     (gmet_break(idir1,1)*gpk1+gmet_break(idir1,2)*gpk2+gmet_break(idir1,3)*gpk3)
     kinetic=dkpg2
   case(2)
     dkpg2=htpisq*2.0_dp*gmet_break(idir1,idir2)
     kinetic=dkpg2
   end select
   if(kpg2>ecut-ecutsm)then
     if(kpg2>ecut-tol12)then
       if(order==0) then
!        Will filter the wavefunction, based on this value, in cgwf.f, getghc.f and precon.f
         kinetic=huge(0.0_dp)*1.d-10
       else
!        The wavefunction has been filtered : no derivative
         kinetic=0
       end if
     else
       if(order==0) then
         xx=max( (ecut-kpg2)*ecutsm_inv , 1.0d-20)
       else
         xx=(ecut-kpg2)*ecutsm_inv
       end if
       if(order==2) then
         d1kpg2=htpisq*2.0_dp*&
&         (gmet_break(idir1,1)*gpk1+gmet_break(idir1,2)*gpk2+gmet_break(idir1,3)*gpk3)
         d2kpg2=htpisq*2.0_dp*&
&         (gmet_break(idir2,1)*gpk1+gmet_break(idir2,2)*gpk2+gmet_break(idir2,3)*gpk3)
       end if
!      This kinetic cutoff smoothing function and its xx derivatives
!      were produced with Mathematica and the fortran code has been
!      numerically checked against Mathematica.
       fsm=1.0_dp/(xx**2*(3+xx*(1+xx*(-6+3*xx))))
       if(order>0) dfsm=-3.0_dp*(-1+xx)**2*xx*(2+5*xx)*fsm**2
       if(order>1) ddfsm=6.0_dp*xx**2*(9+xx*(8+xx*(-52+xx*(-3+xx*(137+xx*(-144+45*xx))))))*fsm**3
       select case (order)
       case(0)
         kinetic=kpg2*fsm
       case(1)
         kinetic=dkpg2*(fsm-ecutsm_inv*kpg2*dfsm)
       case(2)
         kinetic=dkpg2*fsm&
&         -2.0_dp*d1kpg2*dfsm*ecutsm_inv*d2kpg2&
&         +kpg2*ddfsm*(ecutsm_inv**2)*d1kpg2*d2kpg2&
&         -kpg2*dfsm*ecutsm_inv*dkpg2
       end select
     end if
   end if
   kinpw(ig)=kinetic/effmass_free
 end do
!$OMP END PARALLEL DO

end subroutine mkkin
!!***

!!****f* m_kg/kpgio
!! NAME
!! kpgio
!!
!! FUNCTION
!! Do initialization of kg data.
!!
!! INPUTS
!!  ecut=kinetic energy planewave cutoff (hartree)
!!  exchn2n3d=if 1, n2 and n3 are exchanged
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kptns(3,nkpt)=reduced coords of k points
!!  mkmem =number of k points treated by this node.
!!  character(len=4) : mode_paral=either 'COLL' or 'PERS', tells whether
!!   the loop over k points must be done by all processors or not,
!!   in case of parallel execution.
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum number of planewaves as dimensioned in calling routine
!!  nband(nkpt*nsppol)=number of bands at each k point
!!  nkpt=number of k points
!!  nsppol=1 for unpolarized, 2 for polarized
!!
!! OUTPUT
!!  npwarr(nkpt)=array holding npw for each k point, taking into account
!!   the effect of istwfk, and the spreading over processors
!!  npwtot(nkpt)=array holding the total number of plane waves for each k point,
!!  kg(3,mpw*mkmem)=dimensionless coords of G vecs in basis sphere at k point
!!
!! NOTES
!! Note that in case of band parallelism, the number of spin-up
!! and spin-down bands must be equal at each k points
!!
!! PARENTS
!!      m_cut3d,m_dfpt_fef,m_dfpt_looppert,m_dfpt_lw,m_gstate,m_inwffil
!!      m_longwave,m_nonlinear,m_respfn_driver,m_scfcv_core
!!
!! CHILDREN
!!
!! SOURCE

subroutine kpgio(ecut,exchn2n3d,gmet,istwfk,kg,kptns,mkmem,nband,nkpt,&
& mode_paral,mpi_enreg,mpw,npwarr,npwtot,nsppol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: exchn2n3d,mkmem,mpw,nkpt,nsppol
 real(dp),intent(in) :: ecut
 character(len=4),intent(in) :: mode_paral
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: istwfk(nkpt),nband(nkpt*nsppol)
 integer,intent(out) :: kg(3,mpw*mkmem),npwarr(nkpt),npwtot(nkpt)
 real(dp),intent(in) :: gmet(3,3),kptns(3,nkpt)

!Local variables-------------------------------
!scalars
 integer :: ierr,ikg,ikpt,istwf_k,me,nband_down,nband_k,npw1
 logical :: test_npw
 character(len=500) :: message
!arrays
 real(dp) :: kpoint(3)

! *************************************************************************

!Define me
 me=mpi_enreg%me_kpt

 if((mpi_enreg%paralbd==1) .and. (mode_paral=='PERS')) then
   if(nsppol==2)then
     do ikpt=1,nkpt
       nband_k=nband(ikpt)
       nband_down=nband(ikpt+nkpt)
       if(nband_k/=nband_down)then
         write(message,'(a,a,a,a,a,a,a,a,i4,a,i4,a,a,a,i4,a,a,a)')ch10,&
&         ' kpgio : ERROR -',ch10,&
&         '  Band parallel case, one must have same number',ch10,&
&         '  of spin up and spin down bands, but input is :',ch10,&
&         '  nband(up)=',nband_k,', nband(down)=',nband_down,',',ch10,&
&         '  for ikpt=',ikpt,'.',ch10,&
&         '  Action: correct nband in your input file.'
!        MG: Tests v3(10,11,17) and v6(67) fail if this test is enabled
!        call wrtout(std_out,message,mode_paral)
       end if
     end do
   end if
 end if

 npwarr(:)=0
 npwtot(:)=0

 kg=0
 ikg=0
!Find (k+G) sphere for each k.

 do ikpt=1,nkpt

   nband_k = nband(ikpt)

   if(mode_paral=='PERS')then
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,-1,me)) cycle
   end if

   kpoint(:)=kptns(:,ikpt)
   istwf_k=istwfk(ikpt)
   call kpgsph(ecut,exchn2n3d,gmet,ikg,ikpt,istwf_k,kg,kpoint,mkmem,mpi_enreg,mpw,npw1)

   test_npw=.true.
   if (xmpi_paral==1)then
     if (mode_paral=='PERS')then
       test_npw=(minval(mpi_enreg%proc_distrb(ikpt,1:nband_k,1:nsppol))==me)
     end if
   end if
   if (test_npw) npwarr(ikpt)=npw1

!  Make sure npw < nband never happens:
!  if (npw1<nband(ikpt)) then
!  write(message, '(a,a,a,a,i5,a,3f8.4,a,a,i10,a,i10,a,a,a,a)' )ch10,&
!  &   ' kpgio : ERROR -',ch10,&
!  &   '  At k point number',ikpt,' k=',(kptns(mu,ikpt),mu=1,3),ch10,&
!  &   '  npw=',npw1,' < nband=',nband(ikpt),ch10,&
!  &   '  Indicates not enough planewaves for desired number of bands.',ch10,&
!  &   '  Action: change either ecut or nband in input file.'
!  MSG_ERROR(message)
!  end if

!  Find boundary of G sphere for efficient zero padding,
!    Shift to next section of each array kg
   ikg=ikg+npw1
 end do !  End of the loop over k points

 if(mode_paral == 'PERS') then
   call xmpi_sum(npwarr,mpi_enreg%comm_kpt,ierr)
 end if

 !if (mpi_enreg%nproc>1) call wrtout(std_out,' kpgio: loop on k-points done in parallel','COLL')

!XG030513 MPIWF : now, one should sum npwarr over all processors
!of the WF group, to get npwtot (to be spread on all procs of the WF group
 npwtot(:)=npwarr(:)

!Taking into account istwfk
 do ikpt=1,nkpt
   if(istwfk(ikpt)>1)then
     if(istwfk(ikpt)==2)then
       npwtot(ikpt)=2*npwtot(ikpt)-1
     else
       npwtot(ikpt)=2*npwtot(ikpt)
     end if
   end if
 end do

end subroutine kpgio
!!***

!!****f* m_kg/ph1d3d
!! NAME
!! ph1d3d
!!
!! FUNCTION
!! Compute the three-dimensional phase factor $e^{i 2 \pi (k+G) cdot xred}$
!! from the three one-dimensional factors, the k point coordinates,
!! and the atom positions, for all planewaves which fit in the fft box.
!!
!! INPUTS
!!  iatom, jatom= bounds of atom indices in ph1d for which ph3d has to be computed
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  matblk= dimension of ph3d
!!  natom= dimension of ph1d
!!  npw=number of plane waves
!!  n1,n2,n3=dimensions of fft box (ngfft(3)).
!!  phkxred(2,natom)=phase factors exp(2 pi k.xred)
!!  ph1d(2,(2*n1+1)*natom+(2*n2+1)*natom+(2*n3+1)*natom)=exp(2Pi i G xred) for
!!   vectors (Gx,0,0), (0,Gy,0) and (0,0,Gz)
!!   with components ranging from -nj <= Gj <= nj
!!
!! OUTPUT
!!  ph3d(2,npw_k,matblk)=$e^{2 i \pi (k+G) cdot xred}$ for vectors (Gx,Gy,Gz),
!!   and for atoms in the range iatom to jatom with respect to ph1d
!!
!! PARENTS
!!      m_cgprj,m_cgtk,m_cut3d,m_epjdos,m_fock,m_hamiltonian,m_inwffil
!!      m_nonlop_pl,m_nonlop_ylm,m_orbmag,m_suscep_stat,m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine ph1d3d(iatom,jatom,kg_k,matblk,natom,npw_k,n1,n2,n3,phkxred,ph1d,ph3d)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iatom,jatom,matblk,n1,n2,n3,natom,npw_k
!arrays
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(in) :: ph1d(2,(2*n1+1+2*n2+1+2*n3+1)*natom)
 real(dp),intent(in) :: phkxred(2,natom)
 real(dp),intent(out) :: ph3d(2,npw_k,matblk)

!Local variables-------------------------------
!scalars
 integer :: i1,ia,iatblk,ig,shift1,shift2,shift3
 real(dp) :: ph12i,ph12r,ph1i,ph1r,ph2i,ph2r,ph3i,ph3r,phkxi,phkxr
 character(len=500) :: message
!arrays
 real(dp),allocatable :: ph1kxred(:,:)

! *************************************************************************

 if(matblk-1 < jatom-iatom)then
   write(message,'(a,a,a,a,a,i0,a,a,i0,a,i0,a)')&
&   'Input natom-1 must be larger or equal to jatom-iatom,',ch10,&
&   'while their value is : ',ch10,&
&   'natom-1 = ',natom-1,ch10,&
&   'jatom=',jatom,', iatom=',iatom,'.'
   MSG_BUG(message)
 end if

 ABI_ALLOCATE(ph1kxred,(2,-n1:n1))

 ! ia runs from iatom to jatom
 do ia=iatom,jatom

   ! iatblk runs from 1 to matblk
   iatblk=ia-iatom+1
   shift1=1+n1+(ia-1)*(2*n1+1)
   shift2=1+n2+(ia-1)*(2*n2+1)+natom*(2*n1+1)
   shift3=1+n3+(ia-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
   ! Compute product of phkxred by phase for the first component of G vector
   phkxr=phkxred(1,ia)
   phkxi=phkxred(2,ia)
   do i1=-n1,n1
     ph1kxred(1,i1)=ph1d(1,i1+shift1)*phkxr-ph1d(2,i1+shift1)*phkxi
     ph1kxred(2,i1)=ph1d(2,i1+shift1)*phkxr+ph1d(1,i1+shift1)*phkxi
   end do

   ! Compute tri-dimensional phase factor
!$OMP PARALLEL DO PRIVATE(ig,ph1r,ph1i,ph2r,ph2i,ph3r,ph3i,ph12r,ph12i)
   do ig=1,npw_k
     ph1r=ph1kxred(1,kg_k(1,ig))
     ph1i=ph1kxred(2,kg_k(1,ig))
     ph2r=ph1d(1,kg_k(2,ig)+shift2)
     ph2i=ph1d(2,kg_k(2,ig)+shift2)
     ph3r=ph1d(1,kg_k(3,ig)+shift3)
     ph3i=ph1d(2,kg_k(3,ig)+shift3)
     ph12r=ph1r*ph2r-ph1i*ph2i
     ph12i=ph1r*ph2i+ph1i*ph2r
     ph3d(1,ig,iatblk)=ph12r*ph3r-ph12i*ph3i
     ph3d(2,ig,iatblk)=ph12r*ph3i+ph12i*ph3r
   end do
!$OMP END PARALLEL DO
 end do

 ABI_DEALLOCATE(ph1kxred)

end subroutine ph1d3d
!!***

!!****f* m_kg/getph
!!
!! NAME
!! getph
!!
!! FUNCTION
!! Compute three factors of one-dimensional structure factor phase
!! for input atomic coordinates, for all planewaves which fit in fft box.
!! The storage of these atomic factors is made according to the
!! values provided by the index table atindx. This will save time in nonlop.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  natom=number of atoms in cell.
!!  n1,n2,n3=dimensions of fft box (ngfft(3)).
!!  xred(3,natom)=reduced atomic coordinates.
!!
!! OUTPUT
!!  ph1d(2,(2*n1+1)*natom+(2*n2+1)*natom+(2*n3+1)*natom)=exp(2Pi i G.xred) for
!!   integer vector G with components ranging from -nj <= G <= nj.
!!   Real and imag given in usual Fortran convention.
!!
!! PARENTS
!!      m_afterscfloop,m_berryphase_new,m_bethe_salpeter,m_cgtk,m_cut3d
!!      m_dfpt_looppert,m_dfpt_lw,m_epjdos,m_extraprho,m_fock,m_gkk,m_gstate
!!      m_hamiltonian,m_inwffil,m_nonlinear,m_orbmag,m_pead_nl_loop,m_phgamma
!!      m_phpi,m_prcref,m_respfn_driver,m_scfcv_core,m_screening_driver
!!      m_sigma_driver,m_sigmaph,m_wfd
!!
!! CHILDREN
!!
!! SOURCE

subroutine getph(atindx,natom,n1,n2,n3,ph1d,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,natom
!arrays
 integer,intent(in) :: atindx(natom)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(out) :: ph1d(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,ia,ii,ph1d_size1,ph1d_size2,ph1d_sizemin
 !character(len=500) :: msg
 real(dp) :: arg

! *************************************************************************

 ph1d_size1=size(ph1d,1);ph1d_size2=size(ph1d,2)
 ph1d_sizemin=(2*n1+1+2*n2+1+2*n3+1)*natom
 if (ph1d_size1/=2.or.ph1d_size2<ph1d_sizemin) then
   MSG_BUG('Wrong ph1d sizes!')
 end if

 do ia=1,natom

!  Store the phase factor of atom number ia in place atindx(ia)
   i1=(atindx(ia)-1)*(2*n1+1)
   i2=(atindx(ia)-1)*(2*n2+1)+natom*(2*n1+1)
   i3=(atindx(ia)-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)

   do ii=1,2*n1+1
     arg=two_pi*dble(ii-1-n1)*xred(1,ia)
     ph1d(re,ii+i1)=dcos(arg)
     ph1d(im,ii+i1)=dsin(arg)
   end do

   do ii=1,2*n2+1
     arg=two_pi*dble(ii-1-n2)*xred(2,ia)
     ph1d(re,ii+i2)=dcos(arg)
     ph1d(im,ii+i2)=dsin(arg)
   end do

   do ii=1,2*n3+1
     arg=two_pi*dble(ii-1-n3)*xred(3,ia)
     ph1d(re,ii+i3)=dcos(arg)
     ph1d(im,ii+i3)=dsin(arg)
   end do

 end do

!This is to avoid uninitialized ph1d values
 if (ph1d_sizemin<ph1d_size2) then
   ph1d(:,ph1d_sizemin+1:ph1d_size2)=zero
 end if

end subroutine getph
!!***

!!****f* m_kg/kpgstr
!! NAME
!! kpgstr
!!
!! FUNCTION
!! Compute elements of the derivative the kinetic energy operator in reciprocal
!! space at given k point wrt a single cartesian strain component
!!
!! INPUTS
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!  effmass_free=effective mass for electrons (1. in common case)
!!  gmet(3,3) = reciprocal lattice metric tensor (Bohr**-2)
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  istr=1,...6 specifies cartesian strain component 11,22,33,32,31,21
!!  kg(3,npw) = integer coordinates of planewaves in basis sphere.
!!  kpt(3)    = reduced coordinates of k point
!!  npw       = number of plane waves at kpt.
!!
!! OUTPUT
!!  dkinpw(npw)=d/deps(istr) ( (1/2)*(2 pi)**2 * (k+G)**2 )
!!
!! NOTES
!!  Src_6response/kpg3.f
!!
!! PARENTS
!!      m_dfpt_nstwf,m_dfpt_scfcv,m_getgh1c
!!
!! CHILDREN
!!
!! SOURCE

subroutine kpgstr(dkinpw,ecut,ecutsm,effmass_free,gmet,gprimd,istr,kg,kpt,npw)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: istr,npw
 real(dp),intent(in) :: ecut,ecutsm,effmass_free
!arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),kpt(3)
 real(dp),intent(out) :: dkinpw(npw)

!Local variables -------------------------
!scalars
 integer :: ig,ii,ka,kb
 real(dp) :: dfsm,dkinetic,dkpg2,ecutsm_inv,fsm,gpk1,gpk2,gpk3,htpisq
! real(dp) :: d2fsm ! used in commented section below
 real(dp) :: kpg2,xx
 character(len=500) :: message
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 real(dp) :: dgmetds(3,3)

! *********************************************************************

!htpisq is (1/2) (2 Pi) **2:
 htpisq=0.5_dp*(two_pi)**2

 ecutsm_inv=0.0_dp
 if(ecutsm>1.0d-20)ecutsm_inv=1/ecutsm

!Compute derivative of metric tensor wrt strain component istr
 if(istr<1 .or. istr>6)then
   write(message, '(a,i10,a,a,a)' )&
&   'Input istr=',istr,' not allowed.',ch10,&
&   'Possible values are 1,2,3,4,5,6 only.'
   MSG_BUG(message)
 end if

 ka=idx(2*istr-1);kb=idx(2*istr)
 do ii = 1,3
   dgmetds(:,ii)=-(gprimd(ka,:)*gprimd(kb,ii)+gprimd(kb,:)*gprimd(ka,ii))
 end do
!For historical reasons:
 dgmetds(:,:)=0.5_dp*dgmetds(:,:)

 do ig=1,npw
   gpk1=dble(kg(1,ig))+kpt(1)
   gpk2=dble(kg(2,ig))+kpt(2)
   gpk3=dble(kg(3,ig))+kpt(3)
   kpg2=htpisq*&
&   ( gmet(1,1)*gpk1**2+         &
&   gmet(2,2)*gpk2**2+         &
&   gmet(3,3)*gpk3**2          &
&   +2.0_dp*(gpk1*gmet(1,2)*gpk2+  &
&   gpk1*gmet(1,3)*gpk3+  &
&   gpk2*gmet(2,3)*gpk3 )  )
   dkpg2=htpisq*2.0_dp*&
&   (gpk1*(dgmetds(1,1)*gpk1+dgmetds(1,2)*gpk2+dgmetds(1,3)*gpk3)+  &
&   gpk2*(dgmetds(2,1)*gpk1+dgmetds(2,2)*gpk2+dgmetds(2,3)*gpk3)+  &
&   gpk3*(dgmetds(3,1)*gpk1+dgmetds(3,2)*gpk2+dgmetds(3,3)*gpk3) )
   dkinetic=dkpg2
   if(kpg2>ecut-ecutsm)then
     if(kpg2>ecut-tol12)then
!      The wavefunction has been filtered : no derivative
       dkinetic=0.0_dp
     else
       xx=(ecut-kpg2)*ecutsm_inv
!      This kinetic cutoff smoothing function and its xx derivatives
!      were produced with Mathematica and the fortran code has been
!      numerically checked against Mathematica.
       fsm=1.0_dp/(xx**2*(3+xx*(1+xx*(-6+3*xx))))
       dfsm=-3.0_dp*(-1+xx)**2*xx*(2+5*xx)*fsm**2
!      d2fsm=6.0_dp*xx**2*(9+xx*(8+xx*(-52+xx*(-3+xx*(137+xx*&
!      &                        (-144+45*xx))))))*fsm**3
       dkinetic=dkpg2*(fsm-ecutsm_inv*kpg2*dfsm)
     end if
   end if
   dkinpw(ig)=dkinetic/effmass_free
 end do

end subroutine kpgstr
!!***

!!****f* m_kg/mkkpg
!! NAME
!! mkkpg
!!
!! FUNCTION
!! Compute all (k+G) vectors (dp, in reduced coordinates) for given k point,
!! from integer coordinates of G vectors.
!! Eventually compute related data.
!!
!! INPUTS
!!  kg(3,npw)=integer coords of planewaves in basis sphere
!!  kpt(3)=k point in terms of recip. translations
!!  nkpg=second dimension of array kpg
!!  npw=number of plane waves in reciprocal space
!!
!! OUTPUT
!!  kpg(npw,3)= (k+G) components
!!  === if nkpg==9 ===
!!    kpg(npw,4:9)= [(k+G)_a].[(k+G)_b] quantities
!!
!! PARENTS
!!      m_bandfft_kpt,m_cgprj,m_d2frnl,m_dfpt_lwwf,m_dfpt_nstwf,m_dfpt_scfcv
!!      m_dfptnl_pert,m_fock_getghc,m_forstr,m_getgh1c,m_io_kss,m_ksdiago
!!      m_nonlop_test,m_nonlop_ylm,m_orbmag,m_pead_nl_loop,m_vtorho,m_wfd
!!
!! CHILDREN
!!
!! SOURCE


subroutine mkkpg(kg,kpg,kpt,nkpg,npw)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpg,npw
!arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(out) :: kpg(npw,nkpg)

!Local variables-------------------------------
!scalars
 integer :: ipw,mu,mua,mub
 character(len=500) :: message
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)

! *************************************************************************

 DBG_ENTER("COLL")

 if (nkpg==0) return

!-- Test nkpg --
 if (nkpg/=3.and.nkpg/=9) then
   write(message, '(a,i0)' )' Bad value for nkpg !',nkpg
   MSG_BUG(message)
 end if

!-- Compute (k+G) --
!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP PRIVATE(mu,ipw)
 do ipw=1,npw
   do mu=1,3
     kpg(ipw,mu)=kpt(mu)+dble(kg(mu,ipw))
   end do
 end do
!$OMP END PARALLEL DO

!-- Compute [(k+G)_a].[(k+G)_b] --
 if (nkpg==9) then
!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP PRIVATE(ipw,mu,mua,mub)
   do ipw=1,npw
     do mu=4,9
       mua=alpha(mu-3);mub=beta(mu-3)
       kpg(ipw,mu)=kpg(ipw,mua)*kpg(ipw,mub)
     end do
   end do
!$OMP END PARALLEL DO
 end if

 DBG_EXIT("COLL")

end subroutine mkkpg
!!***

!!****f* ABINIT/mkpwind_k
!! NAME
!! mkpwind_k
!!
!! FUNCTION
!! Make plane wave index at k point for basis at second k point,
!! needed to compute overlaps $\langle u_{k,n}|u_{k+b,n}\rangle$
!! as appear in Berry phase derived quantities
!!
!! COPYRIGHT
!! Copyright (C) 2003-2020 ABINIT  group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! dk(3)=real vector difference of ket kpt - bra kpt
!! dtset <type(dataset_type)>=all input variables in this dataset
!! fnkpt=number of kpts in full BZ
!! fkptns=kpts in full BZ
!! gmet(3,3)=metric in reciprocal space
!! indkk_f2ibz(fnkpt,6)=information on folding from FBZ to IBZ (see initberry or initorbmag)
!! ikpt=index of bra k pt in FBZ
!! ikpt1=index of neighbour ket k pt in FBZ
!! mpi_enreg=information about MPI parallelization
!! npwarr(dtset%nkpt)=npw at each kpt
!! symrec(3,3,nsym) = symmetries in reciprocal space in terms of
!!   reciprocal space primitive translations
!!
!! OUTPUT
!! pwind_k1(dtset%mpw)=output index of ikpt1 basis states refered to ikpt
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! NOTES
!!
!! PARENTS
!!      m_orbmag
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkpwind_k(dk,dtset,fnkpt,fkptns,gmet,indkk_f2ibz,ikpt,ikpt1,&
& mpi_enreg,npwarr,pwind_k1,symrec)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: fnkpt,ikpt,ikpt1
  type(dataset_type),intent(in) :: dtset
  type(MPI_type), intent(inout) :: mpi_enreg

  !arrays
  integer,intent(in) :: indkk_f2ibz(fnkpt,6)
  integer,intent(in) :: npwarr(dtset%nkpt)
  integer,intent(in) :: symrec(3,3,dtset%nsym)
  integer,intent(out) :: pwind_k1(dtset%mpw)
  real(dp),intent(in) :: dk(3),fkptns(3,fnkpt),gmet(3,3)

  !Local variables -------------------------
  !scalars
  integer :: exchn2n3d,idum1,ikg1,ikpti,ikpt1i,ipw,istwf_k,isym,isym1,jpw,npw_k,npw_k1
  real(dp) :: ecut_eff

  !arrays
  integer,allocatable :: kg_k(:,:),kg1_k(:,:)
  real(dp) :: dg(3),dum33(3,3),kpt(3),kpt1(3),iadum(3),iadum1(3)

  ! ***********************************************************************

  ikpti = indkk_f2ibz(ikpt,1)
  ikpt1i = indkk_f2ibz(ikpt1,1)


  ecut_eff = dtset%ecut*(dtset%dilatmx)**2
  exchn2n3d = 0 ; istwf_k = 1 ; ikg1 = 0

  ! Build basis sphere of plane waves for the k-point
  ! we avoid using the global kg data because of difficulties in parallel-ism
  ABI_ALLOCATE(kg_k,(3,dtset%mpw))
  kg_k(:,:) = 0
  kpt(:) = dtset%kptns(:,ikpti)
  call kpgsph(ecut_eff,exchn2n3d,gmet,ikg1,ikpt,istwf_k,kg_k,kpt,1,mpi_enreg,dtset%mpw,npw_k)

  ! Build basis sphere of plane waves for the nearest neighbour of the k-point
  ABI_ALLOCATE(kg1_k,(3,dtset%mpw))
  kg1_k(:,:) = 0
  kpt1(:) = dtset%kptns(:,ikpt1i)
  call kpgsph(ecut_eff,exchn2n3d,gmet,ikg1,ikpt,istwf_k,kg1_k,kpt1,1,mpi_enreg,dtset%mpw,npw_k1)

  !        Deal with symmetry transformations
  !

  !        bra k-point k(b) and IBZ k-point kIBZ(b) related by
  !        k(b) = alpha(b) S(b)^t kIBZ(b) + G(b)
  !        where alpha(b), S(b) and G(b) are given by indkk_f2ibz
  !
  !        For the ket k-point:
  !        k(k) = alpha(k) S(k)^t kIBZ(k) + G(k) - GBZ(k)
  !        where GBZ(k) takes k(k) to the BZ
  !

  isym  = indkk_f2ibz(ikpt,2)
  isym1 = indkk_f2ibz(ikpt1,2)

  !        Construct transformed G vector that enters the matching condition:
  !        alpha(k) S(k)^{t,-1} ( -G(b) - GBZ(k) + G(k) )

  dg(:) = -indkk_f2ibz(ikpt,3:5) &
       & - nint(-fkptns(:,ikpt) - dk(:) - tol10 + fkptns(:,ikpt1)) &
       & + indkk_f2ibz(ikpt1,3:5)

  iadum(:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),dg(:))

  dg(:) = iadum(:)

  !        Construct S(k)^{t,-1} S(b)^{t}

  dum33(:,:) = MATMUL(TRANSPOSE(dtset%symrel(:,:,isym1)),symrec(:,:,isym))

  !        Construct alpha(k) alpha(b)

  pwind_k1(:) = 0
  npw_k = npwarr(ikpti)
  do ipw = 1, npw_k

     !          NOTE: the bra G vector is taken for the sym-related IBZ k point,
     !          not for the FBZ k point

     ! original code from initberry
     ! iadum(:) = kg(:,kgindex(ikpti) + ipw)

     iadum(:) = kg_k(:,ipw)

     !          to determine r.l.v. matchings, we transformed the bra vector
     !          Rotation
     iadum1(:)=0
     do idum1=1,3
        iadum1(:)=iadum1(:)+dum33(:,idum1)*iadum(idum1)
     end do
     iadum(:)=iadum1(:)
     iadum(:) = iadum(:) + dg(:)

     do jpw = 1, npw_k1
        iadum1(1:3) = kg1_k(1:3,jpw)
        if ( (iadum(1) == iadum1(1)).and. &
             &     (iadum(2) == iadum1(2)).and. &
             &     (iadum(3) == iadum1(3)) ) then
           pwind_k1(ipw) = jpw
           ! write(std_out,'(a,2i4)')'JWZ debug : bg ipw == jpw ',ipw,jpw
           exit
        end if
     end do
  end do

  ABI_DEALLOCATE(kg_k)
  ABI_DEALLOCATE(kg1_k)

end subroutine mkpwind_k
!!***

!!****f* m_kg/mkkpgcart
!! NAME
!! mkkpgcart
!!
!! FUNCTION
!! Compute all (k+G) vectors (dp, in cartesian coordinates) for given k point,
!! from integer coordinates of G vectors.
!!
!! INPUTS
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  kg(3,npw)=integer coords of planewaves in basis sphere
!!  kpt(3)=k point in terms of recip. translations
!!  nkpg=second dimension of array kpg
!!  npw=number of plane waves in reciprocal space
!!
!! OUTPUT
!!  kpg(npw,3)= (k+G) components
!!
!! PARENTS
!!      m_nonlop_ylm
!!
!! CHILDREN
!!
!! SOURCE


subroutine mkkpgcart(gprimd,kg,kpgcar,kpt,nkpg,npw)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpg,npw
!arrays
 real(dp),intent(in) :: gprimd(3,3)
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(out) :: kpgcar(npw,nkpg)

!Local variables-------------------------------
!scalars
 integer :: ipw,mu
 character(len=500) :: message
!arrays
 real(dp),allocatable :: kpg(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 if (nkpg==0) return

!-- Test nkpg --
 if (nkpg/=3) then
   write(message, '(a,i0)' )' Bad value for nkpg !',nkpg
   MSG_BUG(message)
 end if

!-- Compute (k+G) --
 ABI_ALLOCATE(kpg,(npw,nkpg))
!$OMP PARALLEL DO COLLAPSE(2) &
!$OMP PRIVATE(mu,ipw)
 do ipw=1,npw
   do mu=1,3
     kpg(ipw,mu)=kpt(mu)+dble(kg(mu,ipw))
   end do
 end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO &
!$OMP PRIVATE(ipw)
 do ipw=1,npw
   kpgcar(ipw,1)=kpg(ipw,1)*gprimd(1,1)+kpg(ipw,2)*gprimd(1,2)+kpg(ipw,3)*gprimd(1,3)
   kpgcar(ipw,2)=kpg(ipw,1)*gprimd(2,1)+kpg(ipw,2)*gprimd(2,2)+kpg(ipw,3)*gprimd(2,3)
   kpgcar(ipw,3)=kpg(ipw,1)*gprimd(3,1)+kpg(ipw,2)*gprimd(3,2)+kpg(ipw,3)*gprimd(3,3)
 end do
!$OMP END PARALLEL DO
 ABI_DEALLOCATE(kpg)

 DBG_EXIT("COLL")

end subroutine mkkpgcart
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkkin_metdqdq
!! NAME
!! mkkin_metdqdq
!!
!! FUNCTION
!! Compute elements of the second q-gradient of the metric
!! kinetic energy operator in reciprocal
!! space at given k point wrt cartesian q components.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (DRH, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  effmass=effective mass for electrons (1. in common case)
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  idir= strain perturbation direction
!!  kg(3,npw) = integer coordinates of planewaves in basis sphere.
!!  kpt(3)    = reduced coordinates of k point
!!  npw       = number of plane waves at kpt.
!!  qdir      = direction of the first q-gradient
!!
!! OUTPUT
!!  dqdqkinpw(npw)=d/deps(istr) ( (1/2)*(2 pi)**2 * (k+G)**2 )
!!
!! NOTES
!!
!!  **Since the 2nd derivative w.r.t q-vector is calculated along cartesian
!!    directions, the 1/twopi**2 factor (that in the rest of the code is applied
!!    in the reduced to cartesian derivative conversion process) is here
!!    explicictly included in the formulas.
!!
!!  **Notice that idir=1-9, in contrast to the strain perturbation (idir=1-6),
!!    because this term is not symmetric w.r.t permutations of the two strain
!!    indices.
!!
!!  **A -i factor has been factorized out in all the contributions of the second
!!    q-gradient of the metric Hamiltonian. This is lately included in the contribution
!!    of the corresponing term (T4) to the flexoelectric tensor in dfpt_flexoout.F90
!!
!! PARENTS
!!      m_getgh1c
!!
!! CHILDREN
!!
!! SOURCE

subroutine mkkin_metdqdq(dqdqkinpw,effmass,gprimd,idir,kg,kpt,npw,qdir)

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: idir,npw,qdir
 real(dp),intent(in) :: effmass
!arrays
 integer,intent(in) :: kg(3,npw)
 real(dp),intent(in) :: gprimd(3,3),kpt(3)
 real(dp),intent(out) :: dqdqkinpw(npw)

!Local variables -------------------------
!scalars
 integer :: beta,delta,gamma,ig,ka,kb
 real(dp) :: delbd,delbg,deldg
 real(dp) :: dkinetic,gpk1,gpk2,gpk3,htpi
 real(dp) :: gpkc(3)
!arrays
 integer,save :: idx(18)=(/1,1,2,2,3,3,3,2,3,1,2,1,2,3,1,3,1,2/)

! *********************************************************************

!htpi is (1/2) (2 Pi):
 htpi=0.5_dp*two_pi

 ka=idx(2*idir-1);kb=idx(2*idir)

!For easier formula implementation
 beta=ka
 delta=kb
 gamma=qdir

!Kronecker deltas
 delbd=0.0_dp; delbg=0.0_dp; deldg=0.0_dp
 if (beta==delta) delbd=1.0_dp
 if (beta==gamma) delbg=1.0_dp
 if (delta==gamma) deldg=1.0_dp

 do ig=1,npw
   gpk1=dble(kg(1,ig))+kpt(1)
   gpk2=dble(kg(2,ig))+kpt(2)
   gpk3=dble(kg(3,ig))+kpt(3)

!  Obtain G in cartesian coordinates
   gpkc(1)=gprimd(1,1)*gpk1+gprimd(1,2)*gpk2+gprimd(1,3)*gpk3
   gpkc(2)=gprimd(2,1)*gpk1+gprimd(2,2)*gpk2+gprimd(2,3)*gpk3
   gpkc(3)=gprimd(3,1)*gpk1+gprimd(3,2)*gpk2+gprimd(3,3)*gpk3

   dkinetic=htpi*(2.0_dp*deldg*gpkc(beta)+delbg*gpkc(delta)+ &
 & delbd*gpkc(gamma))

   dqdqkinpw(ig)=dkinetic/effmass
 end do

end subroutine mkkin_metdqdq
!!***

end module m_kg
!!***
