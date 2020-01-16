!!****m* ABINIT/m_nonlop_pl
!! NAME
!!  nonlop_pl
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR, GZ, MT, FF, DRH)
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

module m_nonlop_pl

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_contistr01
 use m_contistr03
 use m_contistr12
 use m_contstr21
 use m_contstr23
 use m_contstr25
 use m_contstr25a
 use m_contstr26
 use m_metstr
 use m_opernl

 use defs_abitypes, only : MPI_type
 use m_geometry,   only : strconv
 use m_kg,         only : ph1d3d
 use m_contract,   only : cont22cso, cont22so, cont24, cont33cso, cont33so, cont35, cont22, cont3, cont13, &
                          metcon, metcon_so, metric_so
 implicit none

 private
!!***

 public :: nonlop_pl
!!***

contains
!!***

!!****f* ABINIT/nonlop_pl
!! NAME
!! nonlop_pl
!!
!! FUNCTION
!! * Compute application of a nonlocal operator Vnl in order to get:
!!    - contracted elements (energy, forces, stresses, ...), if signs=1
!!    - a function in reciprocal space (|out> = Vnl|in>), if signs=2
!!   Operator Vnl, as the following form:
!!    $Vnl=sum_{R,lmn,l''m''n''} {|P_{Rlmn}> Enl^{R}_{lmn,l''m''n''} <P_{Rl''m''n''}|}$
!!   Operator Vnl is -- in the typical case -- the nonlocal potential.
!!   - With norm-conserving pseudopots, $Enl^{R}_{lmn,l''m''n''}$ is the
!!     Kleinmann-Bylander energy $Ekb^{R}_{ln}$.
!!   - The |P_{Rlmn}> are the projector functions.
!! * This routine uses Legendre polynomials Pl to express Vnl.
!!
!! INPUTS
!!  choice: chooses possible output:
!!    choice=1 => a non-local energy contribution
!!          =2 => a gradient with respect to atomic position(s)
!!          =3 => a gradient with respect to strain(s)
!!          =23=> a gradient with respect to atm. pos. and strain(s)
!!          =4 => a gradient and 2nd derivative with respect to atomic pos.
!!          =5 => a gradient with respect to k wavevector
!!          =6 => 2nd derivatives with respect to strain
!!  dimekb1,dimekb2=dimensions of ekb (see ekb)
!!  dimffnlin=second dimension of ffnlin (1+number of derivatives)
!!  dimffnlout=second dimension of ffnlout (1+number of derivatives)
!!  ekb(dimekb1,dimekb2,nspinortot**2)= (Real) Kleinman-Bylander energies (hartree)
!!                                   dimekb1=lmnmax  -  dimekp2=ntypat
!!  ffnlin(npwin,dimffnlin,lmnmax,ntypat)=nonlocal form factors to be used
!!          for the application of the nonlocal operator to the |in> vector
!!  ffnlout(npwout,dimffnlout,lmnmax,ntypat)=nonlocal form factors to be used
!!          for the application of the nonlocal operator to the |out> vector
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!   (bohr^-1)
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2),
!!                        - k point direction in the case (choice=5,signs=2)
!!                        - strain component (1:6) in the case (choice=3,signs=2) or (choice=6,signs=1)
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=ln
!!  istwf_k=option parameter that describes the storage of wfs
!!  kgin(3,npwin)=integer coords of planewaves in basis sphere, for the |in> vector
!!  kgout(3,npwout)=integer coords of planewaves in basis sphere, for the |out> vector
!!  kpgin(npw,npkgin)= (k+G) components and related data, for the |in> vector
!!  kpgout(npw,nkpgout)=(k+G) components and related data, for the |out> vector
!!  kptin(3)=k point in terms of recip. translations, for the |in> vector
!!  kptout(3)=k point in terms of recip. translations, for the |out> vector
!!  lmnmax=max. number of (l,m,n) components over all types of atoms
!!  matblk=dimension of the arrays ph3din and ph3dout
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell
!!  nattyp(ntypat)=number of atoms of each type
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpgin,nkpgout=second sizes of arrays kpgin/kpgout
!!  nloalg(3)=governs the choice of the algorithm for nonlocal operator
!!  nnlout=dimension of enlout: choice=1=>nnlout=1   choice=2=>nnlout=3*natom
!!                              choice=3=>nnlout=6   choice=4=>nnlout=6*natom
!!                              choice=5=>nnlout=1   choice=6=>nnlout=6*(3*natom+6)
!!                              choice=23=>nnlout=6+3*natom
!!  npwin=number of planewaves for given k point, for the |in> vector
!!  npwout=number of planewaves for given k point, for the |out> vector
!!  nspinor=number of spinorial components of the wavefunctions on current proc
!!  nspinortot=total number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in cell
!!  only_SO=flag to calculate only the SO part in nonlop
!!  phkxredin(2,natom)=phase factors exp(2 pi kptin.xred)
!!  phkxredout(2,natom)=phase factors exp(2 pi kptout.xred)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1D structure factors phase information
!!  ph3din(2,npwin,matblk)=3D structure factors, for each atom and plane wave (in)
!!  ph3dout(2,npwout,matblk)=3-dim structure factors, for each atom and plane wave (out)
!!  --- pspso removed in beautification because it was unused ---
!!  pspso(ntypat)=spin-orbit characteristic for each atom type
!!  -------------------------------------------------------------
!!  signs= if 1, get contracted elements (energy, forces, stress, ...)
!!         if 2, applies the non-local operator to a function in reciprocal space
!!  ucvol=unit cell volume (bohr^3)
!!  vectin(2,nspinor*npwin)=input cmplx wavefunction coefficients <G|Cnk>
!!
!! OUTPUT
!!  ==== if (signs==1) ====
!!     enlout(nnlout)= contribution of this state to the nl part
!!                     of the following properties:
!!       if choice=1 : enlout(1)               -> the energy
!!       if choice=2 : enlout(1:3*natom)       -> the forces
!!       if choice=3 : enlout(1:6)             -> the stresses
!!       if choice=23: enlout(1:6+3*natom)     -> the forces and the stresses
!!       if choice=4 : enlout(1:6*natom)       -> the frozen wf part of dynam. matrix
!!       if choice=6 : enlout(1:6*(3*natom+6)) -> the frozen wf part of elastic tensor
!!  ==== if (signs==2) ====
!!     vectout(2,nspinor*npwout)= result of the aplication of the nl operator
!!                                or one of its derivative to the input vect.:
!!       if choice=1 : Vnl |vectin>
!!       if choice=2 : dVnl/d(xred(idir,iatom) |vectin> (xred=reduced atm. pos.)
!!       if choice=3 : dVnl/d(strain(idir)) |vectin>    (symmetric strain =>idir=1...6)
!!       if choice=5 : dVnl/dk(idir) |vectin>           (k wavevector)
!!
!! NOTES
!! In the case signs=1, the array vectout is not used, nor modified
!! so that the same array as vectin can be used as a dummy argument;
!! the same is true for the pairs npwin-npwout, ffnlin-ffnlout,
!! kgin-kgout, ph3din-ph3dout, phkredin-phkxredout).
!!
!! Calculation includes contributions to grads of Etot wrt coord and
!! wrt strains for l=0,1,2,3.
!!
!! WARNINGS
!!  - Warning 1: This routine is in a transient state, during the
!!    time of the implementation of the spin-orbit coupling...
!!    In particular, the OMP parallelisation is still missing,
!!    but it matters here only when nspinor==2.
!!  - Warning 2: the order of atoms is governed by atindx
!!
!! PARENTS
!!      nonlop
!!
!! CHILDREN
!!      cont13,cont22,cont22cso,cont22so,cont24,cont3,cont33cso,cont33so,cont35
!!      contistr01,contistr03,contistr12,contstr21,contstr23,contstr25
!!      contstr25a,contstr26,ddkten,metcon,metcon_so,metric_so,metstr,opernl2
!!      opernl3,opernl4a,opernl4b,ph1d3d,scalewf_nonlop,strconv,strsocv,trace2
!!      xmpi_sum
!!
!! SOURCE

subroutine nonlop_pl(choice,dimekb1,dimekb2,dimffnlin,dimffnlout,ekb,enlout,&
&                     ffnlin,ffnlout,gmet,gprimd,idir,indlmn,istwf_k,kgin,kgout,kpgin,kpgout,&
&                     kptin,kptout,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,&
&                     natom,nattyp,ngfft,nkpgin,nkpgout,nloalg,npwin,npwout,nspinor,nspinortot,&
&                     ntypat,only_SO,phkxredin,phkxredout,ph1d,ph3din,ph3dout,signs,&
&                     ucvol,vectin,vectout)

!Arguments ------------------------------------
!This type is defined in defs_mpi
!The (inout) classification below is misleading; mpi_enreg is temporarily
! changed but reset to its initial condition before exiting.
!scalars
 integer,intent(in) :: choice,dimekb1,dimekb2,dimffnlin,dimffnlout,idir,istwf_k
 integer,intent(in) :: lmnmax,matblk,mgfft,mpsang,mpssoang,natom,nkpgin,nkpgout
 integer,intent(in) :: npwin,npwout,nspinor,nspinortot,ntypat,only_SO,signs
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),kgin(3,npwin),kgout(3,npwout)
 integer,intent(in) :: nattyp(ntypat),ngfft(18),nloalg(3) !,pspso(ntypat) UNUSED
 real(dp),intent(in) :: ekb(dimekb1,dimekb2,nspinortot**2)
 real(dp),intent(in) :: ffnlin(npwin,dimffnlin,lmnmax,ntypat)
 real(dp),intent(in) :: ffnlout(npwout,dimffnlout,lmnmax,ntypat),gmet(3,3)
 real(dp),intent(in) :: gprimd(3,3),kpgin(npwin,nkpgin),kpgout(npwout,nkpgout)
!real(dp),intent(in) :: kptin(3),kptout(3),ph1d(2,3*(2*mgfft+1)*natom) !vz_d
 real(dp),intent(in) :: kptin(3),kptout(3) !vz_d
 real(dp),intent(in) :: ph1d(2,*) !vz_d
 real(dp),intent(in) :: phkxredin(2,natom),phkxredout(2,natom)
 real(dp),intent(inout) :: ph3din(2,npwin,matblk),ph3dout(2,npwout,matblk)
 real(dp),intent(inout) :: vectin(:,:)
 real(dp),intent(out) :: enlout(:) !vz_i
 real(dp),intent(inout) :: vectout(:,:) !vz_i

!Local variables-------------------------------
!mlang is the maximum number of different angular momenta
!(mlang=4 means s,p,d,f)
! Note : in a future version, one should adjust mlang to mpsang.
!mlang2 is the maximum number of unique tensor components for a tensor
!of rank (mlang-1) with index range 1-3
!mlang3 is the maximum number of unique tensor components summed over
!all tensors of rank 0 through mlang-1.
!mlang4 is the total number of additional unique tensor components
!related to strain gradients, ranks 2 through mlang+1.
!mlang6 is the total number of additional unique tensor components
!related to strain 2nd derivaives, ranks 4 through mlang+3.
!mlang1 is the total number of certain additional unique tensor components
!related to internal strain, ranks 1 through mlang
!mlang5 is the total number of other additional unique tensor components
!related to internal strain, ranks 1 through mlang
!scalars
 integer,parameter :: mlang=4
! MG: I tried to use parameters instead of saved variables but [tutorespfn][trf2_1] gets stuck on milou_g95_snofbfpe
 integer,save :: mlang1=((mlang+1)*(mlang+2)*(mlang+3))/6-1
 !integer,save :: mlang2=(mlang*(mlang+1))/2 ! Unused
 integer,save :: mlang3=(mlang*(mlang+1)*(mlang+2))/6
 integer,save :: mlang4=((mlang+2)*(mlang+3)*(mlang+4))/6-4
 integer,save :: mlang5=((mlang+3)*(mlang+4)*(mlang+5))/6-10
 integer,save :: mlang6=((mlang+4)*(mlang+5)*(mlang+6))/6-20
 integer :: compact,ia,ia1,ia2,ia3,ia4,ia5,ierr,iest,ig,ii,ilang,ilang2,ilmn
 integer :: iln,iln0,indx,iproj,ipsang,ishift,isp,ispin,ispinor,ispinor_index,ispinp
 integer :: istr,istr1,istr2,iterm,itypat,jj,jjk,jjs,jjs1,jjs2,jjs3,jjs4,jjstr,jspin
 integer :: mincat,mproj,mu,mumax,n1,n2,n3,ndgxdt,ndgxdtfac,nincat,nlang
 integer :: nproj,nspinso,rank
 integer :: sign,spaceComm,  isft
 real(dp) :: e2nl,e2nldd,enlk
 character(len=500) :: message
!arrays
 integer,allocatable :: indlmn_s(:,:,:),jproj(:)
 real(dp) :: amet(2,3,3,2,2),amet_lo(3,3),e2nl_tmp(6),eisnl(3),rank2(6)
 real(dp) :: rank2c(2,6),strsnl(6),strsnl_out(6),strsso(6,3),strssoc(6),trace(2)!,tsec(2)
 real(dp),allocatable :: d2gxdis(:,:,:,:,:),d2gxdis_s(:,:,:,:)
 real(dp),allocatable :: d2gxds2(:,:,:,:,:),d2gxds2_s(:,:,:,:)
 real(dp),allocatable :: dgxdis(:,:,:,:,:),dgxdis_s(:,:,:,:),dgxds(:,:,:,:,:)
 real(dp),allocatable :: dgxds_s(:,:,:,:),dgxdsfac(:,:,:,:,:)
 real(dp),allocatable :: dgxdt(:,:,:,:,:,:),dgxdt_s(:,:,:,:,:)
 real(dp),allocatable :: dgxdtfac(:,:,:,:,:),ekb_s(:,:),gxa(:,:,:,:,:)
 real(dp),allocatable :: gxa_s(:,:,:,:),gxafac(:,:,:,:),pauli(:,:,:,:)
 real(dp),allocatable :: temp(:,:),tmpfac(:,:),vectin_s(:,:),vectout_s(:,:)
 real(dp),allocatable :: wt(:,:)

! **********************************************************************

 ABI_UNUSED(mgfft)

!Test: spin orbit not allowed for choice=5,6
 if (nspinortot==2 .and. choice==6 ) then
   message = 'nonlop_pl: For nspinortot=2, choice=6 is not yet allowed.'
   MSG_BUG(message)
 end if

 if ((choice<1 .or. choice>6) .and. choice/=23 ) then
   write(message,'(a,i0)')'  Does not presently support this choice=',choice
   MSG_BUG(message)
 end if

!Test: choice 51 and 52 only allowed with nonlop_ylm
!JWZ, 01-Sep-08
 if (choice==51 .or. choice==52) then
   message = 'nonlop_pl: choice 51 or 52 is not yet allowed.'
   MSG_BUG(message)
 end if

!Define dimension of work arrays.
 mincat=min(NLO_MINCAT,maxval(nattyp))
 mproj=maxval(indlmn(3,:,:))
 ABI_ALLOCATE(temp,(2,mlang4))
 ABI_ALLOCATE(tmpfac,(2,mlang4))
 ABI_ALLOCATE(wt,(mlang,mproj))
 ABI_ALLOCATE(jproj,(mlang))
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 ABI_ALLOCATE(ekb_s,(mlang,mproj))
 ABI_ALLOCATE(indlmn_s,(6,lmnmax,ntypat))

!Eventually compute the spin-orbit metric tensor:
 if (mpssoang>mpsang) then
   ABI_ALLOCATE(pauli,(2,2,2,3))
   call metric_so(amet,gprimd,pauli)
 end if

!Allocate array gxa (contains projected scalars).
 ABI_ALLOCATE(gxa,(2,mlang3,mincat,mproj,nspinortot))
 if(nspinor==2)  then
   ABI_ALLOCATE(gxa_s,(2,mlang3,mincat,mproj))
 else
   ABI_ALLOCATE(gxa_s,(0,0,0,0))
 end if

 ABI_ALLOCATE(gxafac,(2,mlang3,mincat,mproj))
 gxa(:,:,:,:,:)=zero

!If choice==2 : first-order atomic displacements
!If signs==2, only one direction is considered
!If signs==1, the three directions are considered
!If choice==4 and signs==1 : second-order atomic displacements,
!the nine components are considered
!If choice==5 and signs==2 : ddk
!component 1 -> from ffnl(:,2,...)
!component 2 -> from ffnl(:,1,...) (too much space is booked for this
!case, since the number of angular momenta is smaller than mlang3, but it is easier)
 ndgxdt=0
 if(signs==2 .and. choice==2) ndgxdt=1
 if(signs==1 .and. (choice==2.or.choice==23)) ndgxdt=3
 if(choice==4) ndgxdt=9
 if(choice==5) ndgxdt=2
!Allocate dgxdt (contains derivatives of gxa with respect to atomic displacements or ddk).
 ABI_ALLOCATE(dgxdt,(2,ndgxdt,mlang3,mincat,mproj,nspinortot))
 dgxdt(:,:,:,:,:,:)=zero
 if(nspinor==2)then
   ABI_ALLOCATE(dgxdt_s,(2,ndgxdt,mlang3,mincat,mproj))
   dgxdt_s(:,:,:,:,:)=zero
 else
   ABI_ALLOCATE(dgxdt_s,(0,0,0,0,0))
 end if
 ndgxdtfac=0
 if(signs==2 .and. choice==2) ndgxdtfac=1
 if(choice==4) ndgxdtfac=3
 if(choice==5) ndgxdtfac=2
 ABI_ALLOCATE(dgxdtfac,(2,ndgxdtfac,mlang3,mincat,mproj))

!Allocate dgxds (contains derivatives of gxa with respect to strains).
 ABI_ALLOCATE(dgxds,(2,mlang4,mincat,mproj,nspinor))
 dgxds(:,:,:,:,:)=zero
 ABI_ALLOCATE(dgxdsfac,(2,mlang4,mincat,mproj,nspinor))
 if(choice==6) then
   ABI_ALLOCATE(dgxdis,(2,mlang1,mincat,mproj,nspinor))
   ABI_ALLOCATE(d2gxdis,(2,mlang5,mincat,mproj,nspinor))
   ABI_ALLOCATE(d2gxds2,(2,mlang6,mincat,mproj,nspinor))
 else
   ABI_ALLOCATE(dgxdis ,(0,0,0,0,0))
   ABI_ALLOCATE(d2gxdis,(0,0,0,0,0))
   ABI_ALLOCATE(d2gxds2,(0,0,0,0,0))
 end if
 ABI_ALLOCATE(dgxds_s  ,(0,0,0,0))
 ABI_ALLOCATE(dgxdis_s ,(0,0,0,0))
 ABI_ALLOCATE(d2gxdis_s,(0,0,0,0))
 ABI_ALLOCATE(d2gxds2_s,(0,0,0,0))
 if(nspinor==2)then
   ABI_DEALLOCATE(dgxds_s)
   ABI_ALLOCATE(dgxds_s,(2,mlang4,mincat,mproj))
   dgxds_s(:,:,:,:)=zero
   if(choice==6) then
     ABI_DEALLOCATE(dgxdis_s)
     ABI_DEALLOCATE(d2gxdis_s)
     ABI_DEALLOCATE(d2gxds2_s)
     ABI_ALLOCATE(dgxdis_s,(2,mlang1,mincat,mproj))
     ABI_ALLOCATE(d2gxdis_s,(2,mlang5,mincat,mproj))
     ABI_ALLOCATE(d2gxds2_s,(2,mlang6,mincat,mproj))
   else
   end if
 end if

!Zero out some arrays
 if(choice==2 .or. choice==4 .or. choice==5 .or. choice==6 .or. choice==23) enlout(:)=0.0d0

 if(signs==2) vectout(:,:)=zero

 if(choice==3.or.choice==23) then
   strsnl(:)=zero
   if(mpssoang>mpsang) strsso(:,:)=zero
 end if
 enlk=zero

!In the case vectin is a spinor, split its second part.
!Also, eventually take into account the storage format of the wavefunction
!(the original vector will be restored before leaving the routine,
!except for the vectin(2,1) component with istwf_k==2,
!that should vanish)
!In sequential, treat the second spinor part first
 if (nspinor==2)then
   ABI_ALLOCATE(vectin_s,(2,npwin))
   ABI_ALLOCATE(vectout_s,(2,npwout))

   isft = npwin;if (mpi_enreg%nproc_spinor>1) isft=0

!  Initialize it
!$OMP PARALLEL DO
   do ig=1,npwin
     vectin_s(1,ig)=vectin(1,ig+isft)
     vectin_s(2,ig)=vectin(2,ig+isft)
   end do

!  Take into account the storage
   if(istwf_k/=1)then
     call scalewf_nonlop(istwf_k,mpi_enreg,npwin,1,vectin_s)
   end if
 end if ! nspinortot==2

!Treat the first spinor part now
 if(istwf_k/=1) then
   call scalewf_nonlop(istwf_k,mpi_enreg,npwin,1,vectin)
 end if


!Big loop on atom types.
 ia1=1
 do itypat=1,ntypat

!  Get atom loop indices for different types:
   ia2=ia1+nattyp(itypat)-1

!  Cut the sum on different atoms in blocks, to allow memory saving.
!  Inner summations on atoms will be done from ia3 to ia4.
!  Note that the maximum range from ia3 to ia4 is mincat (maximum
!  increment of atoms).
   do ia3=ia1,ia2,mincat
     ia4=min(ia2,ia3+mincat-1)
!    Give the increment of number of atoms in this subset.
     nincat=ia4-ia3+1

!    Prepare the phase factors for the atoms between ia3 and ia4 :
!    For nloalg(2)<=0, they were not prepared previously,it is needed to
!    compute them again.
     if(nloalg(2)<=0)then
!      For nloalg(2)==0, it is needed to compute the phase factors.
       if(mincat>matblk)then
         write(message,'(a,a,a,i4,a,i4,a)')&
&         '  With nloc_mem<=0, mincat must be less than matblk.',ch10,&
&         '  Their value is ',mincat,' and ',matblk,'.'
         MSG_BUG(message)
       end if
       call ph1d3d(ia3,ia4,kgin,matblk,natom,npwin,n1,n2,n3,phkxredin,ph1d,ph3din)
     end if

!    Here begins the different treatment for the scalar-relativistic
!    part(s) and the spin-orbit part.
!    Loop on ispinor : 1 for scalar-Relativistic, 2 for spin-orbit
     nspinso=1;if (mpssoang>mpsang) nspinso=2

     ! Change nspinso if collinear run or if nspinor == 2 and SOC is not wanted.
     ! TODO: The last check requires pspso
     if (nspinortot == 1) nspinso = 1

     do ispinor=1,nspinso
       ispinor_index=ispinor
       if (mpi_enreg%paral_spinor==1) ispinor_index=mpi_enreg%me_spinor+1

!
!      mjv 25 6 2008: if only_SO == 1 or 2 skip the scalar relativistic terms
!      only output the spin orbit ones
!
       if (ispinor==1 .and. only_SO>0) cycle

!      Select scalar-relativistic or spin-orbit KB-energies:
       ekb_s(:,:)=zero ; wt(:,:)=zero
!      Loop over (l,n) values (loop over l,m,n and test on l,n)
       iln0=0 ; jproj(:)=0 ; nlang=0
       indlmn_s(:,:,itypat)=0
       do ilmn=1,lmnmax
         if(ispinor/=indlmn(6,ilmn,itypat))cycle
         indlmn_s(:,ilmn,itypat)=indlmn(:,ilmn,itypat)
         iln =indlmn(5,ilmn,itypat)
         if (iln>iln0) then
           iln0=iln
           ipsang=indlmn(1,ilmn,itypat)+1
!          DEBUG
!          write(std_out,*)' nonlop : ipsang,ilmn,itypat,ispinor=',ipsang,ilmn,itypat,ispinor
!          ENDDEBUG
           iproj=indlmn(3,ilmn,itypat)
!          This shift is not needed anymore
!          if (ispinor==2) ipsang=indlmn(1,ilmn,itypat)-mpsang+2
           ekb_s(ipsang,iproj)=ekb(iln,itypat,ispinor)
           wt(ipsang,iproj)=4.d0*pi/ucvol*dble(2*ipsang-1)*ekb_s(ipsang,iproj)
!
!          mjv 27 6 2008: if only_SO == 2 remove the factor of l in the operator
!
           if (only_SO == 2) then
             wt(ipsang,iproj)=4.d0*pi/ucvol*ekb_s(ipsang,iproj)
           end if
           jproj(ipsang)=max(jproj(ipsang),iproj)
           if(iproj>0)nlang=max(nlang,ipsang)
         end if
       end do ! ilmn


!      If nlang is not 0, then some non-local part is to be applied for that type of atom.
       if (nlang/=0) then
!        Operate with the non-local potential on the wavefunction, in order
!        to get projected quantities. Call different routines opernl2,
!        opernl3, opernl4x which corresponds to different writings of the
!        same numerical operations. There is still optimisation left for
!        the case istwf_k/=1 (up to a factor 2 in speed).
!        call timab(74+choice,1,tsec)
         sign=1
         gxa(:,:,:,:,:)=zero
         dgxdt(:,:,:,:,:,:)=zero
         dgxds(:,:,:,:,:)=zero

!        Only the first spinorial component of vectin is taken into account first
         ispin=1;if (mpi_enreg%paral_spinor==1) ispin=ispinor_index
         if(nloalg(1)==2) then
           call opernl2(choice,dgxdis,dgxds,d2gxdis,d2gxds2,dgxdt,&
&           ffnlin,gmet,gxa,ia3,idir,indlmn_s,ispinor,istwf_k,itypat,&
&           jproj,kgin,kpgin,kptin,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&
&           mlang5,mlang6,mproj,ndgxdt,dimffnlin,nincat,nkpgin,nlang,nloalg,npwin,&
&           ntypat,ph3din,sign,vectin)
         else if(nloalg(1)==3) then
           call opernl3(choice,dgxdis,dgxds,d2gxdis,d2gxds2,dgxdt,&
&           ffnlin,gmet,gxa,ia3,idir,indlmn_s,ispinor,istwf_k,itypat,&
&           jproj,kgin,kpgin,kptin,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&
&           mlang5,mlang6,mproj,ndgxdt,dimffnlin,nincat,nkpgin,nlang,nloalg,npwin,&
&           ntypat,ph3din,sign,vectin)
         else if(nloalg(1)==4) then
           call opernl4a(choice,dgxdis,dgxds,d2gxdis,d2gxds2,dgxdt,&
&           ffnlin,gmet,gxa,ia3,idir,indlmn_s,ispinor,istwf_k,itypat,&
&           jproj,kgin,kpgin,kptin,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&
&           mlang5,mlang6,mproj,ndgxdt,dimffnlin,nincat,nkpgin,nlang,nloalg,npwin,&
&           ntypat,ph3din,vectin)
         end if
!        This duplication of the opernl calls is needed to avoid copying
!        vectin, with a detrimental effect on speed.
         if (nspinor==2)then
           ispin=2
           if(nloalg(1)==2) then
             call opernl2(choice,dgxdis_s,dgxds_s,d2gxdis_s,d2gxds2_s,dgxdt_s,&
&             ffnlin,gmet,gxa_s,ia3,idir,indlmn_s,ispinor,istwf_k,itypat,&
&             jproj,kgin,kpgin,kptin,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&
&             mlang5,mlang6,mproj,ndgxdt,dimffnlin,nincat,nkpgin,nlang,nloalg,npwin,&
&             ntypat,ph3din,sign,vectin_s)
           else if(nloalg(1)==3) then
             call opernl3(choice,dgxdis_s,dgxds_s,d2gxdis_s,d2gxds2_s,dgxdt_s,&
&             ffnlin,gmet,gxa_s,ia3,idir,indlmn_s,ispinor,istwf_k,itypat,&
&             jproj,kgin,kpgin,kptin,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&
&             mlang5,mlang6,mproj,ndgxdt,dimffnlin,nincat,nkpgin,nlang,nloalg,npwin,&
&             ntypat,ph3din,sign,vectin_s)
           else if(nloalg(1)==4) then
             call opernl4a(choice,dgxdis_s,dgxds_s,d2gxdis_s,d2gxds2_s,dgxdt_s,&
&             ffnlin,gmet,gxa_s,ia3,idir,indlmn_s,ispinor,istwf_k,itypat,&
&             jproj,kgin,kpgin,kptin,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&
&             mlang5,mlang6,mproj,ndgxdt,dimffnlin,nincat,nkpgin,nlang,nloalg,npwin,&
&             ntypat,ph3din,vectin_s)
           end if
           dgxds(:,:,:,:,ispin)=dgxds_s(:,:,:,:)
           dgxdt(:,:,:,:,:,ispin)=dgxdt_s(:,:,:,:,:)
           gxa(:,:,:,:,ispin)=gxa_s(:,:,:,:)
         end if

!        Parallelism stuff
         spaceComm=mpi_enreg%comm_fft
         call xmpi_sum(dgxds,spaceComm,ierr)
         if (mpi_enreg%paral_spinor==1) then
           spaceComm=mpi_enreg%comm_spinorfft
           jspin=3-ispin
           gxa(:,:,:,:,ispin)=gxa(:,:,:,:,1)
           gxa(:,:,:,:,jspin)=zero
           if ( ndgxdt>0)then
             dgxdt(:,:,:,:,:,ispin)=dgxdt(:,:,:,:,:,1)
             dgxdt(:,:,:,:,:,jspin)=zero
           end if
         end if

         call xmpi_sum(gxa,spaceComm,ierr)
         if(ndgxdt>0) then
           call xmpi_sum(dgxdt,spaceComm,ierr)
         end if

!        XG030513 : MPIWF, at this place, one should perform the reduction
!        and spread of data of gxa, dgxdt and dgxds

!        Loop over spins:
         do isp=1,nspinor
           ispin=isp;if (mpi_enreg%paral_spinor==1) ispin=ispinor_index

!          Perform contractions for the various tensors (d)gx?, producing the
!          contracted tensors (d)gx?fac to be passed back to opernl:
           do ia=1,nincat
             do ilang=1,nlang
               nproj=jproj(ilang)
               if(nproj/=0) then
                 ilang2=(ilang*(ilang+1))/2
                 do iproj=1,nproj

!                  The rank of the tensor gxa equals l:
                   rank=ilang-1
!                  jjs gives the starting address of the relevant components
                   jjs=1+((ilang-1)*ilang*(ilang+1))/6
                   if (ilang>4) then
                     write(message,'(a,i0)')' ilang must fall in range [1..4] but value is ',ilang
                     MSG_BUG(message)
                   end if

!                  Metric & spinorial contraction from gxa to gxafac. The treatment
!                  is different for the scalar-relativistic and spin-orbit parts.
                   if(ispinor==1) then
!                    ------ Scalar-Relativistic ------
                     temp(:,1:((rank+1)*(rank+2))/2)= &
&                     gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin)
                     call metcon(rank,gmet,temp,tmpfac)
                     gxafac(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)= &
&                     wt(ilang,iproj)*tmpfac(:,1:((rank+1)*(rank+2))/2)
                   else
!                    ------ Spin-orbit ------
                     gxafac(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)=zero
!                    Contraction over spins:
                     do ispinp=1,nspinortot
!                      => Imaginary part (multiplying by i, then by the Im of amet):
                       temp(1,1:((rank+1)*(rank+2))/2)= &
&                       -gxa(2,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispinp)
                       temp(2,1:((rank+1)*(rank+2))/2)= &
&                       gxa(1,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispinp)
                       amet_lo(:,:)=amet(2,:,:,ispin,ispinp)
                       call metcon_so(rank,gmet,amet_lo,temp,tmpfac)
                       gxafac(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)= &
&                       gxafac(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)+ &
&                       wt(ilang,iproj)*tmpfac(:,1:((rank+1)*(rank+2))/2)
!                      => Real part:
                       temp(:,1:((rank+1)*(rank+2))/2)= &
&                       gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispinp)
                       amet_lo(:,:)=amet(1,:,:,ispin,ispinp)
                       call metcon_so(rank,gmet,amet_lo,temp,tmpfac)
                       gxafac(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)= &
&                       gxafac(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)+ &
&                       wt(ilang,iproj)*tmpfac(:,1:((rank+1)*(rank+2))/2)
                     end do
                   end if

!                  Eventual tensorial compaction/decompaction for ddk
!                  perturbation:
                   if(choice==5 .and. ilang>=2) then
                     jjk=1+((ilang-2)*(ilang-1)*ilang)/6
                     compact=-1
                     temp(:,1:(rank*(rank+1))/2)= &
&                     dgxdt(:,2,jjk:jjk-1+(rank*(rank+1))/2,ia,iproj,ispin)
                     call ddkten(compact,idir,rank,temp,tmpfac)
                     dgxdt(:,1,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin)= &
&                     dgxdt(:,1,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin)&
&                     +tmpfac(:,1:((rank+1)*(rank+2))/2)
                     compact=1
                     tmpfac(:,1:((rank+1)*(rank+2))/2)= &
&                     gxafac(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)
                     call ddkten(compact,idir,rank,temp,tmpfac)
                     dgxdtfac(:,2,jjk:jjk-1+(rank*(rank+1))/2,ia,iproj)= &
&                     temp(:,1:(rank*(rank+1))/2)
                   end if

!                  Section for strain perturbation
!                  no spin-orbit yet

                   if(choice==3 .and. signs==2) then
                     istr=idir
                     if(ispinor==1) then
!                      ------ Scalar-Relativistic ------
!                      jjstr is the starting address for dgxds and dgxdsfac
                       jjstr=-3+((ilang+1)*(ilang+2)*(ilang+3))/6
!                      diagonal enlk contribution
!                      note sign change (12/05/02)
                       if(istr>3) then
                         gxafac(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)=zero
                       else
                         gxafac(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)=&
&                         -gxafac(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)
                       end if
                       dgxdsfac(:,jjstr:jjstr-1+((rank+3)*(rank+4))/2,ia,iproj,isp)=zero
                       iterm=1
                       temp(:,1:((rank+1)*(rank+2))/2)= &
&                       gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin)
                       call metstr(istr,rank,iterm,gmet,gprimd,temp,tmpfac)
                       dgxdsfac(:,jjstr:jjstr-1+((rank+3)*(rank+4))/2,ia,iproj,isp)= &
&                       wt(ilang,iproj)*tmpfac(:,1:((rank+3)*(rank+4))/2)
                       iterm=2
                       temp(:,1:((rank+1)*(rank+2))/2)= &
&                       gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin)
                       call metstr(istr,rank,iterm,gmet,gprimd,temp,tmpfac)
                       gxafac(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)= &
&                       +gxafac(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)+ &
&                       wt(ilang,iproj)*tmpfac(:,1:((rank+1)*(rank+2))/2)
                       iterm=3
                       temp(:,1:((rank+3)*(rank+4))/2)= &
                       dgxds(:,jjstr:jjstr-1+((rank+3)*(rank+4))/2,ia,iproj,isp)
                       call metstr(istr,rank,iterm,gmet,gprimd,temp,tmpfac)
                       gxafac(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)= &
&                       gxafac(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)+ &
&                       wt(ilang,iproj)*tmpfac(:,1:((rank+1)*(rank+2))/2)
                     end if
!                    end section for strain perturbation
                   end if

!                  Eventual metric & spinorial contraction from dgxdt to dgxdtfac.
!                  either for the dynamical matrix, or for the application of the
!                  gradient of the operator. The treatment is different for the
!                  scalar-relativistic and spin-orbit parts.
                   if ((choice==2.and.signs==2).or. &
&                   (choice==5.and.signs==2).or. &
&                   (choice==4)) then
                     mumax=ndgxdtfac;if (choice==5) mumax=1
                     do mu=1,mumax
                       if(ispinor==1) then
!                        ------ Scalar-Relativistic ------
                         temp(:,1:((rank+1)*(rank+2))/2)= &
&                         dgxdt(:,mu,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin)
                         call metcon(rank,gmet,temp,tmpfac)
                         dgxdtfac(:,mu,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)=&
&                         wt(ilang,iproj)*tmpfac(:,1:((rank+1)*(rank+2))/2)
                       else
!                        ------ Spin-orbit ------
                         dgxdtfac(:,mu,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)=zero
!                        Contraction over spins:
                         do ispinp=1,nspinortot
!                          => Imaginary part (multiplying by i, then by the Im of amet):
                           temp(1,1:((rank+1)*(rank+2))/2)= &
&                           -dgxdt(2,mu,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispinp)
                           temp(2,1:((rank+1)*(rank+2))/2)= &
&                           dgxdt(1,mu,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispinp)
                           amet_lo(:,:)=amet(2,:,:,ispin,ispinp)
                           call metcon_so(rank,gmet,amet_lo,temp,tmpfac)
                           dgxdtfac(:,mu,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)=&
&                           dgxdtfac(:,mu,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)+&
&                           wt(ilang,iproj)*tmpfac(:,1:((rank+1)*(rank+2))/2)
!                          => Real part:
                           temp(:,1:((rank+1)*(rank+2))/2)= &
&                           dgxdt(:,mu,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispinp)
                           amet_lo(:,:)=amet(1,:,:,ispin,ispinp)
                           call metcon_so(rank,gmet,amet_lo,temp,tmpfac)
                           dgxdtfac(:,mu,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)=&
&                           dgxdtfac(:,mu,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj)+&
&                           wt(ilang,iproj)*tmpfac(:,1:((rank+1)*(rank+2))/2)
                         end do
                       end if
                     end do
                   end if


!                  ----  Accumulate the nonlocal energy.
                   do ii=1,ilang2
                     jj=ii-1+jjs
                     enlk=enlk+(gxafac(1,jj,ia,iproj)*gxa(1,jj,ia,iproj,ispin)&
&                     +gxafac(2,jj,ia,iproj)*gxa(2,jj,ia,iproj,ispin) )
                   end do

!                  ----  Accumulate force contributions if requested.
!                  Note that the contraction of gxa with dgxdt can be done like
!                  a cartesian dot product now because the symmetrically-related
!                  terms are accounted for with weights in gxa.
                   if ((choice==2.or.choice==23) .and. signs==1) then
                     ishift=0;if (choice==23) ishift=6
                     ia5=ia+ia3-1
                     do ii=1,ilang2
                       jj=ii-1+jjs
                       do mu=1,3
!                        (includes also factor of 2 from "2*Re[stuff]")
                         indx=mu+3*(ia5-1)+ishift
                         enlout(indx)=enlout(indx)+two*&
&                         ( gxafac(1,jj,ia,iproj)*dgxdt(1,mu,jj,ia,iproj,ispin)&
&                         +gxafac(2,jj,ia,iproj)*dgxdt(2,mu,jj,ia,iproj,ispin))
                       end do
                     end do
                   end if

!                  ----  Accumulate stress tensor contributions if requested.
                   if ((choice==3.or.choice==23).and.signs==1) then
!                    1- Compute contractions involving gxa and dgxds:
!                    --- Same formula for Scalar-relativistic and Spin-orbit ---
                     if (ilang==1) then
                       do ii=1,6
                         rank2(ii)=2.d0*&
&                         (gxafac(1,1,ia,iproj)*dgxds(1,ii,ia,iproj,isp)+ &
&                         gxafac(2,1,ia,iproj)*dgxds(2,ii,ia,iproj,isp) )
                       end do
                     else if (ilang==2) then
                       call cont13(gxafac(:,jjs:jjs+2,ia,iproj),&
&                       dgxds(:, 7:16,ia,iproj,isp),rank2)
                     else if (ilang==3) then
                       call cont24(gxafac(:,jjs:jjs+5,ia,iproj),&
&                       dgxds(:,17:31,ia,iproj,isp),rank2)
                     else if (ilang==4) then
                       call cont35(gxafac(:,jjs:jjs+9,ia,iproj),&
&                       dgxds(:,32:52,ia,iproj,isp),rank2)
                     end if
!                    In all cases add rank2 term into stress tensor
                     strsnl(:)=strsnl(:)-rank2(:)
!                    2- Compute contractions involving gxa and gxa:
                     if(ispinor==1) then
!                      2a ------ Scalar-Relativistic ------
                       if (ilang==2) then
                         strsnl(1)=strsnl(1)-wt(ilang,iproj)*2.d0*&
&                         (gxa(1,jjs  ,ia,iproj,ispin)*gxa(1,jjs  ,ia,iproj,ispin)+&
&                         gxa(2,jjs  ,ia,iproj,ispin)*gxa(2,jjs  ,ia,iproj,ispin))
                         strsnl(2)=strsnl(2)-wt(ilang,iproj)*2.d0*&
&                         (gxa(1,jjs+1,ia,iproj,ispin)*gxa(1,jjs+1,ia,iproj,ispin)+&
&                         gxa(2,jjs+1,ia,iproj,ispin)*gxa(2,jjs+1,ia,iproj,ispin))
                         strsnl(3)=strsnl(3)-wt(ilang,iproj)*2.d0*&
&                         (gxa(1,jjs+2,ia,iproj,ispin)*gxa(1,jjs+2,ia,iproj,ispin)+&
&                         gxa(2,jjs+2,ia,iproj,ispin)*gxa(2,jjs+2,ia,iproj,ispin))
                         strsnl(4)=strsnl(4)-wt(ilang,iproj)*2.d0*&
&                         (gxa(1,jjs+2,ia,iproj,ispin)*gxa(1,jjs+1,ia,iproj,ispin)+&
&                         gxa(2,jjs+2,ia,iproj,ispin)*gxa(2,jjs+1,ia,iproj,ispin))
                         strsnl(5)=strsnl(5)-wt(ilang,iproj)*2.d0*&
&                         (gxa(1,jjs+2,ia,iproj,ispin)*gxa(1,jjs  ,ia,iproj,ispin)+&
&                         gxa(2,jjs+2,ia,iproj,ispin)*gxa(2,jjs  ,ia,iproj,ispin))
                         strsnl(6)=strsnl(6)-wt(ilang,iproj)*2.d0*&
&                         (gxa(1,jjs+1,ia,iproj,ispin)*gxa(1,jjs  ,ia,iproj,ispin)+&
&                         gxa(2,jjs+1,ia,iproj,ispin)*gxa(2,jjs  ,ia,iproj,ispin))
                       else if (ilang==3) then
                         call trace2(gxa(:,jjs:jjs+5,ia,iproj,ispin),gmet,trace)
                         call cont22(gxa(:,jjs:jjs+5,ia,iproj,ispin),gmet,rank2)
                         do ii=1,6
                           strsnl(ii)=strsnl(ii)+wt(ilang,iproj)*&
&                           (2.d0*(trace(1)*gxa(1,jjs+ii-1,ia,iproj,ispin)+&
&                           trace(2)*gxa(2,jjs+ii-1,ia,iproj,ispin))-3.d0*rank2(ii))
                         end do
                       else if (ilang==4) then
                         call cont3(gxa(:,jjs:jjs+9,ia,iproj,ispin),gmet,rank2)
                         strsnl(:)=strsnl(:)-wt(ilang,iproj)*rank2(:)
                       end if
                     else
!                      2b ------ Spin-orbit ------
                       do ispinp=1,nspinortot
                         if (ilang==3) then
                           call cont22so(gxa(:,jjs:jjs+5,ia,iproj,ispin),&
&                           gxa(:,jjs:jjs+5,ia,iproj,ispinp),&
&                           amet(:,:,:,ispin,ispinp),rank2)
                           strsnl(:)=strsnl(:)-wt(ilang,iproj)*3.d0*rank2(:)
                         else if (ilang==4) then
                           call cont33so(gxa(:,jjs:jjs+9,ia,iproj,ispin),&
&                           gxa(:,jjs:jjs+9,ia,iproj,ispinp),&
&                           gmet,amet(:,:,:,ispin,ispinp),rank2)
                           strsnl(:)=strsnl(:)-wt(ilang,iproj)*rank2(:)
                         end if
                       end do
                     end if
!                    3- Compute contractions involving gxa and gxa due to
!                    gradients of antisymmetric tensor (amet):
!                    --- Only in case of Spin-orbit ---
                     if(ispinor==2) then
                       do ispinp=1,nspinortot
!                        Be carefull: no need to compute rank2c(:,1:3) !
                         if (ilang==2) then
                           rank2c(1,4)=gxa(1,jjs+2,ia,iproj,ispin)*gxa(1,jjs+1,ia,iproj,ispinp)&
&                           +gxa(2,jjs+2,ia,iproj,ispin)*gxa(2,jjs+1,ia,iproj,ispinp)
                           rank2c(2,4)=gxa(1,jjs+2,ia,iproj,ispin)*gxa(2,jjs+1,ia,iproj,ispinp)&
&                           -gxa(2,jjs+2,ia,iproj,ispin)*gxa(1,jjs+1,ia,iproj,ispinp)
                           rank2c(1,5)=gxa(1,jjs+2,ia,iproj,ispin)*gxa(1,jjs  ,ia,iproj,ispinp)&
&                           +gxa(2,jjs+2,ia,iproj,ispin)*gxa(2,jjs  ,ia,iproj,ispinp)
                           rank2c(2,5)=gxa(1,jjs+2,ia,iproj,ispin)*gxa(2,jjs  ,ia,iproj,ispinp)&
&                           -gxa(2,jjs+2,ia,iproj,ispin)*gxa(1,jjs  ,ia,iproj,ispinp)
                           rank2c(1,6)=gxa(1,jjs+1,ia,iproj,ispin)*gxa(1,jjs  ,ia,iproj,ispinp)&
&                           +gxa(2,jjs+1,ia,iproj,ispin)*gxa(2,jjs  ,ia,iproj,ispinp)
                           rank2c(2,6)=gxa(1,jjs+1,ia,iproj,ispin)*gxa(2,jjs  ,ia,iproj,ispinp)&
&                           -gxa(2,jjs+1,ia,iproj,ispin)*gxa(1,jjs  ,ia,iproj,ispinp)
                         else if (ilang==3) then
                           call cont22cso(gxa(:,jjs:jjs+5,ia,iproj,ispin),&
&                           gxa(:,jjs:jjs+5,ia,iproj,ispinp),&
&                           gmet,rank2c)
                         else if (ilang==4) then
                           call cont33cso(gxa(:,jjs:jjs+9,ia,iproj,ispin),&
&                           gxa(:,jjs:jjs+9,ia,iproj,ispinp),&
&                           gmet,rank2c)
                         end if
                         if (ilang>1) then
                           do jj=1,3
                             do ii=4,6
                               strsso(ii,jj)=strsso(ii,jj)-2.d0*wt(ilang,iproj)*&
&                               (pauli(1,ispin,ispinp,jj)*rank2c(2,ii)+&
&                               pauli(2,ispin,ispinp,jj)*rank2c(1,ii))
                             end do
                           end do
                         end if
                       end do
                     end if
!                    Enf if (choice==3)
                   end if

!                  ----  Accumulate dynamical matrix contributions if requested.
                   if (choice==4) then
                     ia5=ia+ia3-1
                     do ii=1,ilang2
                       jj=ii-1+jjs
                       do mu=1,6
!                        (includes also factor of 2 from "2*Re[stuff]")
                         enlout(mu+6*(ia5-1))=enlout(mu+6*(ia5-1))+two*&
&                         (gxafac(1,jj,ia,iproj)*dgxdt(1,mu+3,jj,ia,iproj,ispin)&
&                         +gxafac(2,jj,ia,iproj)*dgxdt(2,mu+3,jj,ia,iproj,ispin))
                       end do
                       do mu=1,3
                         enlout(mu+6*(ia5-1))=enlout(mu+6*(ia5-1))+two*&
&                         (dgxdtfac(1,mu,jj,ia,iproj)*dgxdt(1,mu,jj,ia,iproj,ispin)&
&                         +dgxdtfac(2,mu,jj,ia,iproj)*dgxdt(2,mu,jj,ia,iproj,ispin))
                       end do
                       enlout(4+6*(ia5-1))=enlout(4+6*(ia5-1)) +two*&
&                       (dgxdtfac(1,2,jj,ia,iproj)*dgxdt(1,3,jj,ia,iproj,ispin)&
&                       +dgxdtfac(2,2,jj,ia,iproj)*dgxdt(2,3,jj,ia,iproj,ispin))
                       enlout(5+6*(ia5-1))=enlout(5+6*(ia5-1)) +two*&
&                       (dgxdtfac(1,1,jj,ia,iproj)*dgxdt(1,3,jj,ia,iproj,ispin)&
&                       +dgxdtfac(2,1,jj,ia,iproj)*dgxdt(2,3,jj,ia,iproj,ispin))
                       enlout(6+6*(ia5-1))=enlout(6+6*(ia5-1)) +two*&
&                       (dgxdtfac(1,1,jj,ia,iproj)*dgxdt(1,2,jj,ia,iproj,ispin)&
&                       +dgxdtfac(2,1,jj,ia,iproj)*dgxdt(2,2,jj,ia,iproj,ispin))
                     end do
                   end if

!                  ----  Accumulate elastic tensor contributions if requested.

                   if(choice==6) then
!                    XG 081121 : Message to the person who has introduced this CPP option (sorry, I did not have to time to locate who did this ...)
!                    This section of ABINIT should be allowed by default to the user. I have found that on the contrary, the build
!                    system defaults are such that this section is forbidden by default. You might restore this flag if at the same time,
!                    you modify the build system in such a way that by default this section is included, and if the user wants, it can disable it.
!                    #if defined USE_NLSTRAIN_LEGENDRE
                     ia5=ia+ia3-1
                     jjs1=((ilang)*(ilang+1)*(ilang+2))/6
                     jjs2=-3+((ilang+1)*(ilang+2)*(ilang+3))/6
                     jjs3=-9+((ilang+2)*(ilang+3)*(ilang+4))/6
                     jjs4=-19+((ilang+3)*(ilang+4)*(ilang+5))/6

!                    Contribution for two diagonal strains (basically, the nonlocal
!                    enlk)
                     temp(:,1:((rank+1)*(rank+2))/2)= &
&                     gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin)
                     call metcon(rank,gmet,temp,tmpfac)
                     e2nldd=zero
                     do ii=1,((rank+1)*(rank+2))/2
                       e2nldd=e2nldd+&
&                       (gxa(1,jjs-1+ii,ia,iproj,ispin)*tmpfac(1,ii)+&
&                       gxa(2,jjs-1+ii,ia,iproj,ispin)*tmpfac(2,ii))
                     end do

!                    Terms involving one ucvol derivative (diagonal strain only)
!                    and one derivative of nonlocal operator
!                    Loop over strain index
                     do istr2=1,6

!                      rank->rank+2
                       iterm=1
                       temp(:,1:((rank+1)*(rank+2))/2)= &
&                       gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin)
                       call metstr(istr2,rank,iterm,gmet,gprimd,temp,tmpfac)
                       e2nl_tmp(istr2)=0.d0
                       do ii=1,((rank+3)*(rank+4))/2
                         e2nl_tmp(istr2)=e2nl_tmp(istr2)-2.d0*&
&                         (dgxds(1,jjs2-1+ii,ia,iproj,isp)*tmpfac(1,ii)+&
&                         dgxds(2,jjs2-1+ii,ia,iproj,isp)*tmpfac(2,ii))
                       end do
!                      rank->rank

                       iterm=2
                       temp(:,1:((rank+1)*(rank+2))/2)= &
&                       gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin)
                       call metstr(istr2,rank,iterm,gmet,gprimd,temp,tmpfac)
                       do ii=1,((rank+1)*(rank+2))/2
                         e2nl_tmp(istr2)=e2nl_tmp(istr2)-&
&                         (gxa(1,jjs-1+ii,ia,iproj,ispin)*tmpfac(1,ii)+&
&                         gxa(2,jjs-1+ii,ia,iproj,ispin)*tmpfac(2,ii))
                       end do
!                      DEBUG
!                      This and subsequent similar debug sections evaluate the
!                      hermitial conjugate of the contraction immeditely above
!                      and the comparison was useful for development purposes.
!                      rank+2->rank
!                      iterm=3
!                      temp(:,1:((rank+3)*(rank+4))/2)= &
!                      dgxds(:,jjs2:jjs2-1+((rank+3)*(rank+4))/2,ia,iproj,ispin)
!                      call metstr(istr2,rank,iterm,gmet,gprimd,temp,tmpfac)
!                      e2nl_tmp(istr2)=0.d0
!                      do ii=1,((rank+1)*(rank+2))/2
!                      e2nl_tmp(istr2)=e2nl_tmp(istr2)-&
!                      &             (gxa(1,jjs-1+ii,ia,iproj,ispin)*tmpfac(1,ii)+&
!                      &              gxa(2,jjs-1+ii,ia,iproj,ispin)*tmpfac(2,ii))
!                      end do
!                      ENDDEBUG
                     end do

!                    Terms involving two derivatives of the nonlocal operator
!                    Loop over 2nd strain index
                     do istr2=1,6
!                      Loop over 1st strain index, upper triangle only
                       do istr1=1,6
                         iest=istr1+(3*natom+6)*(istr2-1)

!                        Accumulate terms corresponding to two derivatives of ucvol
!                        (simply the nonlocal energy contributin) for both indices
!                        corresponding to diagonal strains

                         if(istr1<=3 .and. istr2<=3) then
                           enlout(iest)= enlout(iest)+wt(ilang,iproj)*e2nldd
                         end if

!                        Accumulate terms computed above from 1st derivatives
!                        when one or more indices corresponds to diagonal strain
                         if(istr2<=3) then
                           enlout(iest)= enlout(iest)+wt(ilang,iproj)*e2nl_tmp(istr1)
                         end if
                         if(istr1<=3) then
                           enlout(iest)= enlout(iest)+wt(ilang,iproj)*e2nl_tmp(istr2)
                         end if

!                        rank->rank+4
                         call contstr21(istr2,istr1,rank,gmet,gprimd,e2nl,&
&                         gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin),&
&                         d2gxds2(:,jjs4:jjs4-1+((rank+5)*(rank+6))/2,ia,iproj,isp))
                         enlout(iest)= enlout(iest)+wt(ilang,iproj)*e2nl

!                        DEBUG
!                        rank+4->rank
!                        call contstr22(istr2,istr1,rank,gmet,gprimd,e2nl,&
!                        &             d2gxds2(:,jjs4:jjs4-1+((rank+5)*(rank+6))/2,ia,iproj,ispin),&
!                        &             gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin))
!                        enlout(iest)= enlout(iest)+wt(ilang,iproj)*e2nl
!                        ENDDEBUG

!                        rank->rank+2
                         call contstr23(istr2,istr1,rank,gmet,gprimd,e2nl,&
&                         gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin),&
&                         dgxds(:,jjs2:jjs2-1+((rank+3)*(rank+4))/2,ia,iproj,isp))
                         enlout(iest)= enlout(iest)+wt(ilang,iproj)*e2nl
!                        DEBUG

!                        rank+2->rank
!                        call contstr24(istr2,istr1,rank,gmet,gprimd,e2nl,&
!                        &             dgxds(:,jjs2:jjs2-1+((rank+3)*(rank+4))/2,ia,iproj,ispin),&
!                        &             gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin))
!                        enlout(iest)= enlout(iest)+wt(ilang,iproj)*e2nl
!                        ENDDEBUG

!                        rank+2->rank+2
                         if(rank<=2) then
                           call contstr25(istr2,istr1,rank,gmet,gprimd,e2nl,&
&                           dgxds(:,jjs2:jjs2-1+((rank+3)*(rank+4))/2,ia,iproj,isp),&
&                           dgxds(:,jjs2:jjs2-1+((rank+3)*(rank+4))/2,ia,iproj,isp))
                         else
                           call contstr25a(istr2,istr1,rank,gmet,gprimd,e2nl,&
&                           dgxds(:,jjs2:jjs2-1+((rank+3)*(rank+4))/2,ia,iproj,isp),&
&                           dgxds(:,jjs2:jjs2-1+((rank+3)*(rank+4))/2,ia,iproj,isp))
                         end if
                         enlout(iest)= enlout(iest)+wt(ilang,iproj)*e2nl

!                        rank->rank
                         call contstr26(istr2,istr1,rank,gmet,gprimd,e2nl,&
&                         gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin),&
&                         gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin))
                         enlout(iest)= enlout(iest)+wt(ilang,iproj)*e2nl

                       end do !istr1

!                      Contributions to internal strain (one cartesian strain and one
!                      reduced-coordinate atomic displacement derivative).
                       iest=7+3*(ia5-1)+(3*natom+6)*(istr2-1)

!                      rank->rank+3
                       call contistr03(istr2,rank,gmet,gprimd,eisnl,&
&                       gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin),&
&                       d2gxdis(:,jjs3:jjs3-1+((rank+4)*(rank+5))/2,ia,iproj,isp))
                       enlout(iest:iest+2)= enlout(iest:iest+2)&
&                       +wt(ilang,iproj)*eisnl(1:3)

!                      DEBUG
!                      rank+3->rank
!                      call contistr30(istr2,rank,gmet,gprimd,eisnl,&
!                      &            d2gxdis(:,jjs3:jjs3-1+((rank+4)*(rank+5))/2,ia,iproj,ispin),&
!                      &            gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin))
!                      enlout(iest:iest+2)= enlout(iest:iest+2)&
!                      &            +wt(ilang,iproj)*eisnl(1:3)
!                      ENDDEBUG

!                      rank+1->rank+2
                       call contistr12(istr2,rank,gmet,gprimd,eisnl,&
&                       dgxdis(:,jjs1:jjs1-1+((rank+2)*(rank+3))/2,ia,iproj,isp),&
&                       dgxds(:,jjs2:jjs2-1+((rank+3)*(rank+4))/2,ia,iproj,isp))
                       enlout(iest:iest+2)= enlout(iest:iest+2)&
&                       +wt(ilang,iproj)*eisnl(1:3)

!                      DEBUG
!                      rank+2->rank+1
!                      call contistr21(istr2,rank,gmet,gprimd,eisnl,&
!                      &            dgxds(:,jjs2:jjs2-1+((rank+3)*(rank+4))/2,ia,iproj,ispin),&
!                      &            dgxdis(:,jjs1:jjs1-1+((rank+2)*(rank+3))/2,ia,iproj,ispin))
!                      enlout(iest:iest+2)= enlout(iest:iest+2)&
!                      &            +wt(ilang,iproj)*eisnl(1:3)
!                      ENDDEBUG

!                      rank->rank+1
                       call contistr01(istr2,rank,gmet,gprimd,eisnl,&
&                       gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin),&
&                       dgxdis(:,jjs1:jjs1-1+((rank+2)*(rank+3))/2,ia,iproj,isp))
                       enlout(iest:iest+2)= enlout(iest:iest+2)&
&                       +wt(ilang,iproj)*eisnl(1:3)
!
!                      DEBUG
!                      rank+1->rank
!                      call contistr10(istr2,rank,gmet,gprimd,eisnl,&
!                      &            dgxdis(:,jjs1:jjs1-1+((rank+2)*(rank+3))/2,ia,iproj,ispin),&
!                      &            gxa(:,jjs:jjs-1+((rank+1)*(rank+2))/2,ia,iproj,ispin))
!                      enlout(iest:iest+2)= enlout(iest:iest+2)&
!                      &            +wt(ilang,iproj)*eisnl(1:3)
!                      ENDDEBUG

                     end do !istr2
                   end if !choice==6

!                  End loop over iproj:
                 end do
!                End condition of existence of a reference state:
               end if

!              End loop over ilang:
             end do

!            End loop over ia:
           end do

!          Operate with the non-local potential on the projected scalars,
!          in order to get matrix element [NOT needed if only force or stress
!          or dynamical matrix is being computed].

           if (signs==2) then
             if(nloalg(2)<=0 .and. choice==2)then
!              Prepare the phase k+q factors for the atoms between ia3 and ia4:
!              they were not prepared previously for nloalg(2)<=0 and choice==2.
               call ph1d3d(ia3,ia4,kgout,matblk,natom,npwout,&
&               n1,n2,n3,phkxredout,ph1d,ph3dout)
             end if

!            call timab(74+choice,1,tsec)
             sign=-1
!            The duplication of the opernl calls has the purpose to avoid
!            copying vectout/vectout_s
             if(ispin==1.or.nspinor==1)then
               if( nloalg(1)==2) then
                 call opernl2(choice,dgxdis,dgxdsfac,d2gxdis,d2gxds2,dgxdtfac,&
&                 ffnlout,gmet,gxafac,ia3,idir,indlmn_s,ispinor,istwf_k,itypat,&
&                 jproj,kgout,kpgout,kptout,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&
&                 mlang5,mlang6,mproj,ndgxdt,dimffnlout,nincat,nkpgout,nlang,nloalg,npwout,&
&                 ntypat,ph3dout,sign,vectout)
               else if( nloalg(1)==3) then
                 call opernl3(choice,dgxdis,dgxdsfac,d2gxdis,d2gxds2,dgxdtfac,&
&                 ffnlout,gmet,gxafac,ia3,idir,indlmn_s,ispinor,istwf_k,itypat,&
&                 jproj,kgout,kpgout,kptout,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&
&                 mlang5,mlang6,mproj,ndgxdt,dimffnlout,nincat,nkpgout,nlang,nloalg,npwout,&
&                 ntypat,ph3dout,sign,vectout)
               else if( nloalg(1)==4) then
                 call opernl4b(choice,dgxdsfac,dgxdtfac,ffnlout,gmet,gxafac,&
&                 ia3,idir,indlmn_s,ispinor,itypat,jproj,kgout,kpgout,kptout,&
&                 lmnmax,matblk,mincat,mlang3,mlang4,mproj,ndgxdt,&
&                 dimffnlout,nincat,nkpgout,nlang,nloalg,npwout,ntypat,ph3dout,vectout)
               end if
             else  ! if ispin == 2
               vectout_s(:,:)=zero
               if( nloalg(1)==2) then
                 call opernl2(choice,dgxdis,dgxdsfac,d2gxdis,d2gxds2,dgxdtfac,&
&                 ffnlout,gmet,gxafac,ia3,idir,indlmn_s,ispinor,istwf_k,itypat,&
&                 jproj,kgout,kpgout,kptout,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&
&                 mlang5,mlang6,mproj,ndgxdt,dimffnlout,nincat,nkpgout,nlang,nloalg,npwout,&
&                 ntypat,ph3dout,sign,vectout_s)
               else if( nloalg(1)==3) then
                 call opernl3(choice,dgxdis,dgxdsfac,d2gxdis,d2gxds2,dgxdtfac,&
&                 ffnlout,gmet,gxafac,ia3,idir,indlmn_s,ispinor,istwf_k,itypat,&
&                 jproj,kgout,kpgout,kptout,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&
&                 mlang5,mlang6,mproj,ndgxdt,dimffnlout,nincat,nkpgout,nlang,nloalg,npwout,&
&                 ntypat,ph3dout,sign,vectout_s)
               else if( nloalg(1)==4) then
                 call opernl4b(choice,dgxds,dgxdtfac,ffnlout,gmet,gxafac,&
&                 ia3,idir,indlmn_s,ispinor,itypat,jproj,kgout,kpgout,kptout,&
&                 lmnmax,matblk,mincat,mlang3,mlang4,mproj,ndgxdt,&
&                 dimffnlout,nincat,nkpgout,nlang,nloalg,npwout,ntypat,ph3dout,vectout_s)
               end if
               vectout(1,1+npwout:2*npwout)=&
&               vectout(1,1+npwout:2*npwout)+vectout_s(1,1:npwout)
               vectout(2,1+npwout:2*npwout)=&
&               vectout(2,1+npwout:2*npwout)+vectout_s(2,1:npwout)

             end if ! end ispin if
           end if ! end signs==2 if

!          End loops over spins:
         end do

!        End condition of existence of a non-local part for that type of atom:
       end if

!      End loop over ispinor:
     end do

!    End sum on atom subset loop, over ia3:
   end do

!  End atom type loop, over itypat:
   ia1=ia2+1
 end do

!De-allocate temporary space.
 ABI_DEALLOCATE(ekb_s)
 ABI_DEALLOCATE(gxa)
 ABI_DEALLOCATE(gxafac)
 ABI_DEALLOCATE(dgxds)
 ABI_DEALLOCATE(dgxdt)
 ABI_DEALLOCATE(dgxdtfac)
 ABI_DEALLOCATE(wt)
 ABI_DEALLOCATE(jproj)
 ABI_DEALLOCATE(temp)
 ABI_DEALLOCATE(tmpfac)
 ABI_DEALLOCATE(dgxdsfac)
 ABI_DEALLOCATE(indlmn_s)
 !if(choice==6)  then
 ABI_DEALLOCATE(dgxdis)
 ABI_DEALLOCATE(d2gxdis)
 ABI_DEALLOCATE(d2gxds2)
 !end if
 !if(nspinor==2) then
 ABI_DEALLOCATE(dgxds_s)
 ABI_DEALLOCATE(dgxdt_s)
 ABI_DEALLOCATE(gxa_s)
 !end if
 !if(nspinor==2.and.choice==6) then
 ABI_DEALLOCATE(dgxdis_s)
 ABI_DEALLOCATE(d2gxdis_s)
 ABI_DEALLOCATE(d2gxds2_s)
 !end if
 if (mpssoang>mpsang)  then
   ABI_DEALLOCATE(pauli)
 end if

!Restore the original content of the vectin array.
!Note that only the first part was modified
 if(istwf_k/=1) then
   call scalewf_nonlop(istwf_k,mpi_enreg,npwin,2,vectin)
 end if

 if (nspinor==2)  then
   ABI_DEALLOCATE(vectin_s)
   ABI_DEALLOCATE(vectout_s)
 end if

 if (mpi_enreg%paral_spinor==1) then
   call xmpi_sum(enlout,mpi_enreg%comm_spinor,ierr)
   call xmpi_sum(strsnl,mpi_enreg%comm_spinor,ierr)
   call xmpi_sum(enlk,mpi_enreg%comm_spinor,ierr)
   call xmpi_sum(strsso,mpi_enreg%comm_spinor,ierr)
 end if

!Save non-local energy
 if((choice==1).and.size(enlout)>0) enlout(1)=enlk ! on test v4/93 size(enlout) can be zero (PMA)

!Do final manipulations to produce strain gradients for
!stress tensor, in cartesian coordinates
 if ((choice==3.or.choice==23) .and. signs==1) then
!  Convert strsnl from reduced to cartesian coordinates
   strsnl_out(:)=0.d0
   call strconv(strsnl,gprimd,strsnl_out)
   strsnl(:) = strsnl_out(:)
!  Add diagonal part (fill up first 6 components of enlout with
!  these gradients; elements 7,8,9 of enlout are not used)
   enlout(1)=strsnl(1)-enlk
   enlout(2)=strsnl(2)-enlk
   enlout(3)=strsnl(3)-enlk
   enlout(4)=strsnl(4)
   enlout(5)=strsnl(5)
   enlout(6)=strsnl(6)
!  Eventually, add spin-orbit part due to gradients of
!  antisymmetric tensor
   if (mpssoang>mpsang) then
     call strsocv(strsso,gprimd,strssoc)
     enlout(1:6)=enlout(1:6)+strssoc(:)
   end if
 end if

!DEBUG
!write(std_out,*)' nonlop_pl: exit '
!ENDDEBUG

contains
!!***

!!****f* ABINIT/trace2
!! NAME
!! trace2
!!
!! FUNCTION
!! Sum indices to compute trace of rank 2 tensor gxa related to l=2
!! $trace=sum_{i,j} {gxa(i,j) gmet(i,j)}$
!!
!! INPUTS
!!  gxa(2,6) = $sum_{G} {e^(2 \pi i G cdot t} {{f_2(|k+G|)} \over {|k+G|^2}} (k+G) cdot (k+G) C(G_{nk})}$
!!  gmet(3,3)=(symmetric) metric tensor in reciprocal space (bohr^-2)
!!
!! OUTPUT
!!  trace(2)=sum_{i,j} {gxa(i,j) gmet(i,j)}$ (Re and Im).
!!
!! NOTES
!! Here index 6 refers to vector components
!! of (k+G) but note tensor is symmetric=>only 6 components.
!! The components are given in the order 11 22 33 32 31 21.
!! The initial 2 handles the Re and Im parts.
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!
!! SOURCE

subroutine trace2(gxa,gmet,trace)

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: gmet(3,3),gxa(2,6)
 real(dp),intent(out) :: trace(2)

!Local variables-------------------------------
!scalars
 integer :: reim

! *************************************************************************

!Write out index summation, Re and Im parts
 do reim=1,2
   trace(reim)=gxa(reim,1)*gmet(1,1)+gxa(reim,2)*gmet(2,2)+&
&   gxa(reim,3)*gmet(3,3)+&
&   2.0d0*(gxa(reim,4)*gmet(3,2)+gxa(reim,5)*gmet(3,1)+&
&   gxa(reim,6)*gmet(2,1))
 end do

end subroutine trace2
!!***

!!****f* ABINIT/strsocv
!! NAME
!! strsocv
!!
!! FUNCTION
!! Convert from antisymmetric storage mode 3x3x3 rank3 tensor in reduced
!! coordinates "red" to symmetric storage mode 3x3 rank2 tensor in
!! cartesian coordinates "cart", using metric tensor "gprimd".
!!
!! INPUTS
!!  red(6,3)=3x3x3 tensor in antisymmetric storage mode,
!!           reduced coordinates
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!
!! OUTPUT
!!  cart(6)=3x3 tensor in symmetric storage mode,
!!          cartesian coordinates
!!
!! NOTES
!! This routine is used to compute spin-orbit stress tensor.
!!
!! red is antisymmetric in first two indices and stored
!!     as 11 22 33 32 31 21.
!! cart is symmetric and stored as 11 22 33 32 31 21.
!!
!!{{\ \begin{eqnarray}
!! cart(1,1) & = &        & red(i,j,2) G(3,i) G(1,j) + red(i,j,3) G(1,i) G(2,j) \nonumber
!! cart(2,2) & = &        & red(i,j,1) G(2,i) G(3,j) + red(i,j,3) G(1,i) G(2,j) \nonumber
!! cart(3,3) & = &        & red(i,j,1) G(2,i) G(3,j) + red(i,j,2) G(3,i) G(1,j) \nonumber
!! cart(3,2) & = &  0.5 ( & red(i,j,3) G(1,i) G(3,j) + red(i,j,2) G(2,i) G(1,j)) \nonumber
!! cart(3,1) & = &  0.5 ( & red(i,j,3) G(3,i) G(2,j) + red(i,j,1) G(2,i) G(1,j)) \nonumber
!! cart(2,1) & = &  0.5 ( & red(i,j,2) G(3,i) G(2,j) + red(i,j,1) G(1,i) G(3,j))
!! \end{eqnarray} }}
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!
!! SOURCE

subroutine strsocv(red,gprimd,cart)

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: gprimd(3,3),red(6,3)
 real(dp),intent(out) :: cart(6)

!Local variables-------------------------------
!scalars
 integer :: ii,jj
!arrays
 real(dp) :: work(3,3,3)

! *************************************************************************

 do ii=1,3
   work(1,1,ii)=0.d0
   work(2,2,ii)=0.d0
   work(3,3,ii)=0.d0
   work(3,2,ii)=red(4,ii) ; work(2,3,ii)=-red(4,ii)
   work(3,1,ii)=red(5,ii) ; work(1,3,ii)=-red(5,ii)
   work(2,1,ii)=red(6,ii) ; work(1,2,ii)=-red(6,ii)
 end do

 cart(:)=0.d0
 do jj=1,3
   do ii=1,3
     cart(1)=cart(1)+work(ii,jj,2)*gprimd(3,ii)*gprimd(1,jj)&
&     +work(ii,jj,3)*gprimd(1,ii)*gprimd(2,jj)
     cart(2)=cart(2)+work(ii,jj,1)*gprimd(2,ii)*gprimd(3,jj)&
&     +work(ii,jj,3)*gprimd(1,ii)*gprimd(2,jj)
     cart(3)=cart(3)+work(ii,jj,1)*gprimd(2,ii)*gprimd(3,jj)&
&     +work(ii,jj,2)*gprimd(3,ii)*gprimd(1,jj)
     cart(4)=cart(4)+work(ii,jj,3)*gprimd(1,ii)*gprimd(3,jj)&
&     +work(ii,jj,2)*gprimd(2,ii)*gprimd(1,jj)
     cart(5)=cart(5)+work(ii,jj,3)*gprimd(3,ii)*gprimd(2,jj)&
&     +work(ii,jj,1)*gprimd(2,ii)*gprimd(1,jj)
     cart(6)=cart(6)+work(ii,jj,2)*gprimd(3,ii)*gprimd(2,jj)&
&     +work(ii,jj,1)*gprimd(1,ii)*gprimd(3,jj)
   end do
 end do
 cart(4)=0.5d0*cart(4)
 cart(5)=0.5d0*cart(5)
 cart(6)=0.5d0*cart(6)

end subroutine strsocv
!!***

!!****f* ABINIT/scalewf_nonlop
!! NAME
!! scalewf_nonlop
!!
!! FUNCTION
!! At the start of nonlop (or similar routines), as well as its end,
!! the wavefunctions, when stored with istwfk/=2,
!! need to be scaled (by a factor of 2 or 1/2),
!! except for the G=0 component.
!! Only the first spinor component is to be modified.
!!
!! INPUTS
!!  istwf_k=storage mode of the vector
!!  mpi_enreg=informations about MPI parallelization
!!  npw=number of planewaves
!!  option=1 multiply by 2
!!        =2 multiply by 1/2
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  vect(2,npw)=vector that is rescaled
!!
!! NOTES
!!  XG030513 : MPIWF One should pay attention to the
!!  G=0 component, that will be only one one proc...
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!
!! SOURCE

subroutine scalewf_nonlop(istwf_k,mpi_enreg,npw,option,vect)

!Arguments ------------------------------------
!This type is defined in defs_mpi
!scalars
 integer,intent(in) :: istwf_k,npw,option
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(inout) :: vect(2,npw)

!Local variables-------------------------------
!scalars
 integer :: ipw
 real(dp) :: scale
 character(len=500) :: message

! *************************************************************************

 DBG_ENTER("COLL")

 if(istwf_k/=1)then

   if(option/=1 .and. option/=2)then
     write(message,'(a,a,a,i0)')&
&     'The argument option should be 1 or 2,',ch10,&
&     'however, option=',option
     MSG_BUG(message)
   end if

   scale=two
   if(option==2)scale=half

!  Storage for the Gamma point. The component of the G=0 vector
!  should not be scaled, and no G=0 imaginary part is allowed.
   if(istwf_k==2)then
     if (mpi_enreg%me_g0==1) then
       vect(2,1)=zero
!$OMP PARALLEL DO
       do ipw=2,npw
         vect(1,ipw)=scale*vect(1,ipw)
         vect(2,ipw)=scale*vect(2,ipw)
       end do
!$OMP END PARALLEL DO
     else
!$OMP PARALLEL DO
       do ipw=1,npw
         vect(1,ipw)=scale*vect(1,ipw)
         vect(2,ipw)=scale*vect(2,ipw)
       end do
!$OMP END PARALLEL DO
     end if
   end if

!  Other storage modes, for k points with time-reversal symmetry.
!  All components should be scaled.
   if(istwf_k>2)then
!$OMP PARALLEL DO
     do ipw=1,npw
       vect(1,ipw)=scale*vect(1,ipw)
       vect(2,ipw)=scale*vect(2,ipw)
     end do
!$OMP END PARALLEL DO
   end if

 end if ! istwf_k/=1

 DBG_EXIT("COLL")

end subroutine scalewf_nonlop
!!***

!!****f* ABINIT/ddkten
!! NAME
!! ddkten
!!
!! FUNCTION
!! Compact or decompact the tensors related to the ffnl(:,1,...)
!! part of the ddk operator, taking into account the direction
!! of the ddk perturbation.
!!
!! INPUTS
!!  compact= if 1, compact from tmpfac
!!  idir=direction of the ddk perturbation
!!  rank=0,1,2, or 3 = rank of tmpfac tensor, also angular momentum (=l)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  temp(2,(rank*(rank+1))/2)=compacted tensor
!!    for l=1, just a scalar
!!    for l=2, a vector
!!  tmpfac(2,(rank+1)*(rank+2)/2)=decompacted tensor
!!    for l=1, a vector
!!    for l=2, a symmetric matrix, stored as
!!     (1 . .)
!!     (6 2 .)
!!     (5 4 3)
!!
!! NOTES
!! For l=0, there is no contribution.
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!
!! SOURCE

subroutine ddkten(compact,idir,rank,temp,tmpfac)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: compact,idir,rank
!arrays
 real(dp),intent(inout) :: temp(2,(rank*(rank+1))/2)
 real(dp),intent(inout) :: tmpfac(2,((rank+1)*(rank+2))/2)

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *************************************************************************

 if(rank/=1 .and. rank/=2 .and. rank/=3)then
   write(message, '(a,i10,a,a,a)' )&
&   'Input rank=',rank,' not allowed.',ch10,&
&   'Possible values are 1,2,3 only.'
   MSG_BUG(message)
 end if

!Take care of p angular momentum
 if(rank==1)then

!  Compaction tmpfac -> temp
   if(compact==1)then
     temp(:,1)=tmpfac(:,idir)

!    Decompaction temp -> tmpfac
   else
     tmpfac(:,1:3)=0.0d0
     tmpfac(:,idir)=temp(:,1)
   end if

!  Take care of d angular momentum
!  rank=2 11->1 22->2 33->3 32->4 31->5 21->6

 else if(rank==2)then

!  Compaction tmpfac -> temp
   if(compact==1)then
     if(idir==1)then
!      Count the number of non-zero derivatives with respect to k(idir)
!      The factor of 2 on the diagonal comes from the derivative with
!      respect to the first K then to the second K
       temp(:,1)=2.0d0*tmpfac(:,1); temp(:,2)=tmpfac(:,6); temp(:,3)=tmpfac(:,5)
     else if(idir==2)then
       temp(:,2)=2.0d0*tmpfac(:,2); temp(:,1)=tmpfac(:,6); temp(:,3)=tmpfac(:,4)
     else if(idir==3)then
       temp(:,3)=2.0d0*tmpfac(:,3); temp(:,1)=tmpfac(:,5); temp(:,2)=tmpfac(:,4)
     end if
!    Decompaction temp -> tmpfac
   else
     tmpfac(:,1:6)=0.0d0
     tmpfac(:,idir)=2.0d0*temp(:,idir)
     if(idir==1)then
       tmpfac(:,5)=temp(:,3); tmpfac(:,6)=temp(:,2)
     else if(idir==2)then
       tmpfac(:,4)=temp(:,3); tmpfac(:,6)=temp(:,1)
     else if(idir==3)then
       tmpfac(:,4)=temp(:,2); tmpfac(:,5)=temp(:,1)
     end if
   end if

!  Take care of f angular momentum
 else if(rank==3)then
!  rank=3 111->1 221->2 331->3 321->4 311->5 211->6 222->7 332->8 322->9 333->10
!  rank=2 11->1 22->2 33->3 32->4 31->5 21->6

!  Compaction tmpfac -> temp
   if(compact==1)then
     if(idir==1)then
!      Count the number of non-zero derivatives with respect to k(idir)
       temp(:,1)=3.0d0*tmpfac(:,1)
       temp(:,2:4)=tmpfac(:,2:4)
       temp(:,5:6)=2.0d0*tmpfac(:,5:6)
     else if(idir==2)then
       temp(:,6)=2.0d0*tmpfac(:,2)
       temp(:,4)=2.0d0*tmpfac(:,9)
       temp(:,5)=tmpfac(:,4)
       temp(:,1)=tmpfac(:,6)
       temp(:,3)=tmpfac(:,8)
       temp(:,2)=3.0d0*tmpfac(:,7)
     else if(idir==3)then
       temp(:,3)=3.0d0*tmpfac(:,10)
       temp(:,5)=2.0d0*tmpfac(:,3)
       temp(:,4)=2.0d0*tmpfac(:,8)
       temp(:,6)=tmpfac(:,4)
       temp(:,1)=tmpfac(:,5)
       temp(:,2)=tmpfac(:,9)
     end if
!    Decompaction temp -> tmpfac
   else
     tmpfac(:,1:10)=0.0d0
     if(idir==1)then
       tmpfac(:,1)=3.0d0*temp(:,1)
       tmpfac(:,2:4)=temp(:,2:4)
       tmpfac(:,5:6)=2.0d0*temp(:,5:6)
     else if(idir==2)then
       tmpfac(:,2)=2.0d0*temp(:,6)
       tmpfac(:,9)=2.0d0*temp(:,4)
       tmpfac(:,4)=temp(:,5)
       tmpfac(:,6)=temp(:,1)
       tmpfac(:,8)=temp(:,3)
       tmpfac(:,7)=3.0d0*temp(:,2)
     else if(idir==3)then
       tmpfac(:,10)=3.0d0*temp(:,3)
       tmpfac(:,3)=2.0d0*temp(:,5)
       tmpfac(:,8)=2.0d0*temp(:,4)
       tmpfac(:,4)=temp(:,6)
       tmpfac(:,5)=temp(:,1)
       tmpfac(:,9)=temp(:,2)
     end if
   end if

 end if

end subroutine ddkten
!!***

end subroutine nonlop_pl
!!***

end module m_nonlop_pl
!!***
