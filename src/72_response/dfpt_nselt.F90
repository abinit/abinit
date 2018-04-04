!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_nselt
!! NAME
!! dfpt_nselt
!!
!! FUNCTION
!! This routine compute the non-stationary expression for the
!! second derivative of the total energy, wrt strain for a whole row of
!! mixed strain derivatives.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DRH, XG, DCA, GMR, MM, AR, MV, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF wavefunctions at k,q.
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!  effmass_free=effective mass for electrons (1. in common case)
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  idir=direction of the perturbation
!!  ipert=type of the perturbation
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the
!!     storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points in the reduced BZ
!!  kxc(nfft,nkxc)=exchange and correlation kernel
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem =number of k points treated by this node.
!!  mk1mem =number of k points treated by this node  (RF data).
!!  mpert =maximum number of ipert
!!  mpi_enreg=information about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  maximum dimension for q points in grids for nonlocal form factors
!!  natom=number of atoms in cell.
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point for each spin
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt_rbz=number of k points in the reduced BZ for this perturbation
!!  nkxc=second dimension of the kxc array. If /=0,
!!   the exchange-correlation kernel must be computed.
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym1=number of symmetry elements in space group consistent with
!!    perturbation
!!  ntypat=number of types of atoms in unit cell.
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band
!!   and k in the reduced Brillouin zone (usually =2)
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  prtbbb=if 1, band-by-band decomposition (also dim of d2bbb)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  qphon(3)=reduced coordinates for the phonon wavelength
!!  rhog(2,nfft)=array for Fourier transform of GS electron density
!!  rhor(nfft,nspden)=GS electron density in electrons/bohr**3.
!!  rhor1(cplex*nfft,nspden)=RF electron density in electrons/bohr**3.
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  typat(natom)=type integer for each atom in cell
!!  ucvol=unit cell volume in bohr**3.
!!  wtk_rbz(nkpt_rbz)=weight assigned to each k point in the reduced BZ
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang)= real spherical harmonics for each G and k+q point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical for each G and k point
!!  ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*useylm)= gradients of real spherical for each G and k+g point
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert)=flags for each element of the 2DTE (=1 if computed)
!!  d2bbb(2,3,3,mpert,mband,mband*prtbbb)=band by band decomposition of some
!!       second order derivatives
!!  d2lo(2,3,mpert,3,mpert)=local contributions to the 2DTEs
!!  d2nl(2,3,mpert,3,mpert)=non-local contributions to the 2DTEs
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      destroy_hamiltonian,dfpt_mkcore,dfpt_mkvxcstr,dfpt_nsteltwf,dotprod_vn
!!      hartrestr,init_hamiltonian,stresssym,vlocalstr,wrtout,xmpi_barrier
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_nselt(blkflg,cg,cg1,cplex,&
& d2bbb,d2lo,d2nl,ecut,ecutsm,effmass_free,&
& gmet,gprimd,gsqcut,idir,&
& ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,mband,mgfft,&
& mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,mpw1,&
& natom,nband_rbz,nfft,ngfft,&
& nkpt_rbz,nkxc,nloalg,npwarr,npwar1,nspden,nspinor,nsppol,&
& nsym1,ntypat,occ_rbz,&
& paral_kgb,ph1d,prtbbb,psps,qphon,rhog,&
& rhor,rhor1,rmet,rprimd,symrc1,typat,ucvol,&
& wtk_rbz,&
& xred,ylm,ylm1,ylmgr,ylmgr1)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_xmpi

 use m_cgtools,    only : dotprod_vn
 use m_pawtab,     only : pawtab_type
 use m_hamiltonian,only : init_hamiltonian,destroy_hamiltonian,gs_hamiltonian_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_nselt'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_72_response, except_this_one => dfpt_nselt
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mband,mgfft,mk1mem
 integer,intent(in) :: mkmem,mpert,mpsang,mpw,mpw1,natom,nfft,nkpt_rbz
 integer,intent(in) :: nkxc,nspden,nspinor,nsppol,nsym1,ntypat
 integer,intent(in) :: paral_kgb,prtbbb
 real(dp),intent(in) :: ecut,ecutsm,effmass_free,gsqcut,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: istwfk_rbz(nkpt_rbz)
 integer,intent(in) :: kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem)
 integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol),ngfft(18)
 integer,intent(in) :: nloalg(3),npwar1(nkpt_rbz),npwarr(nkpt_rbz)
 integer,intent(in) :: symrc1(3,3,nsym1),typat(natom)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: gprimd(3,3),kpt_rbz(3,nkpt_rbz),kxc(nfft,nkxc)
 real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qphon(3),rhog(2,nfft)
 real(dp),intent(in) :: rhor(nfft,nspden)
 real(dp),intent(in) :: rhor1(cplex*nfft,nspden),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz),xred(3,natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*psps%useylm)
 real(dp),intent(out) :: d2bbb(2,3,3,mpert,mband,mband*prtbbb)
 real(dp),intent(inout) :: d2lo(2,3,mpert,3,mpert)
 real(dp),intent(inout) :: d2nl(2,3,mpert,3,mpert)

!Local variables-------------------------------
!scalars
 integer :: ban2tot,bantot,bd2tot_index,bdtot_index
 integer :: icg,icg1,idir1,ifft,ii,ikg,ikg1,ikpt,comm
 integer :: ilm,ipert1,ispden,isppol,istr1,istwf_k
 integer :: mbd2kpsp,mbdkpsp,me,n1,n2,n3,n3xccc,n4,n5,n6
 integer :: nband_k,nfftot,npw1_k,npw_k,option
 real(dp) :: doti,dotr
 real(dp) :: wtk_k
 character(len=500) :: message
 type(gs_hamiltonian_type) :: gs_hamk
!arrays
 integer :: ikpt_fbz(3)
 integer,allocatable :: kg1_k(:,:),kg_k(:,:)
 real(dp) :: kpoint(3),restr(6),dummy(0,0)
 real(dp),allocatable :: d2bbb_k(:,:,:,:),d2nl_k(:,:,:)
 real(dp),allocatable :: occ_k(:)
 real(dp),allocatable :: vhartr01(:),vpsp1(:),vxc1(:,:),xccc3d1(:),ylm1_k(:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr1_k(:,:,:),ylmgr_k(:,:,:)
 type(pawtab_type) :: pawtab_dum(0)

! *********************************************************************

!Init me
 comm = mpi_enreg%comm_cell
 me   = mpi_enreg%me_kpt

!Unit numbers

!Zero only portion of nonlocal matrix to be computed here
 d2nl(:,:,natom+3:natom+4,idir,ipert)=zero
 bdtot_index=0
 bd2tot_index=0
 icg=0
 icg1=0
 mbdkpsp=mband*nkpt_rbz*nsppol
 mbd2kpsp=2*mband**2*nkpt_rbz*nsppol

!Update list of computed matrix elements
 if((ipert==natom+3) .or. (ipert==natom+4)) then
!  Eventually expand when strain coupling to other perturbations is
!  implemented
   do ipert1=natom+3,natom+4
     do idir1=1,3
       blkflg(idir1,ipert1,idir,ipert)=1
     end do
   end do
 end if

!allocate(enl1nk(mbdkpsp))
 ABI_ALLOCATE(d2bbb_k,(2,3,mband,mband*prtbbb))
 ABI_ALLOCATE(d2nl_k,(2,3,mpert))


 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(kg1_k,(3,mpw1))


 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 n4=ngfft(4) ; n5=ngfft(5) ; n6=ngfft(6)
 nfftot=n1*n2*n3

!Initialize Hamiltonian (k-independent terms) - NCPP only
 call init_hamiltonian(gs_hamk,psps,pawtab_dum,nspinor,nsppol,nspden,natom,&
& typat,xred,nfft,mgfft,ngfft,rprimd,nloalg,ph1d=ph1d)

 bantot = 0
 ban2tot = 0

!LOOP OVER SPINS
 do isppol=1,nsppol

   if (nsppol/=1) then
     write(message,*)' ****  In dfpt_nselt for isppol=',isppol
     call wrtout(std_out,message,'COLL')
   end if

   ikg=0
   ikg1=0

   ikpt_fbz(1:3)=0

!  BIG FAT k POINT LOOP
   do ikpt=1,nkpt_rbz

     nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
     istwf_k=istwfk_rbz(ikpt)
     npw_k=npwarr(ikpt)
     npw1_k=npwar1(ikpt)
     kpoint(:)=kpt_rbz(:,ikpt)

     bantot = bantot + nband_k
     ban2tot = ban2tot + 2*nband_k**2

     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       bdtot_index=bdtot_index+nband_k
       bd2tot_index=bd2tot_index+2*nband_k**2
!      Skip the rest of the k-point loop
       cycle
     end if

     ABI_ALLOCATE(occ_k,(nband_k))

     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang))
     ABI_ALLOCATE(ylm1_k,(npw1_k,mpsang*mpsang))
     if (ipert==natom+3.or.ipert==natom+4) then
       ABI_ALLOCATE(ylmgr_k,(npw_k,3,mpsang*mpsang))
       ABI_ALLOCATE(ylmgr1_k,(npw1_k,3,mpsang*mpsang))
     end if

!    enl1_k(:)=zero
     d2nl_k(:,:,:)=zero
     if(prtbbb==1)d2bbb_k(:,:,:,:)=zero
     occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)

     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     if (psps%useylm==1) then
       do ilm=1,mpsang*mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
       if (ipert==natom+3.or.ipert==natom+4) then
         do ilm=1,mpsang*mpsang
           do ii=1,3
             ylmgr_k(1:npw_k,ii,ilm)=ylmgr(1+ikg:npw_k+ikg,ii,ilm)
           end do
         end do
       end if
     end if

     wtk_k=wtk_rbz(ikpt)

     kg1_k(:,:) = 0

     kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
     if (psps%useylm==1) then
       do ilm=1,mpsang*mpsang
         ylm1_k(1:npw1_k,ilm)=ylm1(1+ikg1:npw1_k+ikg1,ilm)
       end do
       if (ipert==natom+3.or.ipert==natom+4) then
         do ilm=1,mpsang*mpsang
           do ii=1,3
             ylmgr1_k(1:npw1_k,ii,ilm)=ylmgr1(1+ikg1:npw1_k+ikg1,ii,ilm)
           end do
         end do
       end if
     end if

!    Compute the eigenvalues, wavefunction,
!    contributions to kinetic energy, nonlocal energy, forces,
!    and update of rhor1 to this k-point and this spin polarization.

!    Note that dfpt_nsteltwf is called with kpoint, while kpt is used inside dfpt_vtowfk
     call dfpt_nsteltwf(cg,cg1,d2nl_k,ecut,ecutsm,effmass_free,gs_hamk,icg,icg1,ikpt,isppol,&
&     istwf_k,kg_k,kg1_k,kpoint,mband,mkmem,mk1mem,mpert,mpi_enreg,mpw,mpw1,natom,nband_k,&
&     npw_k,npw1_k,nspinor,nsppol,ntypat,occ_k,psps,rmet,wtk_k,ylm_k,ylmgr_k)
     d2nl(:,:,:,idir,ipert)=d2nl(:,:,:,idir,ipert)+d2nl_k(:,:,:)
     if(prtbbb==1)then
       d2bbb(:,:,idir,ipert,:,:) = d2bbb(:,:,idir,ipert,:,:) + &
&       d2bbb_k(:,:,:,:)
     end if

     ABI_DEALLOCATE(occ_k)

!    Keep track of total number of bands (all k points so far, even for
!    k points not treated by me)
     bdtot_index=bdtot_index+nband_k
     bd2tot_index=bd2tot_index+2*nband_k**2

!    Shift array memory
     if (mkmem/=0) then
       icg=icg+npw_k*nspinor*nband_k
       ikg=ikg+npw_k
     end if
     if (mk1mem/=0) then
       icg1=icg1+npw1_k*nspinor*nband_k
       ikg1=ikg1+npw1_k
     end if
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylm1_k)
     if (ipert==natom+3.or.ipert==natom+4)  then
       ABI_DEALLOCATE(ylmgr_k)
       ABI_DEALLOCATE(ylmgr1_k)
     end if

   end do ! End big k point loop
 end do ! End loop over spins

 if(xmpi_paral==1)then
   call xmpi_barrier(comm)
   call wrtout(std_out,' dfpt_nselt: loop on k-points and spins done in parallel','COLL')
 end if

!Treat now varying occupation numbers
!if(occopt>=3 .and. occopt <=8) then
!SUPPRESSED metallic coding of vtorho

!Treat fixed occupation numbers
!else

!Accumulation over parallel processed now carried out for all terms
!in dfpt_nstdy.f

!End of test on varying or fixed occupation numbers
!end if

!The imaginary part of d2nl will be must be set to zero here since
!time-reversal symmetry will always be true for the strain peturbation.
!The symmetry-reduced kpt set will leave a non-zero imaginary part.

 d2nl(2,:,natom+3:natom+4,idir,ipert)=zero

!Symmetrize the non-local contributions,
!as was needed for the stresses in a ground-state calculation

 if (nsym1>1) then
!  Pack like symmetric-storage cartesian stress tensor
   ii=0
   do ipert1=natom+3,natom+4
     do idir1=1,3
       ii=ii+1
       restr(ii)=d2nl(1,idir1,ipert1,idir,ipert)
     end do
   end do
!  Do the symmetrization using the ground state routine
   call stresssym(gprimd,nsym1,restr,symrc1)
!  Unpack symmetrized stress tensor
   ii=0
   do ipert1=natom+3,natom+4
     do idir1=1,3
       ii=ii+1
       d2nl(1,idir1,ipert1,idir,ipert)=restr(ii)
     end do
   end do
 end if !nsym>1

!----------------------------------------------------------------------------
!Now, treat the local contribution

 ABI_ALLOCATE(vpsp1,(cplex*nfft))
 n3xccc=0
 if(psps%n1xccc/=0)n3xccc=nfft
 ABI_ALLOCATE(xccc3d1,(cplex*n3xccc))
 ABI_ALLOCATE(vxc1,(cplex*nfft,nspden))
 ABI_ALLOCATE(vhartr01,(nfft))
 xccc3d1(:)=zero

!Double loop over strain perturbations
 do ipert1=natom+3,natom+4
   do idir1=1,3
     if(ipert1==natom+3) then
       istr1=idir1
     else
       istr1=idir1+3
     end if

!    Get first-order local potential.
     call vlocalstr(gmet,gprimd,gsqcut,istr1,mgfft,mpi_enreg,&
&     psps%mqgrid_vl,natom,gs_hamk%nattyp,nfft,ngfft,ntypat,paral_kgb,ph1d,psps%qgrid_vl,&
&     ucvol,psps%vlspl,vpsp1)

!    Get first-order hartree potential.
     call hartrestr(gsqcut,idir1,ipert1,mpi_enreg,natom,nfft,ngfft,&
&     paral_kgb,rhog,rprimd,vhartr01)

!    Get first-order exchange-correlation potential
     if(psps%n1xccc/=0)then
       call dfpt_mkcore(cplex,idir1,ipert1,natom,ntypat,n1,psps%n1xccc,&
&       n2,n3,qphon,rprimd,typat,ucvol,psps%xcccrc,psps%xccc1d,xccc3d1,xred)
     end if ! psps%n1xccc/=0

     option=0
     call dfpt_mkvxcstr(cplex,idir1,ipert1,kxc,mpi_enreg,natom,nfft,ngfft,&
&     dummy,dummy,nkxc,nspden,n3xccc,option,paral_kgb,qphon,rhor,rhor1,rprimd,&
&     0,0,vxc1,xccc3d1)

!    Combines density j2 with local potential j1
     do ispden=1,min(nspden,2)
       do ifft=1,cplex*nfft
         vxc1(ifft,ispden)=vxc1(ifft,ispden)+vpsp1(ifft)+vhartr01(ifft)
       end do
     end do
     call dotprod_vn(cplex,rhor1,dotr,doti,nfft,nfftot,nspden,2,vxc1,ucvol)
     write(std_out,*)
     d2lo(1,idir1,ipert1,idir,ipert)=dotr
     d2lo(2,idir1,ipert1,idir,ipert)=doti
   end do ! istr1
 end do ! ipert1

 call destroy_hamiltonian(gs_hamk)

 ABI_DEALLOCATE(vxc1)
 ABI_DEALLOCATE(xccc3d1)
 ABI_DEALLOCATE(vhartr01)

 ABI_DEALLOCATE(d2bbb_k)
 ABI_DEALLOCATE(d2nl_k)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(kg1_k)
 ABI_DEALLOCATE(vpsp1)

end subroutine dfpt_nselt
!!***
