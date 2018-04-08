!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_nstwf
!! NAME
!! dfpt_nstwf
!!
!! FUNCTION
!! This routine computes the non-local contribution to the
!! 2DTE matrix elements, in the non-stationary formulation
!! Only for norm-conserving pseudopotentials (no PAW)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (XG,AR,MB,MVer,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions at k
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF wavefunctions at k,q.
!!  ddkfil(3)=unit numbers for the three possible ddk files for ipert1
!!       equal to 0 if no dot file is available for this direction
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eig_k(mband*nsppol)=GS eigenvalues at k (hartree)
!!  eig1_k(2*nsppol*mband**2)=matrix of first-order eigenvalues (hartree)
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  icg=shift to be applied on the location of data in the array cg
!!  icg1=shift to be applied on the location of data in the array cg1
!!  idir=direction of the current perturbation
!!  ikpt=number of the k-point
!!  ipert=type of the perturbation
!!  isppol=1 for unpolarized, 2 for spin-polarized
!!  istwf_k=parameter that describes the storage of wfs
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kg1_k(3,npw1_k)=reduced planewave coordinates at k+q, with RF k points
!!  kpt(3)=reduced coordinates of k point
!!  kpq(3)=reduced coordinates of k+q point
!!  mkmem =number of k points treated by this node
!!  mk1mem =number of k points treated by this node (RF data)
!!  mpert =maximum number of ipert
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  nband_k=number of bands at this k point for that spin polarization
!!  npw_k=number of plane waves at this k point
!!  npw1_k=number of plane waves at this k+q point
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rmet(3,3)=real space metric (bohr**2)
!!  ddks(3)<wfk_t>=struct info for for the three possible DDK files for ipert1
!!  wtk_k=weight assigned to the k point.
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k+q point
!!
!! OUTPUT
!!  d2bbb_k(2,3,mband,mband*prtbbb)=band by band decomposition of the second
!!   order derivatives, for the present k point, and perturbation idir, ipert
!!  d2nl_k(2,3,mpert)=non-local contributions to
!!   non-stationary 2DTE, for the present k point, and perturbation idir, ipert
!!
!! TODO
!!  XG 20141103 The localization tensor cannot be defined in the metallic case. It should not be computed.
!!
!! PARENTS
!!      dfpt_nstdy
!!
!! CHILDREN
!!      destroy_rf_hamiltonian,dotprod_g,gaugetransfo,getgh1c
!!      init_rf_hamiltonian,load_k_hamiltonian,load_k_rf_hamiltonian
!!      load_kprime_hamiltonian,mkffnl,mkkpg,timab,wfk_read_bks
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dfpt_nstwf(cg,cg1,ddkfil,dtset,d2bbb_k,d2nl_k,eig_k,eig1_k,gs_hamkq,&
&                 icg,icg1,idir,ikpt,ipert,isppol,istwf_k,kg_k,kg1_k,kpt,kpq,&
&                 mkmem,mk1mem,mpert,mpi_enreg,mpw,mpw1,nband_k,npw_k,npw1_k,nsppol,&
&                 occ_k,psps,rmet,ddks,wtk_k,ylm,ylm1)


 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_profiling_abi
 use m_cgtools
 use m_hamiltonian
 use m_errors
 use m_wfk
 use m_xmpi

 use m_time,    only : timab
 use m_pawcprj, only : pawcprj_type
 use m_kg,      only : mkkpg

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_nstwf'
 use interfaces_66_nonlocal
 use interfaces_66_wfs
 use interfaces_72_response, except_this_one => dfpt_nstwf
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,icg1,idir,ikpt,ipert,isppol,istwf_k
 integer,intent(in) :: mkmem,mk1mem,mpert,mpw,mpw1,nsppol
 integer,intent(inout) :: nband_k,npw1_k,npw_k
 real(dp),intent(in) :: wtk_k
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: ddkfil(3),kg1_k(3,npw1_k)
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*nsppol)
 real(dp),intent(in) :: eig_k(dtset%mband*nsppol),kpt(3),kpq(3),occ_k(nband_k),rmet(3,3)
 real(dp),intent(in) :: ylm(npw_k,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(npw1_k,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout) :: eig1_k(2*nsppol*dtset%mband**2)
 real(dp),intent(out) :: d2bbb_k(2,3,dtset%mband,dtset%mband*dtset%prtbbb)
 real(dp),intent(inout) :: d2nl_k(2,3,mpert)
 type(wfk_t),intent(inout) :: ddks(3)

!Local variables-------------------------------
!scalars
 integer :: berryopt,dimffnl,dimffnl1,dimph3d
 integer :: iband,ider,idir1,ipert1,ipw,jband,nband_kocc,nkpg,nkpg1 !ierr,ii
 integer :: npw_disk,nsp,optlocal,optnl,opt_gvnl1,sij_opt,tim_getgh1c,usevnl
 logical :: ddk
 real(dp) :: aa,dot1i,dot1r,dot2i,dot2r,dot_ndiagi,dot_ndiagr,doti,dotr,lambda
 character(len=500) :: msg
 type(rf_hamiltonian_type) :: rf_hamkq
!arrays
 integer :: ik_ddks(3)
 real(dp) :: dum_grad_berry(1,1),dum_gvnl1(1,1),dum_gs1(1,1),dum_ylmgr(1,3,1),tsec(2)
 real(dp),allocatable :: cg_k(:,:),cwave0(:,:),cwavef(:,:),cwavef_da(:,:)
 real(dp),allocatable :: cwavef_db(:,:),dkinpw(:),eig2_k(:),ffnl1(:,:,:,:),ffnlk(:,:,:,:)
 real(dp),allocatable :: gvnl1(:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:),ph3d(:,:,:)
 type(pawcprj_type),allocatable :: dum_cwaveprj(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Not valid for PAW
 if (psps%usepaw==1) then
   msg='  This routine cannot be used for PAW (use pawnst3 instead) !'
   MSG_BUG(msg)
 end if

!Keep track of total time spent in dfpt_nstwf
 call timab(102,1,tsec)
 tim_getgh1c=2

!Miscelaneous inits
 ABI_DATATYPE_ALLOCATE(dum_cwaveprj,(0,0))
 ddk=(ipert==dtset%natom+1.or.ipert==dtset%natom+10.or.ipert==dtset%natom+11)

!Additional allocations
 if (.not.ddk) then
   ABI_ALLOCATE(dkinpw,(npw_k))
   ABI_ALLOCATE(kinpw1,(npw1_k))
   kinpw1=zero;dkinpw=zero
 else
   ABI_ALLOCATE(dkinpw,(0))
   ABI_ALLOCATE(kinpw1,(0))
 end if
 ABI_ALLOCATE(gvnl1,(2,npw1_k*dtset%nspinor))
 ABI_ALLOCATE(eig2_k,(2*nsppol*dtset%mband**2))
 ABI_ALLOCATE(cwave0,(2,npw_k*dtset%nspinor))
 ABI_ALLOCATE(cwavef,(2,npw1_k*dtset%nspinor))

!Compute (k+G) vectors
 nkpg=0;if (.not.ddk) nkpg=3*gs_hamkq%nloalg(3)
 ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
 if (nkpg>0) then
   call mkkpg(kg_k,kpg_k,kpt,nkpg,npw_k)
 end if

!Compute (k+q+G) vectors
 nkpg1=0;if (.not.ddk) nkpg1=3*gs_hamkq%nloalg(3)
 ABI_ALLOCATE(kpg1_k,(npw1_k,nkpg1))
 if (nkpg1>0) then
   call mkkpg(kg1_k,kpg1_k,kpq,nkpg1,npw1_k)
 end if

!Compute nonlocal form factors ffnl at (k+G)
 dimffnl=0;if (.not.ddk) dimffnl=1
 ABI_ALLOCATE(ffnlk,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
 if (.not.ddk) then
   ider=0
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnlk,psps%ffspl,gs_hamkq%gmet,&
&   gs_hamkq%gprimd,ider,ider,psps%indlmn,kg_k,kpg_k,kpt,psps%lmnmax,&
&   psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,psps%ntypat,psps%pspso,&
&   psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm,dum_ylmgr)
 end if

!Compute nonlocal form factors ffnl1 at (k+q+G)
 dimffnl1=0;if (.not.ddk) dimffnl1=1
 ABI_ALLOCATE(ffnl1,(npw1_k,dimffnl1,psps%lmnmax,psps%ntypat))
 if (.not.ddk) then
   ider=0
   call mkffnl(psps%dimekb,dimffnl1,psps%ekb,ffnl1,psps%ffspl,gs_hamkq%gmet,&
&   gs_hamkq%gprimd,ider,ider,psps%indlmn,kg1_k,kpg1_k,kpq,&
&   psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg1,npw1_k,psps%ntypat,&
&   psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm1,dum_ylmgr)
 end if

!Load k-dependent part in the Hamiltonian datastructure
 call load_k_hamiltonian(gs_hamkq,kpt_k=kpt,npw_k=npw_k,istwf_k=istwf_k,&
& kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnlk)

!Load k+q-dependent part in the Hamiltonian datastructure
!    Note: istwf_k is imposed to 1 for RF calculations (should use istwf_kq instead)
 dimph3d=0;if (.not.ddk) dimph3d=gs_hamkq%matblk
 ABI_ALLOCATE(ph3d,(2,npw1_k,dimph3d))
 call load_kprime_hamiltonian(gs_hamkq,kpt_kp=kpq,npw_kp=npw1_k,istwf_kp=istwf_k,&
& kinpw_kp=kinpw1,kg_kp=kg1_k,kpg_kp=kpg1_k,ffnl_kp=ffnl1,&
& ph3d_kp=ph3d,compute_ph3d=(.not.ddk))

!Load k-dependent part in the 1st-order Hamiltonian datastructure
 call load_k_rf_hamiltonian(rf_hamkq,npw_k=npw_k,dkinpw_k=dkinpw)

!Take care of the npw and kg records
!NOTE : one should be able to modify the rwwf routine to take care
!of the band parallelism, which is not the case yet ...
 ik_ddks = 0
 do idir1=1,3
   if (ddkfil(idir1)/=0)then
!    Read npw record
     nsp=dtset%nspinor
     ik_ddks(idir1) = wfk_findk(ddks(idir1), kpt)
     ABI_CHECK(ik_ddks(idir1) /= -1, "Cannot find kpt")
     npw_disk = ddks(idir1)%hdr%npwarr(ik_ddks(idir1))
     if (npw_k /= npw_disk) then
       write(unit=msg,fmt='(a,i3,a,i5,a,i3,a,a,i5,a,a,i5)')&
&       'For isppol = ',isppol,', ikpt = ',ikpt,' and idir = ',idir,ch10,&
&       'the number of plane waves in the ddk file is equal to', npw_disk,ch10,&
&       'while it should be ',npw_k
       MSG_BUG(msg)
     end if
   end if
 end do

 if (ipert==dtset%natom+1) then
   nband_kocc = 0
   do iband = 1,nband_k
     if (abs(occ_k(iband)) > tol8) nband_kocc = nband_kocc + 1
     nband_kocc = max (nband_kocc, 1)
   end do
 end if

 if(dtset%prtbbb==1)then
   ABI_ALLOCATE(cwavef_da,(2,npw1_k*dtset%nspinor))
   ABI_ALLOCATE(cwavef_db,(2,npw1_k*dtset%nspinor))
   ABI_ALLOCATE(cg_k,(2,npw_k*dtset%nspinor*nband_k))
   if ((ipert == dtset%natom + 1).or.(ipert <= dtset%natom).or. &
&   (ipert == dtset%natom + 2).or.(ipert == dtset%natom + 5)) then
     cg_k(:,:) = cg(:,1+icg:icg+nband_k*npw_k*dtset%nspinor)
   end if
   d2bbb_k(:,:,:,:) = zero
 end if

!Loop over bands
 do iband=1,nband_k

   if(mpi_enreg%proc_distrb(ikpt,iband,isppol) /= mpi_enreg%me_kpt) cycle

!  Read ground-state wavefunctions
   if (dtset%prtbbb==0 .or. ipert==dtset%natom+2) then
     cwave0(:,:)=cg(:,1+(iband-1)*npw_k*dtset%nspinor+icg:iband*npw_k*dtset%nspinor+icg)
   else    ! prtbbb==1 and ipert<=natom , already in cg_k
     cwave0(:,:)=cg_k(:,1+(iband-1)*npw_k*dtset%nspinor:iband*npw_k*dtset%nspinor)
   end if

!  Get first-order wavefunctions
   cwavef(:,:)=cg1(:,1+(iband-1)*npw1_k*dtset%nspinor+icg1:iband*npw1_k*dtset%nspinor+icg1)

!  In case non ddk perturbation
   if (ipert /= dtset%natom + 1) then

     do ipert1=1,mpert

       if( ipert1<=dtset%natom .or. ipert1==dtset%natom+2 )then

!        Initialize data for NL 1st-order hamiltonian
         call init_rf_hamiltonian(1,gs_hamkq,ipert1,rf_hamkq)

         if (((ipert <= dtset%natom).or.(ipert == dtset%natom + 2)) &
&         .and.(ipert1 == dtset%natom+2).and. dtset%prtbbb==1) then
           call gaugetransfo(cg_k,cwavef,cwavef_db,eig_k,eig1_k,iband,nband_k, &
&           dtset%mband,npw_k,npw1_k,dtset%nspinor,nsppol,occ_k)
           cwavef(:,:) = cwavef_db(:,:)
         end if

!        Define the direction along which to move the atom :
!        the polarisation (ipert1,idir1) is refered as j1.
         do idir1=1,3
           if (ipert1<=dtset%natom.or.(ipert1==dtset%natom+2.and.ddkfil(idir1)/=0)) then

!            Get |Vnon-locj^(1)|u0> :
!            First-order non-local, applied to zero-order wavefunction
!            This routine gives MINUS the non-local contribution

!            ==== Atomic displ. perturbation
             if( ipert1<=dtset%natom )then
               lambda=eig_k((isppol-1)*nband_k+iband)
               berryopt=1;optlocal=0;optnl=1;usevnl=0;opt_gvnl1=0;sij_opt=0
               call getgh1c(berryopt,cwave0,dum_cwaveprj,gvnl1,dum_grad_berry,&
&               dum_gs1,gs_hamkq,dum_gvnl1,idir1,ipert1,lambda,mpi_enreg,optlocal,&
&               optnl,opt_gvnl1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)

!              ==== Electric field perturbation
             else if( ipert1==dtset%natom+2 )then
               ! TODO: Several tests fail here ifdef HAVE_MPI_IO_DEFAULT
               ! The problem is somehow related to the use of MPI-IO file views!.
               call wfk_read_bks(ddks(idir1), iband, ik_ddks(idir1), isppol, xmpio_single, cg_bks=gvnl1, &
               eig1_bks=eig2_k(1+(iband-1)*2*nband_k:))
                 !eig1_bks=eig2_k(1+(iband-1)*2*nband_k:2*iband*nband_k))
               !write(777,*)"eig2_k, gvnl1 for band: ",iband,", ikpt: ",ikpt
               !do ii=1,2*nband_k
               !  write(777,*)eig2_k(ii+(iband-1))
               !end do
               !write(777,*)gvnl1

!              In case of band-by-band,
!              construct the first-order wavefunctions in the diagonal gauge
               if (((ipert <= dtset%natom).or.(ipert == dtset%natom + 2)).and.(dtset%prtbbb==1)) then
                 call gaugetransfo(cg_k,gvnl1,cwavef_da,eig_k,eig2_k,iband,nband_k, &
&                 dtset%mband,npw_k,npw1_k,dtset%nspinor,nsppol,occ_k)
                 gvnl1(:,:) = cwavef_da(:,:)
               end if
!              Multiplication by -i
               do ipw=1,npw1_k*dtset%nspinor
                 aa=gvnl1(1,ipw)
                 gvnl1(1,ipw)=gvnl1(2,ipw)
                 gvnl1(2,ipw)=-aa
               end do
             end if

!            MVeithen 021212 :
!            1) Case ipert1 = natom + 2 and ipert = natom + 2:
!            the second derivative of the energy with respect to an electric
!            field is computed from Eq. (38) of X. Gonze, PRB 55 ,10355 (1997).
!            The evaluation of this formula needs the operator $i \frac{d}{dk}.
!            2) Case ipert1 = natom + 2 and ipert < natom:
!            the computation of the Born effective charge tensor uses
!            the operator $-i \frac{d}{dk}.
             if (ipert==dtset%natom+2) gvnl1(:,:) = -gvnl1(:,:)

!            <G|Vnl1|Cnk> is contained in gvnl1
!            construct the matrix element (<uj2|vj1|u0>)complex conjug and add it to the 2nd-order matrix
             call dotprod_g(dotr,doti,istwf_k,npw1_k*dtset%nspinor,2,cwavef,gvnl1,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
             d2nl_k(1,idir1,ipert1)=d2nl_k(1,idir1,ipert1)+wtk_k*occ_k(iband)*two*dotr
             d2nl_k(2,idir1,ipert1)=d2nl_k(2,idir1,ipert1)-wtk_k*occ_k(iband)*two*doti

!            Band by band decomposition of the Born effective charges
!            calculated from a phonon perturbation
             if(dtset%prtbbb==1)then
               d2bbb_k(1,idir1,iband,iband) = wtk_k*occ_k(iband)*two*dotr
               d2bbb_k(2,idir1,iband,iband) = -one*wtk_k*occ_k(iband)*two*doti
             end if

           end if
         end do

         call destroy_rf_hamiltonian(rf_hamkq)
       end if
     end do
   end if     ! ipert /= natom +1

!  Compute the localization tensor

   if (ipert==dtset%natom+1) then

     ipert1=dtset%natom+1
     if(dtset%prtbbb==1)then
       call gaugetransfo(cg_k,cwavef,cwavef_db,eig_k,eig1_k,iband,nband_k, &
&       dtset%mband,npw_k,npw1_k,dtset%nspinor,nsppol,occ_k)
       cwavef(:,:) = cwavef_db(:,:)
     end if

     do idir1 = 1,3
       eig2_k(:) = zero
       gvnl1(:,:) = zero
       if (idir == idir1) then
         gvnl1(:,:) = cwavef(:,:)
         eig2_k(:) = eig1_k(:)
       else
         if (ddkfil(idir1) /= 0) then
           call wfk_read_bks(ddks(idir1), iband, ik_ddks(idir1), isppol, xmpio_single, cg_bks=gvnl1, &
           eig1_bks=eig2_k(1+(iband-1)*2*nband_k:))
             !eig1_bks=eig2_k(1+(iband-1)*2*nband_k:2*iband*nband_k))
           !write(778,*)"eig2_k, gvnl1 for band: ",iband,", ikpt: ",ikpt
           !do ii=1,2*nband_k
           !  write(778,*)eig2_k(ii+(iband-1))
           !end do
           !write(778,*)gvnl1

           if(dtset%prtbbb==1)then
             call gaugetransfo(cg_k,gvnl1,cwavef_da,eig_k,eig2_k,iband,nband_k, &
&             dtset%mband,npw_k,npw1_k,dtset%nspinor,nsppol,occ_k)
             gvnl1(:,:) = cwavef_da(:,:)
           end if

         end if    !ddkfil(idir1)
       end if    !idir == idir1

!      <G|du/dqa> is contained in gvnl1 and <G|du/dqb> in cwavef
!      construct the matrix elements <du/dqa|du/dqb> -> dot
!      <u|du/dqa> -> dot1
!      <du/dqb|u> -> dot2
!      and add them to the 2nd-order matrix

       call dotprod_g(dotr,doti,istwf_k,npw1_k*dtset%nspinor,2,gvnl1,cwavef,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
       d2nl_k(1,idir1,ipert1)=d2nl_k(1,idir1,ipert1)+wtk_k*occ_k(iband)*dotr/(nband_kocc*two)
       d2nl_k(2,idir1,ipert1)=d2nl_k(2,idir1,ipert1)+wtk_k*occ_k(iband)*doti/(nband_kocc*two)


!      XG 020216 : Marek, could you check the next forty lines
!      In the parallel gauge, dot1 and dot2 vanishes
       if(dtset%prtbbb==1)then
         d2bbb_k(1,idir1,iband,iband)=d2bbb_k(1,idir1,iband,iband)+dotr
         d2bbb_k(2,idir1,iband,iband)=d2bbb_k(2,idir1,iband,iband)+doti
         dot_ndiagr=zero ; dot_ndiagi=zero
         do jband = 1,nband_k              !compute dot1 and dot2
           if (abs(occ_k(jband)) > tol8) then
             dot1r=zero ; dot1i=zero
             dot2r=zero ; dot2i=zero
             cwave0(:,:)=cg_k(:,1+(jband-1)*npw_k*dtset%nspinor:jband*npw_k*dtset%nspinor)

             call dotprod_g(dot1r,dot1i,istwf_k,npw1_k*dtset%nspinor,2,cwave0,gvnl1,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
             call dotprod_g(dot2r,dot2i,istwf_k,npw1_k*dtset%nspinor,2,cwavef,cwave0,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

             dot_ndiagr = dot_ndiagr + dot1r*dot2r - dot1i*dot2i
             dot_ndiagi = dot_ndiagi + dot1r*dot2i + dot1i*dot2r
             d2bbb_k(1,idir1,iband,jband) = d2bbb_k(1,idir1,iband,jband) - &
&             (dot1r*dot2r - dot1i*dot2i)
             d2bbb_k(2,idir1,iband,jband) = d2bbb_k(2,idir1,iband,jband) - &
&             (dot1r*dot2i + dot1i*dot2r)
           end if  ! occ_k
         end do !jband
         d2bbb_k(:,idir1,iband,:)=d2bbb_k(:,idir1,iband,:)*wtk_k*occ_k(iband)*half
         d2nl_k(1,idir1,ipert1)= &
&         d2nl_k(1,idir1,ipert1)-wtk_k*occ_k(iband)*dot_ndiagr/(nband_kocc*two)
         d2nl_k(2,idir1,ipert1)=&
&         d2nl_k(2,idir1,ipert1)-wtk_k*occ_k(iband)*dot_ndiagi/(nband_kocc*two)
       end if ! prtbbb==1

     end do  ! idir1
   end if   ! Compute localization tensor, ipert=natom+1

!  End loop over bands
 end do

!Final deallocations
 ABI_DEALLOCATE(cwave0)
 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(eig2_k)
 ABI_DEALLOCATE(gvnl1)
 ABI_DEALLOCATE(ffnlk)
 ABI_DEALLOCATE(ffnl1)
 ABI_DEALLOCATE(dkinpw)
 ABI_DEALLOCATE(kinpw1)
 ABI_DEALLOCATE(kpg_k)
 ABI_DEALLOCATE(kpg1_k)
 ABI_DEALLOCATE(ph3d)
 ABI_DATATYPE_DEALLOCATE(dum_cwaveprj)
 if(dtset%prtbbb==1)  then
   ABI_DEALLOCATE(cg_k)
   ABI_DEALLOCATE(cwavef_da)
   ABI_DEALLOCATE(cwavef_db)
 end if

 call timab(102,2,tsec)

 DBG_EXIT("COLL")

end subroutine dfpt_nstwf
!!***
