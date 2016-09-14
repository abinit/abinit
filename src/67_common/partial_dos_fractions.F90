!{\src2tex{textfont=tt}}
!!****f* ABINIT/partial_dos_fractions
!! NAME
!! partial_dos_fractions
!!
!! FUNCTION
!! calculate partial DOS fractions to feed to the tetrahedron method
!!  1: project states on angular momenta
!!  2: should be able to choose certain atoms or atom types, slabs of space...
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MVer,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  cg(2,mcg)=planewave coefficients of wavefunctions
!!  dtset     structured datatype, from which one uses :
!!   exchn2n3d=if 1, n2 and n3 are exchanged
!!   kpt(3,nkpt)  =irreducible kpoints
!!   kptrlatt(3,3)=lattice vectors for full kpoint grid
!!   mband        =maximum number of bands
!!   mkmem        =number of kpoints in memory
!!   natsph       =number of atoms for which the spherical decomposition must be done
!!   nband        =number of electronic bands for each kpoint
!!   nkpt         =number of irreducible kpoints
!!   nshiftk      =number of kpoint grid shifts
!!   nspinor      =number of spinor components
!!   nsppol       =1 or 2 spin polarization channels
!!   nsym         =number of symmetries
!!   shiftk(3,nshiftk)=kpoint shifts
!!   symrel(3,3,nsym)=symmetry matrices in real space
!!  mbesslang=maximum angular momentum for Bessel function expansion
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mpi_enreg=information about MPI parallelization
!!  m_dos_flag=option for the m-contributions to the partial DOS
!!  ndosfraction=natsph*mbesslang
!!  partial_dos= option for this routine - only 1 is supported at present
!!
!! OUTPUT
!!  dos_fractions(ikpt,iband,isppol,ndosfraction) = percentage of s, p, d..
!!    character on each atom for the wavefunction # ikpt,iband, isppol
!!  == if m_dos_flag==1
!!  dos_fractions_m(ikpt,iband,isppol,ndosfraction*mbesslang) = percentage of s, p, d..
!!    character on each atom for the wavefunction # ikpt,iband, isppol (m-resolved)
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      getkpgnorm,getph,getspin_1state,init_bess_spl,initylmg,kpgio,metric
!!      ph1d3d,recip_ylm,sort_dp,splint,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine partial_dos_fractions(crystal,npwarr,kg,cg,dos_fractions,dos_fractions_m,&
&           dtset,mbesslang,mcg,mpi_enreg,m_dos_flag,ndosfraction,partial_dos)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_splines
 use m_xmpi
 use m_crystal

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'partial_dos_fractions'
 use interfaces_28_numeric_noabirule
 use interfaces_41_geometry
 use interfaces_56_recipspace
 use interfaces_61_occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: m_dos_flag,mbesslang,mcg,ndosfraction,partial_dos
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: crystal
!arrays
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),npwarr(dtset%nkpt)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(out) :: dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction)
 real(dp),intent(out) :: dos_fractions_m(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*mbesslang*m_dos_flag)

!Local variables-------------------------------
!scalars
 integer,parameter :: prtsphere0=0 ! do not output all the band by band details for projections.
 integer :: shift_b,shift_sk,ia,iatom,iband,ierr,ikpt,ilang,ioffkg,is1, is2, isoff
 integer :: ipw,ispinor,isppol,ixint,mbess,mcg_disk,me,shift_cg
 integer :: mgfft,my_nspinor,n1,n2,n3,natsph_tot,nfit,npw_k,nradintmax
 integer :: comm
 real(dp) :: arg,bessarg,bessargmax,bessint_delta,kpgmax,rmax 
 character(len=500) :: msg
!arrays
 integer :: iindex(dtset%mpw)
 integer,allocatable :: iatsph(:),nradint(:),atindx(:)
 real(dp) :: kpoint(3),spin(3)
 real(dp) :: xfit(dtset%mpw),yfit(dtset%mpw)
 real(dp) :: ylm_k(dtset%mpw,mbesslang*mbesslang),ylmgr_dum(1)
 real(dp),allocatable :: bess_fit(:,:,:),bess_spl(:,:),bess_spl_der(:,:)
 real(dp),allocatable :: cg_1band(:,:),cg_1kpt(:,:),kpgnorm(:),ph1d(:,:)
 real(dp),allocatable :: ph3d(:,:,:),ratsph(:),rint(:),sum_1atom_1ll(:,:)
 real(dp),allocatable :: sum_1atom_1lm(:,:),x_bess(:),ylm(:,:)
 real(dp),allocatable :: xred_sph(:,:)
 real(dp),allocatable :: znucl_sph(:)
 real(dp),allocatable :: phkxred(:,:)
 real(dp),allocatable :: cmax(:)
 complex(dpc) :: cgcmat(2,2)

!*************************************************************************

!for the moment, only support projection on angular momenta
 if (partial_dos /= 1 .and. partial_dos /= 2) then
   write(std_out,*) 'Error : partial_dos_fractions only supports angular '
   write(std_out,*) ' momentum projection and spinor components for the moment. return to outscfcv'
   write(std_out,*) ' partial_dos = ', partial_dos
   return
 end if

!impose all kpoints have same number of bands
 do isppol=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
     if (dtset%nband((isppol-1)*dtset%nkpt + ikpt) /= dtset%mband) then
       write(std_out,*) 'Error : partial_dos_fractions wants same number of',' bands at each kpoint'
       write(std_out,*) ' isppol, ikpt = ', isppol,ikpt, dtset%nband((isppol-1)*dtset%nkpt + ikpt), dtset%mband
       write(std_out,*) ' all nband = ', dtset%nband
       return
     end if
   end do
 end do

!initialize dos_fractions
 dos_fractions(:,:,:,:) = zero
 if (m_dos_flag==1) dos_fractions_m(:,:,:,:) = zero

 if(mpi_enreg%paral_kgb==1) then
   comm = mpi_enreg%comm_kpt
   me = mpi_enreg%me_kpt
 else
   comm = mpi_enreg%comm_cell
   me=xmpi_comm_rank(comm)
 end if
 
 my_nspinor = max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 mcg_disk = dtset%mpw*my_nspinor*dtset%mband

 n1 = dtset%ngfft(1); n2 = dtset%ngfft(2); n3 = dtset%ngfft(3)
 mgfft = maxval(dtset%ngfft(1:3))

!##############################################################
!FIRST CASE: project on angular momenta to get dos parts
!##############################################################

 if (partial_dos == 1) then

   natsph_tot = dtset%natsph + dtset%natsph_extra
   ABI_ALLOCATE(iatsph,(natsph_tot))
   ABI_ALLOCATE(ratsph,(natsph_tot))
   ABI_ALLOCATE(znucl_sph,(natsph_tot))
   ABI_ALLOCATE(nradint,(natsph_tot))
   ABI_ALLOCATE(atindx,(natsph_tot))
   ABI_ALLOCATE(phkxred,(2,natsph_tot))
   ABI_ALLOCATE(cmax,(natsph_tot))

!  initialize atindx
   do iatom=1,natsph_tot
     atindx(iatom) = iatom
   end do

   iatsph(1:dtset%natsph)=dtset%iatsph(1:dtset%natsph)
!  random choice to set index of type for fictitious atoms to atom 1. Will become input variable as:
!  dtset%natsph_extra : int
!  dtset%ratsph_extra : real
!  dtset%xredsph_extra(3, dtset%natsph_extra) : real
   iatsph(dtset%natsph+1:natsph_tot)=1

   do iatom = 1, dtset%natsph
     ratsph(iatom) = dtset%ratsph(dtset%typat(iatsph(iatom)))
     znucl_sph(iatom) = dtset%znucl(dtset%typat(iatsph(iatom)))
   end do
   do iatom = 1, dtset%natsph_extra
     ratsph(iatom+dtset%natsph) = dtset%ratsph_extra
     znucl_sph(iatom+dtset%natsph) = 0
   end do

!  init bessel function integral for recip_ylm max ang mom + 1
   ABI_ALLOCATE(sum_1atom_1ll,(mbesslang,natsph_tot))
   ABI_ALLOCATE(sum_1atom_1lm,(mbesslang**2,natsph_tot))

   bessint_delta = 0.1_dp
   ! TODO: ecuteff instead of ecut?
   kpgmax = sqrt(dtset%ecut)
   rmax = zero; bessargmax = zero; nradintmax = 0
   do iatom=1,natsph_tot
     rmax=max(rmax,ratsph(iatom))
     bessarg=ratsph(iatom)*two_pi*kpgmax
     bessargmax=max(bessargmax,bessarg)
     nradint(iatom) = int (bessarg / bessint_delta) + 1
     nradintmax=max(nradintmax,nradint(iatom))
   end do
   !write(std_out,*) ' partial_dos_fractions :  rmax = ', rmax
!  use same number of grid points to calculate Bessel function and to do the integration later on r
   mbess = nradintmax
!  make sure bessargmax is a multiple of bessint_delta
   bessargmax = bessint_delta*mbess

   ABI_ALLOCATE(bess_spl,(mbess,mbesslang))
   ABI_ALLOCATE(bess_spl_der,(mbess,mbesslang))
   ABI_ALLOCATE(x_bess,(nradintmax))
   ABI_ALLOCATE(rint,(nradintmax))
   ABI_ALLOCATE(bess_fit,(dtset%mpw,nradintmax,mbesslang))

!  initialize general Bessel function array on uniform grid
!  x_bess, from 0 to (2 \pi |k+G|_{max} |r_{max}|)
!  TODO: Replace with: use m_paw_numeric, only : paw_jbessel_4spline, paw_spline
   call init_bess_spl(mbess,bessint_delta,mbesslang,bess_spl,bess_spl_der,x_bess)

!  get kg matrix of the positions of G vectors in recip space
!  for each electronic state, get corresponding wavefunction and project on Ylm
!  remember: cg(2,dtset%mpw*my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)

   ABI_ALLOCATE(xred_sph, (3, natsph_tot))
   do iatom=1,dtset%natsph
     xred_sph(:,iatom) = crystal%xred(:,iatsph(iatom))
   end do
   do iatom=1,dtset%natsph_extra
     xred_sph(:,dtset%natsph+iatom) = dtset%xredsph_extra(:, iatom)
   end do

   ABI_ALLOCATE(ph1d,(2,(2*n1+1 + 2*n2+1 + 2*n3+1)*natsph_tot))
   call getph(atindx,natsph_tot,n1,n2,n3,ph1d,xred_sph)

!  Now get Ylm factors: returns "real Ylms", which are real (+m) and
!  imaginary (-m) parts of actual complex Ylm. Yl-m = Ylm*
!  Single call to initylmg for all kg (all mkmem are in memory)
   ABI_ALLOCATE(ylm,(dtset%mpw*dtset%mkmem,mbesslang*mbesslang))
   call initylmg(crystal%gprimd,kg,dtset%kpt,dtset%mkmem,mpi_enreg,mbesslang,&
&   dtset%mpw,dtset%nband,dtset%nkpt,npwarr,dtset%nsppol,0,crystal%rprimd,ylm,ylmgr_dum)

!  kpgnorm contains norms only for kpoints used by this processor
   ABI_ALLOCATE(kpgnorm,(dtset%mpw*dtset%mkmem))
   kpgnorm (:) = zero
   ioffkg = 0
   do ikpt = 1, dtset%nkpt
     if (all(mpi_enreg%proc_distrb(ikpt,1,1:dtset%nsppol)/=me)) cycle
     npw_k = npwarr(ikpt)
     call getkpgnorm(crystal%gprimd, dtset%kpt(:,ikpt), kg(:,ioffkg+1:ioffkg+npw_k),&
&     kpgnorm(ioffkg+1:ioffkg+npw_k), npw_k)

     ioffkg = ioffkg + npw_k
   end do !ikpt

   shift_sk = 0
   do isppol=1,dtset%nsppol
!    kg array is the same for both sppol ?????
     ioffkg = 0

     do ikpt=1,dtset%nkpt
       if (all(mpi_enreg%proc_distrb(ikpt,:,isppol) /= mpi_enreg%me)) cycle
       kpoint(:) = dtset%kpt(:,ikpt)
       npw_k = npwarr(ikpt)

!      for each kpoint set up the phase factors, ylm factors
       do ilang=1,mbesslang*mbesslang
         do ipw=1,npw_k
           ylm_k(ipw,ilang) = ylm(ioffkg+ipw,ilang)
         end do
       end do

!      make phkred for all atoms
       do ia=1,natsph_tot
         arg=two_pi*( kpoint(1)*xred_sph(1,ia) + kpoint(2)*xred_sph(2,ia) + kpoint(3)*xred_sph(3,ia) )
         phkxred(1,ia)=cos(arg)
         phkxred(2,ia)=sin(arg)
       end do

!      get full phases for the following
       ABI_ALLOCATE(ph3d,(2,npw_k,natsph_tot))
       ph3d = zero
       call ph1d3d(1,natsph_tot,kg(:,ioffkg+1:ioffkg+npw_k), &
        natsph_tot,natsph_tot,npw_k,n1,n2,n3,phkxred,ph1d,ph3d)
!      phases exp (2 pi i (k+G).x_tau) are now in ph3d

!      get Bessel function factors on array of |k+G|*r distances
!      since we need many r distances and have a large number of different
!      |k+G|, get j_l on uniform grid (above, in array gen_besj),
!      and spline it for each kpt Gvector set.
!      TODO: Precompute (k+G) integrals here and pass them to recip_ylm)
!      
       nfit = npw_k
       do ixint=1,nradintmax
         rint(ixint) = (ixint-1)*rmax / (nradintmax-1)
         do ipw=1,npw_k
           xfit(ipw) = two_pi * kpgnorm(ipw+ioffkg) * rint(ixint)
           iindex(ipw) = ipw
         end do

         call sort_dp(npw_k,xfit,iindex,tol14)
         do ilang=1,mbesslang
           call splint(mbess,x_bess,bess_spl(:,ilang),bess_spl_der(:,ilang),nfit,xfit,yfit)
!          re-order results for different G vectors
           do ipw=1,npw_k
             bess_fit(iindex(ipw),ixint,ilang) = yfit(ipw)
           end do
         end do
       end do

       shift_b = 0
       do iband=1,dtset%mband
         if (mpi_enreg%proc_distrb(ikpt,iband,isppol)/=mpi_enreg%me) cycle

         do ispinor=1,my_nspinor
           ! Select wavefunction in cg array
           shift_cg = shift_sk + shift_b

           call recip_ylm (bess_fit, cg(:,shift_cg+1:shift_cg+npw_k), dtset%istwfk(ikpt),&
&           nradint, nradintmax,mbesslang, mpi_enreg, dtset%mpw, natsph_tot, npw_k,&
&           ph3d, prtsphere0, rint, ratsph, sum_1atom_1ll, sum_1atom_1lm,&
&           crystal%ucvol, ylm_k, znucl_sph)

           do iatom=1,natsph_tot
             do ilang=1,mbesslang
               dos_fractions(ikpt,iband,isppol,mbesslang*(iatom-1) + ilang) &
&               = dos_fractions(ikpt,iband,isppol,mbesslang*(iatom-1) + ilang) &
&               + sum_1atom_1ll(ilang,iatom)
             end do
           end do

           if (m_dos_flag==1) then
             do iatom=1,natsph_tot
               do ilang=1,mbesslang**2
                 dos_fractions_m(ikpt,iband,isppol,mbesslang**2*(iatom-1) + ilang) &
&                 = dos_fractions_m(ikpt,iband,isppol,mbesslang**2*(iatom-1) + ilang) &
&                 + sum_1atom_1lm(ilang,iatom)
               end do
             end do
           end if

           ! Increment band, spinor shift
           shift_b = shift_b + npw_k
         end do ! spinor
       end do ! band

       ! Increment kpt and (spin, kpt) shifts
       ioffkg = ioffkg + npw_k
       shift_sk = shift_sk + dtset%mband*my_nspinor*npw_k

       ABI_DEALLOCATE(ph3d)
     end do ! end kpt
   end do ! end sppol

!  Gather all contributions from different processors
   call xmpi_sum(dos_fractions,comm,ierr)
   if (m_dos_flag==1) then
     call xmpi_sum(dos_fractions_m,comm,ierr)
   end if
   if (mpi_enreg%paral_spinor == 1)then
     call xmpi_sum(dos_fractions,mpi_enreg%comm_spinor,ierr)
     if (m_dos_flag==1) then
       call xmpi_sum(dos_fractions_m,mpi_enreg%comm_spinor,ierr)
     end if
   end if

   ABI_DEALLOCATE(atindx)
   ABI_DEALLOCATE(bess_fit)
   ABI_DEALLOCATE(bess_spl)
   ABI_DEALLOCATE(bess_spl_der)
   ABI_DEALLOCATE(cmax)
   ABI_DEALLOCATE(iatsph)
   ABI_DEALLOCATE(kpgnorm)
   ABI_DEALLOCATE(nradint)
   ABI_DEALLOCATE(ph1d)
   ABI_DEALLOCATE(phkxred)
   ABI_DEALLOCATE(ratsph)
   ABI_DEALLOCATE(rint)
   ABI_DEALLOCATE(sum_1atom_1ll)
   ABI_DEALLOCATE(sum_1atom_1lm)
   ABI_DEALLOCATE(x_bess)
   ABI_DEALLOCATE(xred_sph)
   ABI_DEALLOCATE(ylm)
   ABI_DEALLOCATE(znucl_sph)

 else if (partial_dos == 2) then

   if (dtset%nsppol /= 1 .or. dtset%nspinor /= 2) then
     msg = "spinor projection of DOS is meaningless if nsppol==2 or nspinor/=2. Not calculating projections."
     MSG_WARNING(msg)
     dos_fractions = zero
     return
   end if
   if (my_nspinor /= 2) then
     msg = "spinor parallelization for this projection of DOS is not coded. Not calculating projections."
     MSG_WARNING(msg)
     dos_fractions = zero
     return
   end if

   ! FIXME: WHAT THE FUCK!
   ABI_ALLOCATE(cg_1kpt,(2,mcg_disk))
   shift_sk = 0
   isppol = 1

   do ikpt=1,dtset%nkpt
     if (all(mpi_enreg%proc_distrb(ikpt,:,isppol)/=mpi_enreg%me)) cycle

     cg_1kpt(:,:) = cg(:,shift_sk+1:shift_sk+mcg_disk)
     npw_k = npwarr(ikpt)
     ABI_ALLOCATE(cg_1band,(2,2*npw_k))
     shift_b=0
     do iband=1,dtset%mband
       if (mpi_enreg%proc_distrb(ikpt,iband,isppol)/=mpi_enreg%me) cycle

       cg_1band(:,:) = cg_1kpt(:,shift_b+1:shift_b+2*npw_k)
       call getspin_1state(cg_1band, npw_k, spin, cgcmat)

       do is1 = 1, 2
         do is2 = 1, 2
           isoff = is2 + (is1-1)*2
           dos_fractions(ikpt,iband,isppol,isoff) = dos_fractions(ikpt,iband,isppol,isoff) &
&           + real(cgcmat(is1,is2))
         end do
       end do

       dos_fractions(ikpt,iband,isppol,5) = dos_fractions(ikpt,iband,isppol,5) + spin(1)
       dos_fractions(ikpt,iband,isppol,6) = dos_fractions(ikpt,iband,isppol,6) + spin(2)
       dos_fractions(ikpt,iband,isppol,7) = dos_fractions(ikpt,iband,isppol,7) + spin(3)

       shift_b = shift_b + 2*npw_k
     end do
     ABI_DEALLOCATE(cg_1band)
     shift_sk = shift_sk + dtset%mband*2*npw_k
   end do
   ABI_DEALLOCATE(cg_1kpt)

!  Gather all contributions from different processors
   call xmpi_sum(dos_fractions,comm,ierr)
! below for future use - spinors should not be parallelized for the moment
!   if (mpi_enreg%paral_spinor == 1)then
!     call xmpi_sum(dos_fractions,mpi_enreg%comm_spinor,ierr)
!   end if

 else
   MSG_WARNING('only partial_dos==1 is coded')
 end if

end subroutine partial_dos_fractions
!!***
