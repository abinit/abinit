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
!! Copyright (C) 1998-2018 ABINIT group (MVer,MB,MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  crystal<crystal_t>= data type gathering info on symmetries and unit cell
!!  dtset<type(dataset_type)>=all input variables for this dataset
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  cg(2,mcg)=planewave coefficients of wavefunctions
!!  mcg=size of wave-functions array (cg) =mpw*my_nspinor*mband*mkmem*nsppol
!!  collect=1 if fractions should be MPI collected at the end, 0 otherwise.
!!  mpi_enreg=information about MPI parallelization
!!
!! SIDE EFFECTS
!!  dos%fractions(ikpt,iband,isppol,ndosfraction) = percentage of s, p, d..
!!    character on each atom for the wavefunction # ikpt,iband, isppol
!!  == if prtdosm /= 0
!!  dos%fractions_m(ikpt,iband,isppol,ndosfraction*mbesslang) = percentage of s, p, d..
!!    character on each atom for the wavefunction # ikpt,iband, isppol (m-resolved)
!!
!! NOTES
!!
!!   psi(r) = (4pi/sqrt(ucvol)) \sum_{LMG} i**l u(G) e^{i(k+G).Ra} x Y_{LM}^*(k+G) Y_{LM}(r-Ra) j_L(|k+G||r-Ra|)
!!
!!   int_(ratsph) |psi(r)|**2 = \sum_LM rho(LM)
!!
!!   where
!!
!!   rho_{LM} = 4pi \int_o^{rc} dr r**2 ||\sum_G u(G) Y_LM^*(k+G) e^{i(k+G).Ra} j_L(|k+G| r)||**2
!!
!!   where S is a RSH. The final expression is:
!!
!!   When k = G0/2, we have u_{G0/2}(G) = u_{G0/2}(-G-G0)^* and P can be rewritten as
!!
!!     P = (4pi i^L}/sqrt(ucvol) \sum^'_G w(G) S_{LM}(k+G) \int_0^ratsph dr r^2 j_L(|k+G|r) x
!!                                  2 Re[u_k(G) e^{i(k+G).R_atom}]  if L = 2n
!!                                  2 Im[u_k(G) e^{i(k+G).R_atom}]  if L = 2n + 1
!!
!!  where the sum over G is done on the reduced G-sphere and w(G) = 1/2 if G=G0 else 1.
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      cg_getspin,cwtime,destroy_mpi_enreg,getkpgnorm,getph,initmpi_seq
!!      initylmg,jlspline_free,ph1d3d,recip_ylm,sort_dp,splint,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine partial_dos_fractions(dos,crystal,dtset,eigen,occ,npwarr,kg,cg,mcg,collect,mpi_enreg)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors
 use m_splines
 use m_xmpi
 use m_mpinfo
 use m_crystal
 use m_sort

 use m_time,          only : cwtime
 use m_numeric_tools, only : simpson
 use m_special_funcs, only : jlspline_t, jlspline_new, jlspline_free, jlspline_integral
 use m_cgtools,       only : cg_getspin
 use m_epjdos,        only : recip_ylm, epjdos_t
 use m_io_tools,      only : get_unit
 use m_fstrings,      only : int2char4
 use m_gsphere,       only : getkpgnorm
 use m_kg,            only : ph1d3d, getph

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'partial_dos_fractions'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_51_manage_mpi
 use interfaces_56_recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcg,collect
 type(epjdos_t),intent(inout) :: dos
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: crystal
!arrays
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),npwarr(dtset%nkpt)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: prtsphere0=0 ! do not output all the band by band details for projections.
 integer :: shift_b,shift_sk,iat,iatom,iband,ierr,ikpt,ilang,ioffkg,is1, is2, isoff
 integer :: ipw,isppol,ixint,mbess,mcg_disk,me_kpt,shift_cg
 integer :: mgfft,my_nspinor,n1,n2,n3,natsph_tot,npw_k,nradintmax
 integer :: rc_ylm,itypat,nband_k
 integer :: abs_shift_b
 integer :: unit_procar, ipauli
 real(dp),parameter :: bessint_delta = 0.1_dp
 real(dp) :: arg,bessarg,bessargmax,kpgmax,rmax
 real(dp) :: cpu,wall,gflops
 character(len=500) :: msg
 character(len=4) :: ikproc_str
 character(len=fnlen) :: filename
 type(jlspline_t) :: jlspl
 type(MPI_type) :: mpi_enreg_seq
!arrays
 integer :: iindex(dtset%mpw),nband_tmp(1),npwarr_tmp(1)
 integer,allocatable :: iatsph(:),nradint(:),atindx(:),typat_extra(:),kg_k(:,:)
 real(dp) :: kpoint(3),spin(3),ylmgr_dum(1)
 real(dp) :: xfit(dtset%mpw),yfit(dtset%mpw)
 real(dp),allocatable :: ylm_k(:,:)
 real(dp),allocatable :: bess_fit(:,:,:)
 real(dp),allocatable :: cg_1band(:,:),cg_1kpt(:,:),kpgnorm(:),ph1d(:,:)
 real(dp),allocatable :: ph3d(:,:,:),ratsph(:),rint(:),sum_1atom_1ll(:,:,:)
 real(dp),allocatable :: sum_1atom_1lm(:,:,:)
 real(dp),allocatable :: xred_sph(:,:),znucl_sph(:),phkxred(:,:)
 complex(dpc) :: cgcmat(2,2)

!*************************************************************************

 ! for the moment, only support projection on angular momenta
 if (dos%partial_dos_flag /= 1 .and. dos%partial_dos_flag /= 2) then
   write(std_out,*) 'Error: partial_dos_fractions only supports angular '
   write(std_out,*) ' momentum projection and spinor components for the moment. return to outscfcv'
   write(std_out,*) ' partial_dos = ', dos%partial_dos_flag
   return
 end if

 ! impose all kpoints have same number of bands
 do isppol=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
     if (dtset%nband((isppol-1)*dtset%nkpt + ikpt) /= dtset%mband) then
       write(std_out,*) 'Error: partial_dos_fractions wants same number of bands at each kpoint'
       write(std_out,*) ' isppol, ikpt = ', isppol,ikpt, dtset%nband((isppol-1)*dtset%nkpt + ikpt), dtset%mband
       write(std_out,*) ' all nband = ', dtset%nband
       return
     end if
   end do
 end do

 ! Real or complex spherical harmonics?
 rc_ylm = 2; if (dos%prtdosm == 2) rc_ylm = 1

 me_kpt = mpi_enreg%me_kpt
 my_nspinor = max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

 n1 = dtset%ngfft(1); n2 = dtset%ngfft(2); n3 = dtset%ngfft(3)
 mgfft = maxval(dtset%ngfft(1:3))

 call cwtime(cpu, wall, gflops, "start")

! open file for each proc, and print header for master node
 unit_procar=get_unit()
 call int2char4(me_kpt, ikproc_str)
 filename = 'PROCAR_'//ikproc_str
 open(file=filename, unit=unit_procar)
 if(mpi_enreg%me==0) then
   write (unit_procar,'(a)') 'PROCAR lm decomposed - need to concatenate files in parallel case'
   write (unit_procar,'(a,I10,a,I10,a,I10,a)') '# of k-points: ', dtset%nkpt, &
&     ' # of bands:', dtset%mband, ' # of ions:', dtset%natom, ch10
 end if

!##############################################################
!FIRST CASE: project on angular momenta to get dos parts
!##############################################################

 if (dos%partial_dos_flag == 1) then
   natsph_tot = dtset%natsph + dtset%natsph_extra

   ABI_ALLOCATE(iatsph, (natsph_tot))
   ABI_ALLOCATE(typat_extra, (natsph_tot))
   ABI_ALLOCATE(ratsph, (natsph_tot))
   ABI_ALLOCATE(znucl_sph, (natsph_tot))
   ABI_ALLOCATE(nradint, (natsph_tot))
   ABI_ALLOCATE(atindx, (natsph_tot))
   ABI_ALLOCATE(phkxred, (2,natsph_tot))

   ! initialize atindx
   do iatom=1,natsph_tot
     atindx(iatom) = iatom
   end do

   iatsph(1:dtset%natsph) = dtset%iatsph(1:dtset%natsph)
   do iatom=1,dtset%natsph
     itypat = dtset%typat(iatsph(iatom))
     typat_extra(iatom) = itypat
     ratsph(iatom) = dtset%ratsph(itypat)
     znucl_sph(iatom) = dtset%znucl(itypat)
   end do

   ! fictitious atoms are declared with
   ! %natsph_extra, %ratsph_extra and %xredsph_extra(3, dtset%natsph_extra)
   ! they have atom index (natom + ii) and itype = ntypat + 1
   do iatom=1,dtset%natsph_extra
     typat_extra(iatom+dtset%natsph) = dtset%ntypat + 1
     ratsph(iatom+dtset%natsph) = dtset%ratsph_extra
     znucl_sph(iatom+dtset%natsph) = zero
     iatsph(iatom+dtset%natsph) = dtset%natom + iatom
   end do

   ! init bessel function integral for recip_ylm max ang mom + 1
   ABI_ALLOCATE(sum_1atom_1ll,(dtset%nspinor**2,dos%mbesslang,natsph_tot))
   ABI_ALLOCATE(sum_1atom_1lm,(dtset%nspinor**2,dos%mbesslang**2,natsph_tot))

   ! Note ecuteff instead of ecut.
   kpgmax = sqrt(dtset%ecut * dtset%dilatmx**2)
   rmax = zero; bessargmax = zero; nradintmax = 0
   do iatom=1,natsph_tot
     rmax = max(rmax, ratsph(iatom))
     bessarg = ratsph(iatom)*two_pi*kpgmax
     bessargmax = max(bessargmax, bessarg)
     nradint(iatom) = int (bessarg / bessint_delta) + 1
     nradintmax = max(nradintmax,nradint(iatom))
   end do
   !write(std_out,*)' partial_dos_fractions: rmax=', rmax,' nradintmax: ", nradintmax
!  use same number of grid points to calculate Bessel function and to do the integration later on r
!  and make sure bessargmax is a multiple of bessint_delta
   mbess = nradintmax
   bessargmax = bessint_delta*mbess

   ABI_ALLOCATE(rint,(nradintmax))
   ABI_ALLOCATE(bess_fit,(dtset%mpw,nradintmax,dos%mbesslang))

   ! initialize general Bessel function array on uniform grid xx, from 0 to (2 \pi |k+G|_{max} |r_{max}|)
   jlspl = jlspline_new(mbess, bessint_delta, dos%mbesslang)

   ABI_ALLOCATE(xred_sph, (3, natsph_tot))
   do iatom=1,dtset%natsph
     xred_sph(:,iatom) = crystal%xred(:,iatsph(iatom))
   end do
   do iatom=1,dtset%natsph_extra
     xred_sph(:,dtset%natsph+iatom) = dtset%xredsph_extra(:, iatom)
   end do

   ABI_ALLOCATE(ph1d,(2,(2*n1+1 + 2*n2+1 + 2*n3+1)*natsph_tot))
   call getph(atindx,natsph_tot,n1,n2,n3,ph1d,xred_sph)

   ! Fake MPI data to be used for sequential call to initylmg.
   call initmpi_seq(mpi_enreg_seq)
   mpi_enreg_seq%my_natom = dtset%natom

   shift_sk = 0
   abs_shift_b =  0 ! offset to allow for automatic update with +1 below
   do isppol=1,dtset%nsppol
     ioffkg = 0
     do ikpt=1,dtset%nkpt
       nband_k = dtset%nband((isppol-1)*dtset%nkpt + ikpt)
       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_kpt)) then
         abs_shift_b = abs_shift_b + nband_k ! jump the whole kpt in the eig and occ arrays
         cycle
       end if
       npw_k = npwarr(ikpt)
       kpoint(:) = dtset%kpt(:,ikpt)
   
       write (unit_procar,'(a,I7,a,3F12.6,a,F12.6,a)') ' k-point ', ikpt, ' : ', kpoint(:), ' weight = ', dtset%wtk(ikpt), ch10

       ! make phkred for all atoms
       do iat=1,natsph_tot
         arg=two_pi*( kpoint(1)*xred_sph(1,iat) + kpoint(2)*xred_sph(2,iat) + kpoint(3)*xred_sph(3,iat) )
         phkxred(1,iat)=cos(arg)
         phkxred(2,iat)=sin(arg)
       end do

       ABI_MALLOC(kg_k, (3, npw_k))
       kg_k = kg(:,ioffkg+1:ioffkg+npw_k)

       ! kpgnorm contains norms only for kpoints used by this processor
       ABI_ALLOCATE(kpgnorm, (npw_k))
       call getkpgnorm(crystal%gprimd, kpoint, kg_k, kpgnorm, npw_k)

       ! Now get Ylm(k, G) factors: returns "real Ylms", which are real (+m) and
       ! imaginary (-m) parts of actual complex Ylm. Yl-m = Ylm*
       ! Call initylmg for a single k-point (mind mpi_enreg_seq).
       ABI_MALLOC(ylm_k, (npw_k, dos%mbesslang**2))
       npwarr_tmp(1) = npw_k; nband_tmp(1) = nband_k
       call initylmg(crystal%gprimd,kg_k,kpoint,1,mpi_enreg_seq,dos%mbesslang,&
       npw_k,nband_tmp,1,npwarr_tmp,1,0,crystal%rprimd,ylm_k,ylmgr_dum)

       ! get phases exp (2 pi i (k+G).x_tau) in ph3d
       ABI_ALLOCATE(ph3d,(2,npw_k,natsph_tot))
       call ph1d3d(1,natsph_tot,kg_k,natsph_tot,natsph_tot,npw_k,n1,n2,n3,phkxred,ph1d,ph3d)

       ! get Bessel function factors on array of |k+G|*r distances
       ! since we need many r distances and have a large number of different
       ! |k+G|, get j_l on uniform grid (above, in array gen_besj),
       ! and spline it for each kpt Gvector set.
       ! Note that we use the same step based on rmax, this can lead to (hopefully small)
       ! inaccuracies when we integrate from 0 up to rmax(iatom)
       do ixint=1,nradintmax
         rint(ixint) = (ixint-1)*rmax / (nradintmax-1)
         do ipw=1,npw_k
           xfit(ipw) = two_pi * kpgnorm(ipw) * rint(ixint)
           iindex(ipw) = ipw
         end do

         call sort_dp(npw_k,xfit,iindex,tol14)
         do ilang=1,dos%mbesslang
           call splint(mbess, jlspl%xx, jlspl%bess_spl(:,ilang), jlspl%bess_spl_der(:,ilang), npw_k, xfit, yfit)
           ! re-order results for different G vectors
           do ipw=1,npw_k
             bess_fit(iindex(ipw),ixint,ilang) = yfit(ipw)
           end do
         end do
       end do ! ixint

       shift_b = 0
       do iband=1,nband_k
         abs_shift_b = abs_shift_b + 1
         if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband,isppol,me_kpt)) cycle
         !write(std_out,*)"in band:",iband
         ! TODO: eventually import eig and occ down to here - a pain, but printing outside would imply saving a huge array in memory
         write (unit_procar,'(a,I7,a,F12.6,a,F12.6,a)') 'band ', iband, ' # energy ', &
&          eigen(abs_shift_b), ' # occ. ', occ(abs_shift_b), ch10

         ! Select wavefunction in cg array
         shift_cg = shift_sk + shift_b

         call recip_ylm(bess_fit, cg(:,shift_cg+1:shift_cg+my_nspinor*npw_k), dtset%istwfk(ikpt),&
&          mpi_enreg, nradint, nradintmax, dos%mbesslang , dtset%mpw, natsph_tot, typat_extra, dos%mlang_type,&
&          npw_k, dtset%nspinor, ph3d, prtsphere0, rint, ratsph, rc_ylm, sum_1atom_1ll, sum_1atom_1lm,&
&          crystal%ucvol, ylm_k, znucl_sph)
         ! on exit the sum_1atom_* have both spinors counted

         ! Accumulate
         do iatom=1,natsph_tot
           do ilang=1,dos%mbesslang
             dos%fractions(ikpt,iband,isppol,dos%mbesslang*(iatom-1) + ilang) &
&             = dos%fractions(ikpt,iband,isppol,dos%mbesslang*(iatom-1) + ilang) &
&             + sum_1atom_1ll(1,ilang,iatom)
           end do
         end do

         if (dos%prtdosm /= 0) then
           do iatom=1,natsph_tot
             do ilang=1,dos%mbesslang**2
               dos%fractions_m(ikpt,iband,isppol,dos%mbesslang**2*(iatom-1) + ilang) &
&               = dos%fractions_m(ikpt,iband,isppol,dos%mbesslang**2*(iatom-1) + ilang) &
&               + sum_1atom_1lm(1,ilang,iatom)
             end do
           end do
         end if

         ! Increment band, spinor shift
         !shift_b = shift_b + npw_k
         shift_b = shift_b + my_nspinor*npw_k

! now we have both spinor components.
         write (unit_procar,'(a)') 'ion      s     py     pz     px    dxy    dyz    dz2    dxz    dx2    tot'
         do ipauli= 1,dtset%nspinor**2
!Contract with Pauli matrices to get projections for this k and band, all atoms and ilang
           do iatom = 1, natsph_tot
             write (unit_procar, '(I3)', advance='no') iatom
             do ilang=1,min(dos%mbesslang**2,9)
               write (unit_procar, '(F7.3)',advance='no') sum_1atom_1lm(ipauli,ilang,iatom)
             end do
             write (unit_procar, '(F7.3)',advance='yes') sum(sum_1atom_1lm(ipauli,:,iatom))
           end do
         end do
         write (unit_procar,*)


       end do ! band

       ! Increment kpt and (spin, kpt) shifts
       ioffkg = ioffkg + npw_k
       shift_sk = shift_sk + nband_k*my_nspinor*npw_k

       ABI_FREE(kg_k)
       ABI_FREE(kpgnorm)
       ABI_FREE(ylm_k)
       ABI_DEALLOCATE(ph3d)
     end do ! ikpt
   end do ! isppol

   ! collect = 1 ==> gather all contributions from different processors
   if (collect == 1) then
     call xmpi_sum(dos%fractions,mpi_enreg%comm_kpt,ierr)
     if (dos%prtdosm /= 0) call xmpi_sum(dos%fractions_m,mpi_enreg%comm_kpt,ierr)

! this is now done inside recip_ylm
!     if (mpi_enreg%paral_spinor == 1)then
!       call xmpi_sum(dos%fractions,mpi_enreg%comm_spinor,ierr)
!       if (dos%prtdosm /= 0) call xmpi_sum(dos%fractions_m,mpi_enreg%comm_spinor,ierr)
!     end if
   end if

   ABI_DEALLOCATE(atindx)
   ABI_DEALLOCATE(bess_fit)
   ABI_DEALLOCATE(iatsph)
   ABI_DEALLOCATE(typat_extra)
   ABI_DEALLOCATE(nradint)
   ABI_DEALLOCATE(ph1d)
   ABI_DEALLOCATE(phkxred)
   ABI_DEALLOCATE(ratsph)
   ABI_DEALLOCATE(rint)
   ABI_DEALLOCATE(sum_1atom_1ll)
   ABI_DEALLOCATE(sum_1atom_1lm)
   ABI_DEALLOCATE(xred_sph)
   ABI_DEALLOCATE(znucl_sph)

   call jlspline_free(jlspl)
   call destroy_mpi_enreg(mpi_enreg_seq)

 !##############################################################
 !2ND CASE: project on spinors
 !##############################################################

 else if (dos%partial_dos_flag == 2) then

   if (dtset%nsppol /= 1 .or. dtset%nspinor /= 2) then
     MSG_WARNING("spinor projection is meaningless if nsppol==2 or nspinor/=2. Not calculating projections.")
     return
   end if
   if (my_nspinor /= 2) then
     MSG_WARNING("spinor projection with spinor parallelization is not coded. Not calculating projections.")
     return
   end if
   ABI_CHECK(mpi_enreg%paral_spinor == 0, "prtdos 5 does not support spinor parallelism")

   ! FIXME: We should not allocate such a large chunk of memory!
   mcg_disk = dtset%mpw*my_nspinor*dtset%mband
   ABI_ALLOCATE(cg_1kpt,(2,mcg_disk))
   shift_sk = 0
   isppol = 1

   do ikpt=1,dtset%nkpt
     nband_k = dtset%nband((isppol-1)*dtset%nkpt + ikpt)
     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_kpt)) cycle
     npw_k = npwarr(ikpt)

     cg_1kpt(:,:) = cg(:,shift_sk+1:shift_sk+mcg_disk)
     ABI_ALLOCATE(cg_1band,(2,2*npw_k))
     shift_b=0
     do iband=1,nband_k
       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband,isppol,me_kpt)) cycle

       ! Select wavefunction in cg array
       !shift_cg = shift_sk + shift_b
       cg_1band(:,:) = cg_1kpt(:,shift_b+1:shift_b+2*npw_k)
       call cg_getspin(cg_1band, npw_k, spin, cgcmat=cgcmat)

       ! MG: TODO: imag part of off-diagonal terms is missing.
       ! I will add them later on.
       do is1=1,2
         do is2=1,2
           isoff = is2 + (is1-1)*2
           dos%fractions(ikpt,iband,isppol,isoff) = dos%fractions(ikpt,iband,isppol,isoff) &
&           + real(cgcmat(is1,is2))
         end do
       end do

       dos%fractions(ikpt,iband,isppol,5) = dos%fractions(ikpt,iband,isppol,5) + spin(1)
       dos%fractions(ikpt,iband,isppol,6) = dos%fractions(ikpt,iband,isppol,6) + spin(2)
       dos%fractions(ikpt,iband,isppol,7) = dos%fractions(ikpt,iband,isppol,7) + spin(3)

       shift_b = shift_b + 2*npw_k
     end do
     ABI_DEALLOCATE(cg_1band)
     shift_sk = shift_sk + nband_k*2*npw_k
   end do
   ABI_DEALLOCATE(cg_1kpt)

   ! Gather all contributions from different processors
   if (collect == 1) then
     call xmpi_sum(dos%fractions,mpi_enreg%comm_kpt,ierr)
     call xmpi_sum(dos%fractions,mpi_enreg%comm_bandfft,ierr)
     !below for future use - spinors should not be parallelized for the moment
     !if (mpi_enreg%paral_spinor == 1)then
     !  call xmpi_sum(dos%fractions,mpi_enreg%comm_spinor,ierr)
     !end if
   end if

 else
   MSG_WARNING('only partial_dos==1 or 2 is coded')
 end if

 close(unit_procar)

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2),a)')" partial_dos_fractions: cpu_time: ",cpu,"[s], walltime: ",wall," [s]"
 call wrtout(std_out,msg,"PERS")

end subroutine partial_dos_fractions
!!***
