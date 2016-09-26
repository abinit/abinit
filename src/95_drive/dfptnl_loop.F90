!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfptnl_loop
!! NAME
!! dfptnl_loop
!!
!! FUNCTION
!! Loop over the perturbations j1, j2 and j3
!!
!! COPYRIGHT
!! Copyright (C) 2016-2016 ABINIT group (LB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol) = array for planewave coefficients of wavefunctions
!!  cgindex(nkpt,nsppol) = for each k-point, cgindex tores the location
!!                         of the WF in the cg array
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  etotal = new total energy (no meaning at output)
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space(bohr^-1)
!!  gsqcut=Fourier cutoff on G^2 for "large sphere" of radius double
!!   that of the basis sphere--appropriate for charge density rho(G),
!!   Hartree potential, and pseudopotentials
!!  kg(3,mpw*mkmem)=reduced planewave coordinates
!!  kneigh(30,nkpt) = index of the neighbours of each k-point
!!  kg_neigh(30,nkpt,3) = necessary to construct the vector joining a k-point
!!                         to its nearest neighbour in case of a single k-point,
!!                         a line of k-points or a plane of k-points.
!!  kptindex(2,nkpt3)= index of the k-points in the reduced BZ
!!                     related to a k-point in the full BZ
!!  kpt3(3,nkpt3) = reduced coordinates of k-points in the full BZ
!!  kxc(nfftf,nkxc)=exchange-correlation kernel
!!  k3xc(nfftf,nk3xc)=third-order exchange-correlation kernel
!!  mband = maximum number of bands
!!  mgfft = maximum single fft dimension
!!  mkmem = Number of k points treated by this node.
!!  mkmem_max = maximal number of k-points on each processor (MPI //)
!!  mk1mem = Number of k points for first-order WF treated by this node.
!!  mpert =maximum number of ipert
!!  mpi_enreg=MPI-parallelisation information
!!  mpw   = maximum number of planewaves in basis sphere (large number)
!!  mvwtk(30,nkpt) = weights to compute the finite difference ddk
!!  natom = number of atoms in unit cell
!!  nfftf  = (effective) number of FFT grid points (for this processor)
!!  nkpt  = number of k points
!!  nkpt3 = number of k-points in the full BZ
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  nneigh  = total number of neighbours required to evaluate the finite
!!          difference formula
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  npwarr(nkpt) = array holding npw for each k point
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  pwind(mpw,nneigh,mkmem) = array used to compute the overlap matrix smat
!!                           between k-points
!!  rfpert(3,mpert,3,mpert,3,mpert) = array defining the type of perturbations
!!       that have to be computed
!!       1   ->   element has to be computed explicitely
!!      -1   ->   use symmetry operations to obtain the corresponding element
!!  rprimd(3,3)=dimensional primitive translations (bohr)
!!  ucvol = unit cell volume (bohr^3)
!!  xred(3,natom) = reduced atomic coordinates
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert) = flags for each element of the 3DTE
!!                             (=1 if computed)
!!  d3etot(2,3,mpert,3,mpert,3,mpert) = matrix of the 3DTEs
!!
!! SIDE EFFECTS
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!
!! PARENTS
!!      nonlinear
!!
!! CHILDREN
!!      appdig,dfpt_mkcore,dfpt_mkvxc,dfpt_vlocal,dfptnl_mv,dfptnl_resp
!!      dotprod_vn,fourdp,getph,hartre,initylmg,inwffil,read_rhor,status,timab
!!      wffclose,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfptnl_loop(blkflg,cg,cgindex,dtfil,dtset,d3etot,eigen0,etotal,gmet,gprimd,gsqcut, &
& hdr,kg,kneigh,kg_neigh,kptindex,kpt3,kxc,k3xc,mband,mgfft,mkmem,mkmem_max,mk1mem,&
& mpert,mpi_enreg,mpw,mvwtk,natom,nfftf,nkpt,nkpt3,nkxc,nk3xc,nneigh,nspinor,nsppol,&
& npwarr,occ,pawfgr,pawtab,psps,pwind,rfpert,rhog,rhor,rprimd,ucvol,usecprj,vtrial,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes

 use m_errors
 use m_profiling_abi
 use m_hdr
 use m_nctk
 use m_wffile
 use m_wfk

 use m_io_tools,    only : file_exists
 use m_ioarr,       only : read_rhor
 use m_hamiltonian, only : destroy_hamiltonian,destroy_rf_hamiltonian,gs_hamiltonian_type,&
                           init_hamiltonian,init_rf_hamiltonian,rf_hamiltonian_type
 use m_pawfgr,      only : pawfgr_type
 use m_pawrhoij,    only : pawrhoij_type
 use m_pawtab,      only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfptnl_loop'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_72_response
 use interfaces_79_seqpar_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mgfft,mk1mem,mkmem,mkmem_max,mpert,mpw,natom,nfftf
 integer,intent(in) :: nk3xc,nkpt,nkpt3,nkxc,nneigh,nspinor,nsppol,usecprj
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(inout) :: etotal
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps

!arrays
 integer,intent(in) :: cgindex(nkpt,nsppol),kg(3,mk1mem*mpw),kneigh(30,nkpt)
 integer,intent(in) :: kg_neigh(30,nkpt,3)
 integer,intent(in) :: kptindex(2,nkpt3),npwarr(nkpt),pwind(mpw,nneigh,mkmem)
 integer,intent(in) :: rfpert(3,mpert,3,mpert,3,mpert)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert,3,mpert) !vz_i
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),gmet(3,3)
 real(dp),intent(in) :: eigen0(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: gprimd(3,3),k3xc(nfftf,nk3xc),kpt3(3,nkpt3)
 real(dp),intent(in) :: kxc(nfftf,nkxc),mvwtk(30,nkpt),rhog(2,nfftf),rhor(nfftf,dtset%nspden),rprimd(3,3)
 real(dp),intent(in) :: vtrial(nfftf,dtset%nspden),xred(3,natom)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
 real(dp),intent(inout) :: d3etot(2,3,mpert,3,mpert,3,mpert) !vz_i
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=51
 integer :: ask_accurate,comm_cell,counter,cplex,formeig,i1dir
 integer :: i1pert,i2dir,i2pert,i3dir,i3pert,iatom,idir_dkde,ierr,iexit,ifft,ii,index,ir
 integer :: ireadwf,itypat,mcg,mpsang,n1,n2,n3,n3xccc,nfftotf,nhat1grdim,nspden,nwffile
 integer :: option,optene,optorth,pert1case,pert2case,pert3case
 integer :: rdwrpaw,timrev,usexcnhat
 real(dp) :: dummy_real,ecut_eff
 character(len=500) :: message
 character(len=fnlen) :: fiden1i,fiwf1i,fiwf3i,fiwfddk,fnamewff(3)
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 type(wffile_type) :: wff1,wff2,wff3,wfft1,wfft2,wfft3,wffddk(3)
 type(wfk_t) :: ddk_f(3)
 type(wvl_data) :: wvl
 type(hdr_type) :: hdr_den
!arrays
 integer,allocatable :: atindx(:),atindx1(:),nattyp(:)
 integer :: file_index(3)
 real(dp) :: rho_dum(1),tsec(2)
 real(dp),allocatable :: cg1(:,:),cg2(:,:),cg3(:,:),eigen1(:),eigen2(:),eigen3(:)
 real(dp),allocatable :: nhat(:,:),nhat1(:,:),nhat1gr(:,:,:),vresid_dum(:,:)
 real(dp),allocatable :: ph1d(:,:),rho1r1(:,:)
 real(dp),allocatable :: rho2g1(:,:),rho2r1(:,:),rho3r1(:,:),vhartr1(:)
 real(dp),allocatable :: vpsp1(:),vtrial1(:,:),vxc1(:,:),work(:),xc_tmp(:,:)
 real(dp),allocatable :: xccc3d1(:),xccc3d2(:),xccc3d3(:)
 type(pawrhoij_type),allocatable :: rhoij_dum(:)

! ***********************************************************************

 call timab(502,1,tsec)
 call status(0,dtfil%filstat,iexit,level,'enter         ')

 comm_cell = mpi_enreg%comm_cell

 timrev = 1
 cplex = 2 - timrev
 nspden = dtset%nspden
 ecut_eff = (dtset%ecut)*(dtset%dilatmx)**2
 mpsang = psps%mpsang
 optorth=1;if (psps%usepaw==1) optorth=0

 ABI_ALLOCATE(cg1,(2,dtset%mpw*dtset%nspinor*mband*dtset%mk1mem*dtset%nsppol))
 ABI_ALLOCATE(cg2,(2,dtset%mpw*dtset%nspinor*mband*dtset%mk1mem*dtset%nsppol))
 ABI_ALLOCATE(cg3,(2,dtset%mpw*dtset%nspinor*mband*dtset%mk1mem*dtset%nsppol))
 ABI_ALLOCATE(eigen1,(2*dtset%mband*dtset%mband*dtset%nkpt*dtset%nsppol))
 ABI_ALLOCATE(eigen2,(2*dtset%mband*dtset%mband*dtset%nkpt*dtset%nsppol))
 ABI_ALLOCATE(eigen3,(2*dtset%mband*dtset%mband*dtset%nkpt*dtset%nsppol))
 ABI_ALLOCATE(rho1r1,(cplex*nfftf,dtset%nspden))
 ABI_ALLOCATE(rho2r1,(cplex*nfftf,dtset%nspden))
 ABI_ALLOCATE(rho2g1,(2,nfftf))
 ABI_ALLOCATE(rho3r1,(cplex*nfftf,dtset%nspden))

 ask_accurate=1 ; formeig = 1 ; ireadwf = 1
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 nfftotf=n1*n2*n3

!Generate an index table of atoms, in order for them to be used
!type after type.
 ABI_ALLOCATE(atindx,(natom))
 ABI_ALLOCATE(atindx1,(natom))
 ABI_ALLOCATE(nattyp,(psps%ntypat))
 index=1
 do itypat=1,psps%ntypat
   nattyp(itypat)=0
   do iatom=1,natom
     if(dtset%typat(iatom)==itypat)then
       atindx(iatom)=index
       atindx1(index)=iatom
       index=index+1
       nattyp(itypat)=nattyp(itypat)+1
     end if
   end do
 end do

!Generate the 1-dimensional phases
 ABI_ALLOCATE(ph1d,(2,3*(2*mgfft+1)*natom))
 call getph(atindx,natom,n1,n2,n3,ph1d,xred)

!==== Initialize most of the Hamiltonian (and derivative) ====
!1) Allocate all arrays and initialize quantities that do not depend on k and spin.
!2) Perform the setup needed for the non-local factors:
!* Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
!* PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 call init_hamiltonian(gs_hamkq,psps,pawtab,dtset%nspinor,dtset%nspden,natom,&
& dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,&
& usecprj=usecprj,ph1d=ph1d,nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda)

 ABI_ALLOCATE(vpsp1,(cplex*nfftf))
 ABI_ALLOCATE(xccc3d1,(cplex*nfftf))
 ABI_ALLOCATE(xccc3d2,(cplex*nfftf))
 ABI_ALLOCATE(xccc3d3,(cplex*nfftf))
 ABI_ALLOCATE(vhartr1,(cplex*nfftf))
 ABI_ALLOCATE(vxc1,(cplex*nfftf,dtset%nspden))
 ABI_ALLOCATE(vtrial1,(cplex*nfftf,dtset%nspden))

 ABI_ALLOCATE(vresid_dum,(0,0))
! PAW stuff
 nhat1grdim = 0
 usexcnhat = 0
 ABI_ALLOCATE(nhat,(0,0))
 ABI_ALLOCATE(nhat1,(0,0))
 ABI_ALLOCATE(nhat1gr,(0,0,0))

 mcg=mpw*nspinor*mband*mkmem*nsppol

!Loop over the perturbations j1, j2, j3
 pert1case = 0 ; pert2case = 0 ; pert3case = 0
 do i1pert = 1, mpert
   do i1dir = 1, 3

     if ((maxval(rfpert(i1dir,i1pert,:,:,:,:))==1)) then

       pert1case = i1dir + (i1pert-1)*3
       counter = pert1case
       call appdig(pert1case,dtfil%fnamewff1,fiwf1i)

       call status(counter,dtfil%filstat,iexit,level,'call inwffil  ')

       call inwffil(ask_accurate,cg1,dtset,dtset%ecut,ecut_eff,eigen1,dtset%exchn2n3d,&
&       formeig,gmet,hdr,&
&       ireadwf,dtset%istwfk,kg,dtset%kptns,dtset%localrdwf,&
&       dtset%mband,mcg,dtset%mk1mem,mpi_enreg,mpw,&
&       dtset%nband,dtset%ngfft,dtset%nkpt,npwarr,&
&       dtset%nsppol,dtset%nsym,&
&       occ,optorth,rprimd,&
&       dtset%symafm,dtset%symrel,dtset%tnons,&
&       dtfil%unkg1,wff1,wfft1,dtfil%unwff1,fiwf1i,wvl)

       if (ireadwf==1) then
         call WffClose (wff1,ierr)
       end if

       rho1r1(:,:) = zero
       if (dtset%get1den /= 0 .or. dtset%ird1den /= 0) then
         rdwrpaw=0
         call appdig(pert1case,dtfil%fildens1in,fiden1i)
         call status(counter,dtfil%filstat,iexit,level,'call ioarr    ')

         call read_rhor(fiden1i, cplex, dtset%nspden, nfftf, dtset%ngfft, rdwrpaw, mpi_enreg, rho1r1, &
         hdr_den, rhoij_dum, comm_cell, check_hdr=hdr)
         etotal = hdr_den%etot; call hdr_free(hdr_den)
       end if

       xccc3d1(:) = zero
       if ((psps%n1xccc/=0).and.(i1pert <= natom)) then
         call status(counter,dtfil%filstat,iexit,level,'call dfpt_mkcore   ')
         call dfpt_mkcore(cplex,i1dir,i1pert,natom,psps%ntypat,n1,psps%n1xccc,&
&         n2,n3,dtset%qptn,rprimd,dtset%typat,ucvol,&
&         psps%xcccrc,psps%xccc1d,xccc3d1,xred)
       end if ! psps%n1xccc/=0

       do i3pert = 1, mpert
         do i3dir = 1, 3

           if ((maxval(rfpert(i1dir,i1pert,:,:,i3dir,i3pert))==1)) then

             pert3case = i3dir + (i3pert-1)*3
             counter = 100*pert3case + pert1case
             call appdig(pert3case,dtfil%fnamewff1,fiwf3i)

             call status(counter,dtfil%filstat,iexit,level,'call inwffil  ')
             call inwffil(ask_accurate,cg3,dtset,dtset%ecut,ecut_eff,eigen3,dtset%exchn2n3d,&
&             formeig,gmet,hdr,&
&             ireadwf,dtset%istwfk,kg,dtset%kptns,dtset%localrdwf,&
&             dtset%mband,mcg,dtset%mk1mem,mpi_enreg,mpw,&
&             dtset%nband,dtset%ngfft,dtset%nkpt,npwarr,&
&             dtset%nsppol,dtset%nsym,&
&             occ,optorth,rprimd,&
&             dtset%symafm,dtset%symrel,dtset%tnons,&
&             dtfil%unkg1,wff3,wfft3,dtfil%unwff3,&
&             fiwf3i,wvl)
             if (ireadwf==1) then
               call WffClose (wff3,ierr)
             end if

             rho3r1(:,:) = zero
             if (dtset%get1den /= 0 .or. dtset%ird1den /= 0) then
               rdwrpaw=0
               call appdig(pert3case,dtfil%fildens1in,fiden1i)
               call status(counter,dtfil%filstat,iexit,level,'call ioarr    ')

               call read_rhor(fiden1i, cplex, dtset%nspden, nfftf, dtset%ngfft, rdwrpaw, mpi_enreg, rho3r1, &
               hdr_den, rhoij_dum, comm_cell, check_hdr=hdr)
               etotal = hdr_den%etot; call hdr_free(hdr_den)
             end if

             xccc3d3(:) = zero
             if ((psps%n1xccc/=0).and.(i3pert <= natom)) then
               call status(counter,dtfil%filstat,iexit,level,'call dfpt_mkcore   ')
               call dfpt_mkcore(cplex,i3dir,i3pert,natom,psps%ntypat,n1,psps%n1xccc,&
&               n2,n3,dtset%qptn,rprimd,dtset%typat,ucvol,&
&               psps%xcccrc,psps%xccc1d,xccc3d3,xred)
             end if ! psps%n1xccc/=0

             do i2pert = 1, mpert

!              In case of electric field perturbation, evaluate the ddk
!              using the finite difference expression of
!              Marzari and Vanderbilt PRB 56, 12847 (1997).

!               d3_berry(:,:) = zero

!               if ((i2pert==dtset%natom+2).and.&
!&               (maxval(rfpert(i1dir,i1pert,:,i2pert,i3dir,i3pert)) == 1)) then

!                 call timab(511,1,tsec)
!                 call status(counter,dtfil%filstat,iexit,level,'call dfptnl_mv  ')
!                 call dfptnl_mv(cg,cgindex,cg1,cg3,dtset,dtfil,d3_berry,gmet,&
!&                 i1pert,i3pert,i1dir,i3dir,&
!&                 kneigh,kg_neigh,kptindex,kpt3,mband,mkmem,mkmem_max,mk1mem,&
!&                 mpi_enreg,mpw,mvwtk,natom,nkpt,nkpt3,nneigh,npwarr,nspinor,nsppol,pwind)
!                 call timab(511,2,tsec)

!               end if

               if (mpi_enreg%me == 0) then

                 if(sum(rfpert(i1dir,i1pert,:,i2pert,i3dir,i3pert))>0)then
                   write(message,'(a,a,a,a,a,a)')ch10,ch10,&
&                   ' Decomposition of the third-order energy for the set of perturbations',ch10
                   call wrtout(std_out,message,'COLL')
                   call wrtout(ab_out,message,'COLL')
                   if (i1pert < natom + 1) then
                     write(message,'(a,i3,a,i3)') &
&                     ' j1 : displacement of atom ',i1pert,' along direction ', i1dir
                   end if
                   if (i1pert == dtset%natom + 2) then
                     write(message,'(a,i4)')' j1 : homogeneous electric field along direction ',i1dir
                   end if
                   call wrtout(std_out,message,'COLL')
                   call wrtout(ab_out,message,'COLL')
                   if (i3pert < natom + 1) then
                     write(message,'(a,i3,a,i3,a)') &
&                     ' j3 : displacement of atom ',i3pert,' along direction ', i3dir,ch10
                   end if
                   if (i3pert == dtset%natom + 2) then
                     write(message,'(a,i4,a)')' j3 : homogeneous electric field along direction ',i3dir,ch10
                   end if
                   call wrtout(std_out,message,'COLL')
                   call wrtout(ab_out,message,'COLL')
                 end if

               end if ! mpi_enreg%me == 0

               do i2dir = 1, 3

                 if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
                   pert2case = i2dir + (i2pert-1)*3

                   blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 1

                   call status(counter,dtfil%filstat,iexit,level,'call inwffil  ')
                   call inwffil(ask_accurate,cg2,dtset,dtset%ecut,ecut_eff,eigen2,dtset%exchn2n3d,&
&                   formeig,gmet,hdr,&
&                   ireadwf,dtset%istwfk,kg,dtset%kptns,dtset%localrdwf,&
&                   dtset%mband,mcg,dtset%mk1mem,mpi_enreg,mpw,&
&                   dtset%nband,dtset%ngfft,dtset%nkpt,npwarr,&
&                   dtset%nsppol,dtset%nsym,&
&                   occ,optorth,rprimd,&
&                   dtset%symafm,dtset%symrel,dtset%tnons,&
&                   dtfil%unkg1,wff2,wfft2,dtfil%unwff2,&
&                   fiwf3i,wvl)
                   if (ireadwf==1) then
                     call WffClose (wff2,ierr)
                   end if

!                  Read the first-order densities from disk-files
                   rho2r1(:,:) = zero ; rho2g1(:,:) = zero

                   if (dtset%get1den /= 0 .or. dtset%ird1den /= 0) then
                     rdwrpaw=0
                     call appdig(pert2case,dtfil%fildens1in,fiden1i)
                     call status(counter,dtfil%filstat,iexit,level,'call ioarr    ')

                     call read_rhor(fiden1i, cplex, dtset%nspden, nfftf, dtset%ngfft, rdwrpaw, mpi_enreg, rho2r1, &
                     hdr_den, rhoij_dum, comm_cell, check_hdr=hdr)
                     etotal = hdr_den%etot; call hdr_free(hdr_den)

!                    Compute up+down rho1(G) by fft
                     ABI_ALLOCATE(work,(cplex*nfftf))
                     work(:)=rho2r1(:,1)
                     call status(counter,dtfil%filstat,iexit,level,'call fourdp   ')
                     call fourdp(cplex,rho2g1,work,-1,mpi_enreg,nfftf,dtset%ngfft,dtset%paral_kgb,0)
                     ABI_DEALLOCATE(work)

                   end if

!                  Compute first-order local potentials
!                  (hartree, xc and pseudopotential)

                   n3xccc=0; if(psps%n1xccc/=0)n3xccc=nfftf
                   xccc3d2(:)=zero ; vpsp1(:)=zero

                   if (i2pert <= natom) then

                     call status(counter,dtfil%filstat,iexit,level,'call dfpt_vlocal   ')
                     call dfpt_vlocal(atindx,cplex,gmet,gsqcut,i2dir,i2pert,mpi_enreg,psps%mqgrid_vl,natom,&
&                     nattyp,nfftf,dtset%ngfft,psps%ntypat,n1,n2,n3,dtset%paral_kgb,ph1d,psps%qgrid_vl,&
&                     dtset%qptn,ucvol,psps%vlspl,vpsp1,xred)

                     if (psps%n1xccc/=0) then
                       call status(counter,dtfil%filstat,iexit,level,'call dfpt_mkcore   ')
                       call dfpt_mkcore(cplex,i2dir,i2pert,natom,psps%ntypat,n1,psps%n1xccc,&
&                       n2,n3,dtset%qptn,rprimd,dtset%typat,ucvol,&
&                       psps%xcccrc,psps%xccc1d,xccc3d2,xred)
                     end if ! psps%n1xccc/=0

                   end if  ! i2pert <= natom

                   call status(counter,dtfil%filstat,iexit,level,'get vtrial1   ')
                   option=1;optene=0
                   call dfpt_rhotov(cplex,dummy_real,dummy_real,dummy_real,dummy_real,gmet,gprimd,&
&                   gsqcut,i2dir,i2pert,dtset%ixc,kxc,mpi_enreg,dtset%natom,nfftf,dtset%ngfft,nhat,&
&                   nhat1,nhat1gr,nhat1grdim,nkxc,nspden,n3xccc,optene,option,dtset%paral_kgb,&
&                   dtset%qptn,rhog,rho2g1,rhor,rho2r1,rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,&
&                   vpsp1,vresid_dum,dummy_real,vtrial1,vxc1,xccc3d2)

                   call init_rf_hamiltonian(cplex,gs_hamkq,i2pert,rf_hamkq,has_e1kbsc=1)

                   nwffile = 1
                   file_index(1) = i2dir + 3*(i2pert-1)
                   fnamewff(1) = dtfil%fnamewff1

                   if (i2pert==natom+2) then

                     nwffile = 3
                     file_index(2) = i2dir+natom*3
                     idir_dkde = i2dir
                     if (i3dir/=i2dir) then ! see m_rf2.F90 => getidirs
                       if (i2dir==2.and.i3dir==3) idir_dkde = 4
                       if (i2dir==1.and.i3dir==3) idir_dkde = 5
                       if (i2dir==1.and.i3dir==2) idir_dkde = 6
                       if (i2dir==3.and.i3dir==2) idir_dkde = 7
                       if (i2dir==3.and.i3dir==1) idir_dkde = 8
                       if (i2dir==2.and.i3dir==1) idir_dkde = 9
                     end if
                     file_index(3) = idir_dkde+9+(dtset%natom+6)*3
                     fnamewff(2) = dtfil%fnamewffddk
                     fnamewff(3) = dtfil%fnamewffdkde

                   end if

                   do ii=1,nwffile
                     call appdig(file_index(ii),fnamewff(ii),fiwfddk)
                     ! Checking the existence of data file
                     if (.not. file_exists(fiwfddk)) then
                       ! Trick needed to run Abinit test suite in netcdf mode.
                       if (file_exists(nctk_ncify(fiwfddk))) then
                         write(message,"(3a)")"- File: ",trim(fiwfddk),&
                         " does not exist but found netcdf file with similar name."
                         call wrtout(std_out,message,'COLL')
                         fiwfddk = nctk_ncify(fiwfddk)
                       end if
                       if (.not. file_exists(fiwfddk)) then
                         MSG_ERROR('Missing file: '//TRIM(fiwfddk))
                       end if
                     end if
                     write(message,'(2a)')'-dfptnl_loop : read the wavefunctions from file: ',trim(fiwfddk)
                     call wrtout(std_out,message,'COLL')
                     call wrtout(ab_out,message,'COLL')
#ifdef DEV_MG_WFK
!                    Note that the unit number for these files is 50,51,52 or 53 (dtfil%unddk=50)
                     call wfk_open_read(ddk_f(ii),fiwfddk,1,dtset%iomode,dtfil%unddk+(ii-1),mpi_enreg%comm_cell)
#else
                     call WffOpen(dtset%iomode,mpi_enreg%comm_cell,fiwfddk,ierr,wffddk(ii),master,me,dtfil%unddk+(ii-1))
#endif
                   end do

!                  Perform DFPT part of the 3dte calculation
                   call timab(512,1,tsec)
                   call status(counter,dtfil%filstat,iexit,level,'call dfptnl_resp ')
!                  NOTE : eigen2 equals zero here
                   call dfptnl_pert(cg,cg1,cg3,cplex,dtfil,dtset,d3etot,eigen0,gs_hamkq,k3xc,i1dir,&
&                   i2dir,i3dir,i1pert,i2pert,i3pert,kg,mband,mgfft,mkmem,mk1mem,mpert,mpi_enreg,&
&                   mpsang,mpw,natom,nfftf,nfftotf,nkpt,nk3xc,nspden,nspinor,nsppol,npwarr,occ,pawfgr,ph1d,psps,&
&                   rf_hamkq,rho1r1,rho2r1,rho3r1,rprimd,ucvol,vtrial,vtrial1,wffddk,ddk_f,&
&                   xccc3d1,xccc3d2,xccc3d3,xred)
                   call timab(512,2,tsec)

                   call status(counter,dtfil%filstat,iexit,level,'after dfptnl_resp')

!!                  Describe the perturbation and write out the result
!                   if (mpi_enreg%me == 0) then
!                     if (i2pert < natom + 1) then
!                       write(message,'(a,i3,a,i3)') &
!&                       ' j2 : displacement of atom ',i2pert,&
!&                       ' along direction ', i2dir
!                     end if
!                     if (i2pert == dtset%natom + 2) then
!                       write(message,'(a,i4)') &
!&                       ' j2 : homogeneous electric field along direction ',&
!&                       i2dir
!                     end if
!                     call wrtout(std_out,message,'COLL')
!                     call wrtout(ab_out,message,'COLL')
!                     write(ab_out,'(20x,a,13x,a)')'real part','imaginary part'
!                     write(ab_out,'(5x,a2,1x,f22.10,3x,f22.10)')'xc',exc3*sixth,zero
!                     write(ab_out,'(5x,a3,f22.10,3x,f22.10)')'dft',&
!&                     d3etot(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert),&
!&                     d3etot(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
!                     write(ab_out,*)
!                     write(std_out,'(18x,a,11x,a)')'real part','imaginary part'
!                     write(std_out,'(5x,a2,1x,f20.10,3x,f20.10)')'xc',exc3*sixth,zero
!                     write(std_out,'(5x,a3,f22.10,3x,f22.10)')'dft',&
!&                     d3etot(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert),&
!&                     d3etot(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
!                     write(std_out,*)
!                   end if  ! mpi_enreg%me == 0

!                   d3etot(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = &
!&                   d3etot(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) + exc3*sixth
!                   d3etot(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = &
!&                   d3etot(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)

!                  Eventually close the dot file
#ifdef DEV_MG_WFK
                   call wfk_close(ddk_f(1))
#else
                   call WffClose(wffddk(1),ierr)
#endif
                   if (i2pert==dtset%natom+2) then
#ifdef DEV_MG_WFK
                     call wfk_close(ddk_f(2))
                     call wfk_close(ddk_f(3)) ! TO CHANGE
#else
                     call WffClose(wffddk(2),ierr)
                     call WffClose(wffddk(3),ierr) ! TO CHANGE
#endif
                   end if

                 end if   !rfpert
               end do    !i2dir
             end do     ! i2pert

           end if   ! rfpert
         end do    ! i3dir
       end do     ! i3pert

     end if   ! rfpert
   end do    ! i1dir
 end do     ! i1pert

 call status(0,dtfil%filstat,iexit,level,'exit          ')

!More memory cleaning
 call destroy_hamiltonian(gs_hamkq)
 call destroy_rf_hamiltonian(rf_hamkq)
 
 ABI_DEALLOCATE(cg1)
 ABI_DEALLOCATE(cg2)
 ABI_DEALLOCATE(cg3)
 ABI_DEALLOCATE(eigen1)
 ABI_DEALLOCATE(eigen2)
 ABI_DEALLOCATE(eigen3)
 ABI_DEALLOCATE(rho1r1)
 ABI_DEALLOCATE(rho2r1)
 ABI_DEALLOCATE(rho2g1)
 ABI_DEALLOCATE(rho3r1)
 ABI_DEALLOCATE(atindx1)
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(nattyp)
 ABI_DEALLOCATE(nhat)
 ABI_DEALLOCATE(nhat1)
 ABI_DEALLOCATE(nhat1gr)
 ABI_DEALLOCATE(ph1d)
 ABI_DEALLOCATE(vresid_dum)
 ABI_DEALLOCATE(vtrial1)
 ABI_DEALLOCATE(vxc1)
 ABI_DEALLOCATE(vhartr1)
 ABI_DEALLOCATE(vpsp1)
 ABI_DEALLOCATE(xccc3d1)
 ABI_DEALLOCATE(xccc3d2)
 ABI_DEALLOCATE(xccc3d3)


 call timab(502,2,tsec)

end subroutine dfptnl_loop
!!***
