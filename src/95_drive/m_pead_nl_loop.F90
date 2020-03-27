!!****m* ABINIT/m_pead_nl_loop
!! NAME
!!  m_pead_nl_loop
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2002-2020 ABINIT group (MVeithen,MB)
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

module m_pead_nl_loop

 use defs_basis
 use defs_wvltypes
 use m_wffile
 use m_abicore
 use m_xmpi
 use m_hdr
 use m_dtset
 use m_dtfil
#if defined HAVE_MPI2
 use mpi
#endif

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes, only : MPI_type
 use m_time,     only : timab
 use m_kg,       only : getph, mkkpg
 use m_cgtools,  only : dotprod_vn, dotprod_g
 use m_fft,      only : fourdp, fftpac, fourwf
 use m_ioarr,    only : read_rhor
 use m_pawtab,   only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type
 use m_pawcprj,    only : pawcprj_type
 use m_inwffil,  only : inwffil
 use m_spacepar, only : hartre
 use m_initylmg, only : initylmg
 use m_dfpt_mkvxc, only : dfpt_mkvxc
 use m_mkcore,     only : dfpt_mkcore
 use m_mklocl,     only : dfpt_vlocal
 use m_hamiltonian,only : init_hamiltonian, gs_hamiltonian_type
 use m_mkffnl,     only : mkffnl
 use m_mpinfo,     only : proc_distrb_cycle
 use m_nonlop,     only : nonlop

 implicit none

 private

#if defined HAVE_MPI1
 include 'mpif.h'
#endif
!!***

 public :: pead_nl_loop
!!***

contains
!!***

!!****f* ABINIT/pead_nl_loop
!! NAME
!! pead_nl_loop
!!
!! FUNCTION
!! Loop over the perturbations j1, j2 and j3
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol) = array for planewave coefficients of wavefunctions
!!  cgindex(nkpt,nsppol) = for each k-point, cgindex tores the location
!!                         of the WF in the cg array
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
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
!!  kxc(nfft,nkxc)=exchange-correlation kernel
!!  k3xc(nfft,nk3xc)=third-order exchange-correlation kernel
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
!!  nfft  = (effective) number of FFT grid points (for this processor)
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
!!  d3lo(2,3,mpert,3,mpert,3,mpert) = matrix of the 3DTEs
!!
!! SIDE EFFECTS
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!
!! PARENTS
!!      nonlinear
!!
!! CHILDREN
!!      appdig,dfpt_mkcore,dfpt_mkvxc,dfpt_vlocal,pead_nl_mv,pead_nl_resp
!!      dotprod_vn,fourdp,getph,hartre,initylmg,inwffil,read_rhor,status,timab
!!      wffclose,wrtout
!!
!! SOURCE

subroutine pead_nl_loop(blkflg,cg,cgindex,dtfil,dtset,d3lo,&
& gmet,gprimd,gsqcut, &
& hdr,kg,kneigh,kg_neigh,kptindex,kpt3,kxc,k3xc,mband,mgfft,mkmem,mkmem_max,mk1mem,&
& mpert,mpi_enreg,mpw,mvwtk,natom,nfft,nkpt,nkpt3,nkxc,nk3xc,nneigh,nspinor,nsppol,&
& npwarr,occ,psps,pwind,&
& rfpert,rprimd,ucvol,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mgfft,mk1mem,mkmem,mkmem_max,mpert,mpw,natom,nfft
 integer,intent(in) :: nk3xc,nkpt,nkpt3,nkxc,nneigh,nspinor,nsppol
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(hdr_type),intent(inout) :: hdr
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: cgindex(nkpt,nsppol),kg(3,mk1mem*mpw),kneigh(30,nkpt)
 integer,intent(in) :: kg_neigh(30,nkpt,3)
 integer,intent(in) :: kptindex(2,nkpt3),npwarr(nkpt),pwind(mpw,nneigh,mkmem)
 integer,intent(in) :: rfpert(3,mpert,3,mpert,3,mpert)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert,3,mpert) !vz_i
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),gmet(3,3)
 real(dp),intent(in) :: gprimd(3,3),k3xc(nfft,nk3xc),kpt3(3,nkpt3)
 real(dp),intent(in) :: kxc(nfft,nkxc),mvwtk(30,nkpt),rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol)
 real(dp),intent(inout) :: d3lo(2,3,mpert,3,mpert,3,mpert) !vz_i

!Local variables-------------------------------
!scalars
 integer,parameter :: level=51
 integer :: ask_accurate,counter,cplex,formeig,i1dir
 integer :: i1pert,i2dir,i2pert,i3dir,i3pert,iatom,ierr,ifft,index,ir
 integer :: ireadwf,itypat,mcg,mpsang,n1,n2,n3,n3xccc,nfftot,nspden,option,optorth
 integer :: pert1case,pert2case,pert3case,rdwrpaw,timrev,comm_cell
 logical :: nmxc
 real(dp) :: ecut_eff,exc3,valuei
 character(len=500) :: message
 character(len=fnlen) :: fiden1i,fiwf1i,fiwf3i
 type(wffile_type) :: wff1,wff2,wfft1,wfft2
 type(wvl_data) :: wvl
 type(hdr_type) :: hdr_den
!arrays
 integer,allocatable :: atindx(:),atindx1(:),nattyp(:)
 real(dp) :: d3_berry(2,3),rho_dum(1),tsec(2),ylmgr_dum(1)
 real(dp),allocatable :: cg1(:,:),cg3(:,:),eigen1(:),ph1d(:,:),rho1r1(:,:)
 real(dp),allocatable :: rho2g1(:,:),rho2r1(:,:),rho3r1(:,:),vhartr1(:)
 real(dp),allocatable :: vpsp1(:),vtrial1(:,:),vxc1(:,:),work(:),xc_tmp(:,:)
 real(dp),allocatable :: xccc3d1(:),xccc3d2(:),xccc3d3(:),ylm(:,:,:)
 type(pawrhoij_type),allocatable :: rhoij_dum(:)

! ***********************************************************************

 call timab(502,1,tsec)

 comm_cell = mpi_enreg%comm_cell

 timrev = 1
 cplex = 2 - timrev
 nspden = dtset%nspden
 ecut_eff = (dtset%ecut)*(dtset%dilatmx)**2
 mpsang = psps%mpsang
 optorth=1;if (psps%usepaw==1) optorth=0

 ABI_ALLOCATE(cg1,(2,dtset%mpw*dtset%nspinor*mband*dtset%mk1mem*dtset%nsppol))
 ABI_ALLOCATE(cg3,(2,dtset%mpw*dtset%nspinor*mband*dtset%mk1mem*dtset%nsppol))
 ABI_ALLOCATE(eigen1,(2*dtset%mband*dtset%mband*dtset%nkpt*dtset%nsppol))
 ABI_ALLOCATE(rho1r1,(cplex*nfft,dtset%nspden))
 ABI_ALLOCATE(rho2r1,(cplex*nfft,dtset%nspden))
 ABI_ALLOCATE(rho2g1,(2,nfft))
 ABI_ALLOCATE(rho3r1,(cplex*nfft,dtset%nspden))
 ABI_ALLOCATE(ylm,(2,dtset%mpw*dtset%mkmem,mpsang*mpsang*psps%useylm))

 ask_accurate=1 ; formeig = 1 ; ireadwf = 1
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 nfftot=n1*n2*n3

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

!Set up the Ylm for each k point
 if (psps%useylm==1) then
   call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,psps%mpsang,&
&   dtset%mpw,dtset%nband,dtset%nkpt,&
&   npwarr,dtset%nsppol,0,rprimd,ylm,ylmgr_dum)
 end if

 ABI_ALLOCATE(vpsp1,(cplex*nfft))
 ABI_ALLOCATE(xccc3d1,(cplex*nfft))
 ABI_ALLOCATE(xccc3d2,(cplex*nfft))
 ABI_ALLOCATE(xccc3d3,(cplex*nfft))
 ABI_ALLOCATE(vhartr1,(cplex*nfft))
 ABI_ALLOCATE(vxc1,(cplex*nfft,dtset%nspden))
 ABI_ALLOCATE(vtrial1,(cplex*nfft,dtset%nspden))

!Loop over the perturbations j1, j2, j3

 pert1case = 0 ; pert2case = 0 ; pert3case = 0

 do i1pert = 1, mpert
   do i1dir = 1, 3

     if ((maxval(rfpert(i1dir,i1pert,:,:,:,:))==1)) then

       pert1case = i1dir + (i1pert-1)*3
       counter = pert1case
       call appdig(pert1case,dtfil%fnamewff1,fiwf1i)

       mcg=mpw*nspinor*mband*mkmem*nsppol
       call inwffil(ask_accurate,cg1,dtset,dtset%ecut,ecut_eff,eigen1,dtset%exchn2n3d,&
&       formeig,hdr,ireadwf,dtset%istwfk,kg,dtset%kptns,dtset%localrdwf,&
&       dtset%mband,mcg,dtset%mk1mem,mpi_enreg,mpw,&
&       dtset%nband,dtset%ngfft,dtset%nkpt,npwarr,&
&       dtset%nsppol,dtset%nsym,&
&       occ,optorth,dtset%symafm,dtset%symrel,dtset%tnons,&
&       dtfil%unkg1,wff1,wfft1,dtfil%unwff1,fiwf1i,wvl)

       if (ireadwf==1) then
         call WffClose (wff1,ierr)
       end if

       rho1r1(:,:) = 0._dp
       if (dtset%get1den /= 0 .or. dtset%ird1den /= 0) then
         rdwrpaw=0
         call appdig(pert1case,dtfil%fildens1in,fiden1i)

         call read_rhor(fiden1i, cplex, dtset%nspden, nfft, dtset%ngfft, rdwrpaw, mpi_enreg, rho1r1, &
         hdr_den, rhoij_dum, comm_cell, check_hdr=hdr)
         call hdr_den%free()
       end if

       xccc3d1(:) = 0._dp
       if ((psps%n1xccc/=0).and.(i1pert <= natom)) then
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

             mcg=mpw*nspinor*mband*mkmem*nsppol
             call inwffil(ask_accurate,cg3,dtset,dtset%ecut,ecut_eff,eigen1,dtset%exchn2n3d,&
&             formeig,hdr,ireadwf,dtset%istwfk,kg,dtset%kptns,dtset%localrdwf,&
&             dtset%mband,mcg,dtset%mk1mem,mpi_enreg,mpw,&
&             dtset%nband,dtset%ngfft,dtset%nkpt,npwarr,&
&             dtset%nsppol,dtset%nsym,&
&             occ,optorth,dtset%symafm,dtset%symrel,dtset%tnons,&
&             dtfil%unkg1,wff2,wfft2,dtfil%unwff2,&
&             fiwf3i,wvl)
             if (ireadwf==1) then
               call WffClose (wff2,ierr)
             end if

             rho3r1(:,:) = 0._dp
             if (dtset%get1den /= 0 .or. dtset%ird1den /= 0) then
               rdwrpaw=0
               call appdig(pert3case,dtfil%fildens1in,fiden1i)

               call read_rhor(fiden1i, cplex, dtset%nspden, nfft, dtset%ngfft, rdwrpaw, mpi_enreg, rho3r1, &
               hdr_den, rhoij_dum, comm_cell, check_hdr=hdr)
               call hdr_den%free()
             end if

             xccc3d3(:) = 0._dp
             if ((psps%n1xccc/=0).and.(i3pert <= natom)) then
               call dfpt_mkcore(cplex,i3dir,i3pert,natom,psps%ntypat,n1,psps%n1xccc,&
&               n2,n3,dtset%qptn,rprimd,dtset%typat,ucvol,&
&               psps%xcccrc,psps%xccc1d,xccc3d3,xred)
             end if ! psps%n1xccc/=0

             do i2pert = 1, mpert

!              In case of electric field perturbation, evaluate the ddk
!              using the finite difference expression of
!              Marzari and Vanderbilt PRB 56, 12847 (1997) [[cite:Marzari1997]].

               d3_berry(:,:) = 0._dp

               if ((i2pert==dtset%natom+2).and.&
&               (maxval(rfpert(i1dir,i1pert,:,i2pert,i3dir,i3pert)) == 1)) then

                 call timab(511,1,tsec)
                 call pead_nl_mv(cg,cgindex,cg1,cg3,dtset,dtfil,d3_berry,gmet,&
&                 i1pert,i3pert,i1dir,i3dir,&
&                 kneigh,kg_neigh,kptindex,kpt3,mband,mkmem,mkmem_max,mk1mem,&
&                 mpi_enreg,mpw,mvwtk,natom,nkpt,nkpt3,nneigh,npwarr,nspinor,nsppol,pwind)
                 call timab(511,2,tsec)

               end if

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

!                  Read the first-order densities from disk-files
                   rho2r1(:,:) = 0._dp ; rho2g1(:,:) = 0._dp

                   if (dtset%get1den /= 0 .or. dtset%ird1den /= 0) then
                     rdwrpaw=0
                     call appdig(pert2case,dtfil%fildens1in,fiden1i)

                     call read_rhor(fiden1i, cplex, dtset%nspden, nfft, dtset%ngfft, rdwrpaw, mpi_enreg, rho2r1, &
                     hdr_den, rhoij_dum, comm_cell, check_hdr=hdr)
                     call hdr_den%free()

!                    Compute up+down rho1(G) by fft
                     ABI_ALLOCATE(work,(cplex*nfft))
                     work(:)=rho2r1(:,1)
                     call fourdp(cplex,rho2g1,work,-1,mpi_enreg,nfft,1,dtset%ngfft,0)
                     ABI_DEALLOCATE(work)

                   end if

!                  Compute first-order local potentials
!                  (hartree, xc and pseudopotential)

                   n3xccc=0; if(psps%n1xccc/=0)n3xccc=nfft
                   xccc3d2(:)=0._dp ; vpsp1(:)=0._dp

                   if (i2pert <= natom) then

                     call dfpt_vlocal(atindx,cplex,gmet,gsqcut,i2dir,i2pert,mpi_enreg,psps%mqgrid_vl,natom,&
&                     nattyp,nfft,dtset%ngfft,psps%ntypat,n1,n2,n3,ph1d,psps%qgrid_vl,&
&                     dtset%qptn,ucvol,psps%vlspl,vpsp1,xred)

                     if (psps%n1xccc/=0) then
                       call dfpt_mkcore(cplex,i2dir,i2pert,natom,psps%ntypat,n1,psps%n1xccc,&
&                       n2,n3,dtset%qptn,rprimd,dtset%typat,ucvol,&
&                       psps%xcccrc,psps%xccc1d,xccc3d2,xred)
                     end if ! psps%n1xccc/=0

                   end if  ! i2pert <= natom

                   call hartre(cplex,gsqcut,0,mpi_enreg,nfft,dtset%ngfft,rho2g1,rprimd,vhartr1)
                   option=1 ; nmxc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)
                   call dfpt_mkvxc(cplex,dtset%ixc,kxc,mpi_enreg,nfft,dtset%ngfft,&
&                   rho_dum,0,rho_dum,0,nkxc,nmxc,dtset%nspden,n3xccc,option,&
&                   dtset%qptn,rho2r1,rprimd,0,vxc1,xccc3d2)

                   if(dtset%nsppol==1)then
                     if(cplex==1)then
                       do ir=1,nfft
                         vtrial1(ir,1)=vpsp1(ir)+vhartr1(ir)+vxc1(ir,1)
                       end do
                     else
                       do ir=1,nfft
                         vtrial1(2*ir-1,1)=vpsp1(2*ir-1)+vhartr1(2*ir-1)+vxc1(2*ir-1,1)
                         vtrial1(2*ir  ,1)=vpsp1(2*ir  )+vhartr1(2*ir  )+vxc1(2*ir  ,1)
                       end do
                     end if
                   else
                     if(cplex==1)then
                       do ir=1,nfft
                         vtrial1(ir,1)=vpsp1(ir)+vhartr1(ir)+vxc1(ir,1)
                         vtrial1(ir,2)=vpsp1(ir)+vhartr1(ir)+vxc1(ir,2)
                       end do
                     else
!                      fab: I think there was an error in the definition of  vtrial1(2*ir-1,2); I have corrected it...
                       do ir=1,nfft
                         vtrial1(2*ir-1,1)=vpsp1(2*ir-1)+vhartr1(2*ir-1)+vxc1(2*ir-1,1)
                         vtrial1(2*ir  ,1)=vpsp1(2*ir  )+vhartr1(2*ir  )+vxc1(2*ir  ,1)
                         vtrial1(2*ir-1,2)=vpsp1(2*ir-1)+vhartr1(2*ir-1)+vxc1(2*ir-1  ,2)
                         vtrial1(2*ir  ,2)=vpsp1(2*ir  )+vhartr1(2*ir  )+vxc1(2*ir  ,2)
                       end do
                     end if
                   end if

!                  Compute the third-order xc energy
!                  take into account the contribution of the term
!$
!                  \frac{d}{d \lambda}
!                  \frac{\delta^2 E_{Hxc}}{\delta n(r) \delta n(r\prim)}
!$
!                  (seventh term of Eq. (110) of X. Gonze, PRA 52, 1096 (1995) [[cite:Gonze1995]]).

!                  the following are essentially the 4th and the 3rd terms of PRB 71,125107 [[cite:Veithen2005]], but the
!                  multiplication for rho1 will be done by dotprod_vn later

!                  in the non spin polarized case xc_tmp has only 1 component
                   if (nspden==1)then

                     ABI_ALLOCATE(xc_tmp,(cplex*nfft,1))

                     if (cplex==1) then
!                      This, and the next lines, have to be changed in case cplex=2
                       do ifft=1,nfft
                         xc_tmp(ifft,1)= k3xc(ifft,1)*(rho2r1(ifft,1)+3*xccc3d2(ifft))*rho3r1(ifft,1)
                       end do
                     else
                       do ifft=1,nfft   ! 2*ifft-1 denotes the real part, 2*ifft the imaginary part
                         xc_tmp(2*ifft-1,1)= k3xc(ifft,1)*( (rho2r1(2*ifft-1,1)+3*xccc3d2(2*ifft-1))*rho3r1(2*ifft-1,1) &
&                         -( rho2r1(2*ifft,1)+3*xccc3d2(2*ifft))*rho3r1(2*ifft,1))

                         xc_tmp(2*ifft,1)= k3xc(ifft,1)*( (rho2r1(2*ifft-1,1)+3*xccc3d2(2*ifft-1))*rho3r1(2*ifft,1) &
&                         +( rho2r1(2*ifft,1)+3*xccc3d2(2*ifft))*rho3r1(2*ifft-1,1))
                       end do

                     end if

                   end if

!                  fab: modifications for the spin polarized raman part:
!                  in the spin polarized case xc_tmp has 2 components
!                  note that now the non linear core correction is divided by 2
                   if (nspden==2) then

                     ABI_ALLOCATE(xc_tmp,(cplex*nfft,2))

                     if (cplex==1) then
                       do ifft=1,nfft
                         xc_tmp(ifft,1)= k3xc(ifft,1)*(rho2r1(ifft,2)+(3._dp/2._dp)*xccc3d2(ifft))*rho3r1(ifft,2)+ &
&                         k3xc(ifft,2)*(rho2r1(ifft,2)+(3._dp/2._dp)*xccc3d2(ifft))*(rho3r1(ifft,1)-rho3r1(ifft,2))+ &
&                         k3xc(ifft,2)*((rho2r1(ifft,1)-rho2r1(ifft,2))+(3._dp/2._dp)*xccc3d2(ifft))*rho3r1(ifft,2)+ &
&                         k3xc(ifft,3)*((rho2r1(ifft,1)-rho2r1(ifft,2))+(3._dp/2._dp)*xccc3d2(ifft))*(rho3r1(ifft,1)-rho3r1(ifft,2))
                         xc_tmp(ifft,2)= k3xc(ifft,2)*(rho2r1(ifft,2)+(3._dp/2._dp)*xccc3d2(ifft))*rho3r1(ifft,2)+ &
&                         k3xc(ifft,3)*(rho2r1(ifft,2)+(3._dp/2._dp)*xccc3d2(ifft))*(rho3r1(ifft,1)-rho3r1(ifft,2))+ &
&                         k3xc(ifft,3)*((rho2r1(ifft,1)-rho2r1(ifft,2))+(3._dp/2._dp)*xccc3d2(ifft))*rho3r1(ifft,2)+ &
&                         k3xc(ifft,4)*((rho2r1(ifft,1)-rho2r1(ifft,2))+(3._dp/2._dp)*xccc3d2(ifft))*(rho3r1(ifft,1)-rho3r1(ifft,2))
                       end do

                     else
                       do ifft=1,nfft
!                        These sections should be rewritten, to be easier to read ... (defining intermediate scalars)
                         xc_tmp(2*ifft-1,1)= k3xc(ifft,1)*&
&                         ( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*rho3r1(2*ifft-1,2)- &
&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft,2))+   &
&                         k3xc(ifft,2)*&
&                         ( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*(rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2))- &
&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*(rho3r1(2*ifft,1)-rho3r1(2*ifft,2)))+ &
&                         k3xc(ifft,2)*&
&                         ( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*rho3r1(2*ifft-1,2)- &
&                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft,2))+ &
&                         k3xc(ifft,3)*&
&                         ( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
&                         (rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2))- &
&                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*&
&                         (rho3r1(2*ifft,1)-rho3r1(2*ifft,2)))
                         xc_tmp(2*ifft,1)=k3xc(ifft,1)*&
&                         ( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*rho3r1(2*ifft,2)+ &
&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft-1,2))+   &
&                         k3xc(ifft,2)*&
&                         ( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*(rho3r1(2*ifft,1)-rho3r1(2*ifft,2))+ &
&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*(rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2)))+ &
&                         k3xc(ifft,2)*&
&                         ( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*rho3r1(2*ifft,2)+ &
&                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft-1,2))+ &
&                         k3xc(ifft,3)*&
&                         ( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
&                         (rho3r1(2*ifft,1)-rho3r1(2*ifft,2))+ &
&                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*&
&                         (rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2)))
!                        fab: now the spin down component
                         xc_tmp(2*ifft-1,2)= k3xc(ifft,2)*&
&                         ( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*rho3r1(2*ifft-1,2)- &
&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft,2))+   &
&                         k3xc(ifft,3)*( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
&                         (rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2))- &
&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*(rho3r1(2*ifft,1)-rho3r1(2*ifft,2)))+ &
&                         k3xc(ifft,3)*( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
&                         rho3r1(2*ifft-1,2)- &
&                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft,2))+ &
&                         k3xc(ifft,4)*( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
&                         (rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2))- &
                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*&
&                         (rho3r1(2*ifft,1)-rho3r1(2*ifft,2)))
                         xc_tmp(2*ifft,2)=k3xc(ifft,1)*( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
&                         rho3r1(2*ifft,2)+ &
&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft-1,2))+   &
&                         k3xc(ifft,3)*( (rho2r1(2*ifft-1,2)+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
&                         (rho3r1(2*ifft,1)-rho3r1(2*ifft,2))+ &
&                         (rho2r1(2*ifft,2)+(3._dp/2._dp)*xccc3d2(2*ifft))*(rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2)))+ &
&                         k3xc(ifft,3)*( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
&                         rho3r1(2*ifft,2)+ &
&                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*rho3r1(2*ifft-1,2))+ &
&                         k3xc(ifft,4)*( ((rho2r1(2*ifft-1,1)-rho2r1(2*ifft-1,2))+(3._dp/2._dp)*xccc3d2(2*ifft-1))*&
&                         (rho3r1(2*ifft,1)-rho3r1(2*ifft,2))+ &
&                         ((rho2r1(2*ifft,1)-rho2r1(2*ifft,2))+(3._dp/2._dp)*xccc3d2(2*ifft))*&
&                         (rho3r1(2*ifft-1,1)-rho3r1(2*ifft-1,2)))
                       end do

!                      fab: this is the end if over cplex
                     end if
!                    fab: this is the enf if over nspden
                   end if

                   call dotprod_vn(1,rho1r1,exc3,valuei,nfft,nfftot,nspden,1,xc_tmp,ucvol,mpi_comm_sphgrid=mpi_enreg%comm_fft)
                   ABI_DEALLOCATE(xc_tmp)

!                  Perform DFPT part of the 3dte calculation

                   call timab(512,1,tsec)
                   call pead_nl_resp(cg,cg1,cg3,cplex,dtfil,dtset,d3lo,i1dir,i2dir,i3dir,i1pert,i2pert,i3pert,&
&                   kg,mband,mgfft,mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,natom,nfft,nkpt,nspden,&
&                   nspinor,nsppol,npwarr,occ,ph1d,psps,rprimd,vtrial1,xred,ylm)
                   call timab(512,2,tsec)


!                  Describe the perturbation and write out the result
                   if (mpi_enreg%me == 0) then
                     if (i2pert < natom + 1) then
                       write(message,'(a,i3,a,i3)') &
&                       ' j2 : displacement of atom ',i2pert,&
&                       ' along direction ', i2dir
                     end if
                     if (i2pert == dtset%natom + 2) then
                       write(message,'(a,i4)') &
&                       ' j2 : homogeneous electric field along direction ',&
&                       i2dir
                     end if
                     call wrtout(std_out,message,'COLL')
                     call wrtout(ab_out,message,'COLL')
                     write(ab_out,'(20x,a,13x,a)')'real part','imaginary part'
                     write(ab_out,'(5x,a2,1x,f22.10,3x,f22.10)')'xc',exc3*sixth,zero
                     if (i2pert == natom + 2) then
                       write(ab_out,'(5x,a3,f22.10,3x,f22.10)')'ddk',&
&                       d3_berry(1,i2dir),d3_berry(2,i2dir)
                     end if
                     write(ab_out,'(5x,a3,f22.10,3x,f22.10)')'dft',&
&                     d3lo(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert),&
&                     d3lo(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
                     write(ab_out,*)
                     write(std_out,'(18x,a,11x,a)')'real part','imaginary part'
                     write(std_out,'(5x,a2,1x,f20.10,3x,f20.10)')'xc',exc3*sixth,zero
                     write(std_out,'(5x,a3,f22.10,3x,f22.10)')'ddk',&
&                     d3_berry(1,i2dir),d3_berry(2,i2dir)
                     write(std_out,'(5x,a3,f22.10,3x,f22.10)')'dft',&
&                     d3lo(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert),&
&                     d3lo(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)
                     write(std_out,*)
                   end if  ! mpi_enreg%me == 0

                   d3lo(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = &
&                   d3lo(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) + exc3*sixth
                   d3lo(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = &
&                   d3lo(:,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) + d3_berry(:,i2dir)

                 end if   !rfpert
               end do    !i2dir
             end do     ! i2pert

           end if   ! rfpert
         end do    ! i3dir
       end do     ! i3pert

     end if   ! rfpert
   end do    ! i1dir
 end do     ! i1pert


 ABI_DEALLOCATE(cg1)
 ABI_DEALLOCATE(cg3)
 ABI_DEALLOCATE(eigen1)
 ABI_DEALLOCATE(rho1r1)
 ABI_DEALLOCATE(rho2r1)
 ABI_DEALLOCATE(rho2g1)
 ABI_DEALLOCATE(rho3r1)
 ABI_DEALLOCATE(atindx1)
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(nattyp)
 ABI_DEALLOCATE(ph1d)
 ABI_DEALLOCATE(ylm)
 ABI_DEALLOCATE(vtrial1)
 ABI_DEALLOCATE(vxc1)
 ABI_DEALLOCATE(vhartr1)
 ABI_DEALLOCATE(vpsp1)
 ABI_DEALLOCATE(xccc3d1)
 ABI_DEALLOCATE(xccc3d2)
 ABI_DEALLOCATE(xccc3d3)

 call timab(502,2,tsec)

end subroutine pead_nl_loop
!!***

!!****f* ABINIT/pead_nl_resp
!! NAME
!! pead_nl_resp
!!
!! FUNCTION
!! Compute the linear response part to the 3dte
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol) = array for planewave
!!                                          coefficients of wavefunctions
!!  cg1 = first-order wavefunction relative to the perturbations i1pert
!!  cg3 = first-order wavefunction relative to the perturbations i3pert
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!          if 2, COMPLEX
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  i1dir,i2dir,i3dir=directions of the corresponding perturbations
!!  i1pert,i2pert,i3pert = type of perturbation that has to be computed
!!  kg(3,mpw*mkmem)=reduced planewave coordinates
!!  mband = maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem = maximum number of k points which can fit in core memory
!!  mk1mem = maximum number of k points for first-order WF
!!           which can fit in core memory
!!  mpert =maximum number of ipert
!!  mpi_enreg=MPI-parallelisation information
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw   = maximum number of planewaves in basis sphere (large number)
!!  natom = number of atoms in unit cell
!!  nfft  = (effective) number of FFT grid points (for this processor)
!!  nkpt  = number of k points
!!  nspden = number of spin-density components
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  npwarr(nkpt) = array holding npw for each k point
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
!!  rprimd(3,3) = dimensional primitive translations (bohr)
!!  vtrial1(cplex*nfft,nspden)=firs-order local potential
!!  xred(3,natom) = reduced atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= spherical harmonics for
!!       each G and k point
!!
!! OUTPUT
!!  d3lo(2,3,mpert,3,mpert,3,mpert) = matrix of the 3DTEs
!!
!! PARENTS
!!      pead_nl_loop
!!
!! CHILDREN
!!      dotprod_g,fftpac,fourwf,init_hamiltonian
!!      mkffnl,mkkpg,nonlop,status,xmpi_sum
!!
!! SOURCE

subroutine pead_nl_resp(cg,cg1,cg3,cplex,dtfil,dtset,d3lo,&
& i1dir,i2dir,i3dir,i1pert,i2pert,i3pert,&
& kg,mband,mgfft,mkmem,mk1mem,&
& mpert,mpi_enreg,mpsang,mpw,natom,nfft,nkpt,nspden,nspinor,nsppol,&
& npwarr,occ,ph1d,psps,rprimd,vtrial1,xred,ylm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert,mband,mgfft
 integer,intent(in) :: mk1mem,mkmem,mpert,mpsang,mpw,natom,nfft,nkpt,nspden
 integer,intent(in) :: nspinor,nsppol
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),npwarr(nkpt)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: cg3(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom),rprimd(3,3)
 real(dp),intent(in) :: xred(3,natom),ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(inout) :: vtrial1(cplex*nfft,nspden)
 real(dp),intent(inout) :: d3lo(2,3,mpert,3,mpert,3,mpert)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=52
 integer :: bantot,choice,counter,cpopt,dimffnl,iband,icg0,ider,ierr
 integer :: ii,ikg,ikpt,ilm,ipw,isppol,istwf_k,jband,jj
 integer :: me,n1,n2,n3,n4,n5,n6,nband_k,nkpg,nnlout,npw_k
 integer :: option,paw_opt,signs,spaceComm,tim_fourwf,tim_nonlop
 real(dp) :: dot1i,dot1r,dot2i,dot2r,doti,dotr,lagi,lagr,sumi,sumr,weight
 type(gs_hamiltonian_type) :: gs_hamk
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp) :: buffer(2),enlout(3),kpq(3),kpt(3)
 real(dp) :: dum_svectout(1,1),dum(1),rmet(3,3),ylmgr_dum(1,1,1)
 real(dp),allocatable :: cwave0(:,:),cwavef3(:,:),ffnlk(:,:,:,:)
 real(dp),allocatable :: gh0(:,:),gh1(:,:),gvnl(:,:),kpg_k(:,:)
 real(dp),allocatable :: vlocal1(:,:,:),wfraug(:,:,:,:),ylm_k(:,:)
 type(pawcprj_type) :: cprj_dum(1,1)
 type(pawtab_type) :: pawtab_dum(0)

!***********************************************************************

 ABI_UNUSED(dtfil%ireadwf)

 me = mpi_enreg%me
 spaceComm=mpi_enreg%comm_cell

 bantot = 0
 icg0 = 0
 option = 2
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)

 ABI_ALLOCATE(vlocal1,(cplex*n4,n5,n6))
 ABI_ALLOCATE(wfraug,(2,n4,n5,n6))

!Initialize Hamiltonian (k-independent terms) - NCPP only
 call init_hamiltonian(gs_hamk,psps,pawtab_dum,nspinor,nsppol,nspden,natom,&
& dtset%typat,xred,nfft,mgfft,dtset%ngfft,rprimd,dtset%nloalg,ph1d=ph1d,&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab)
!& paw_ij=paw_ij)
 rmet = MATMUL(TRANSPOSE(rprimd),rprimd)

 sumr = zero ; sumi = zero

!Loop over spins

 do isppol = 1, nsppol

   call fftpac(isppol,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,vtrial1,vlocal1,option)

!  Loop over k-points

   ikg = 0
   do ikpt = 1, nkpt

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,-1,mpi_enreg%me))cycle

     counter = 100*ikpt

     nband_k = dtset%nband(ikpt+(isppol-1)*nkpt)
     npw_k = npwarr(ikpt)
     istwf_k = dtset%istwfk(ikpt)

     kpt(:) = dtset%kptns(:,ikpt)
     kpq(:) = dtset%kptns(:,ikpt) ! In case of non zero q, kpt = kpt + q

     ABI_ALLOCATE(cwave0,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(cwavef3,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(gh0,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(gvnl,(2,npw_k*dtset%nspinor))
     ABI_ALLOCATE(gh1,(2,npw_k*dtset%nspinor))

     ABI_ALLOCATE(kg_k,(3,npw_k))
     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     kg_k(:,1:npw_k) = kg(:,1+ikg:npw_k+ikg)
     if (psps%useylm==1) then
       do ilm=1,mpsang*mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
     end if

!    Compute (k+G) and (k+q+G) vectors (only if useylm=1)
     nkpg=0;if (i2pert<natom+1) nkpg=3*dtset%nloalg(3)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpt,nkpg,npw_k)
     end if

!    Compute nonlocal form factors ffnl at (k+G), for all atoms
     dimffnl=1
     ABI_ALLOCATE(ffnlk,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
     if (i2pert<natom+1) then
       ider=0
       call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnlk,psps%ffspl,gs_hamk%gmet,gs_hamk%gprimd,&
&       ider,ider,psps%indlmn,kg_k,kpg_k,kpt,psps%lmnmax,psps%lnmax,psps%mpsang,&
&       psps%mqgrid_ff,nkpg,npw_k,psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&       psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
     end if

!    Load k-dependent part in the Hamiltonian datastructure
     call gs_hamk%load_k(kpt_k=kpt,npw_k=npw_k,istwf_k=istwf_k,&
&     kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnlk,compute_gbound=.true.)
!    Load k+q-dependent part in the Hamiltonian datastructure
!    call load_kprime_hamiltonian...  !! To be activated when q/=0

!    Loop over bands

     do iband = 1,nband_k

       cwave0(:,:)=cg(:,1+(iband - 1)*npw_k*dtset%nspinor+icg0:&
&       iband*npw_k*dtset%nspinor+icg0)
       cwavef3(:,:)=cg3(:,1+(iband-1)*npw_k*dtset%nspinor+icg0:&
&       iband*npw_k*dtset%nspinor+icg0)

!      Compute vtrial1 | cwafef3 >
       tim_fourwf = 0 ; weight = one
       call fourwf(cplex,vlocal1,cwavef3,gh1,wfraug,gs_hamk%gbound_k,gs_hamk%gbound_k,&
&       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,dtset%ngfft,npw_k,npw_k,n4,n5,n6,option,&
&       tim_fourwf,weight,weight,&
&       use_gpu_cuda=dtset%use_gpu_cuda)

!      In case i2pert = phonon-type perturbation
!      add first-order change in the nonlocal potential
       if (i2pert<natom+1) then
         signs=2 ; choice=2 ; nnlout=3 ; tim_nonlop = 0 ; paw_opt=0 ; cpopt=-1
         call nonlop(choice,cpopt,cprj_dum,enlout,gs_hamk,i2dir,dum,mpi_enreg,1,nnlout,paw_opt,&
&         signs,dum_svectout,tim_nonlop,cwavef3,gvnl,iatom_only=i2pert)
         gh1(:,:) = gh1(:,:) + gvnl(:,:)
       end if

       ii = (iband-1)*npw_k*dtset%nspinor + icg0
       call dotprod_g(dotr,doti,istwf_k,npw_k,2,cg1(:,ii+1:ii+npw_k),gh1,mpi_enreg%me_g0,xmpi_comm_self)

!      Compute vtrial1 | cwave0 >
       tim_fourwf = 0 ; weight = one
       call fourwf(cplex,vlocal1,cwave0,gh0,wfraug,gs_hamk%gbound_k,gs_hamk%gbound_k,&
&       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,dtset%ngfft,npw_k,npw_k,n4,n5,n6,option,&
&       tim_fourwf,weight,weight,use_gpu_cuda=dtset%use_gpu_cuda)

!      In case i2pert = phonon-type perturbation
!      add first-order change in the nonlocal potential
       if (i2pert<natom+1) then
         signs=2 ; choice=2 ; nnlout=3 ; tim_nonlop = 0 ; paw_opt=0 ; cpopt=-1
         call nonlop(choice,cpopt,cprj_dum,enlout,gs_hamk,i2dir,dum,mpi_enreg,1,nnlout,paw_opt,&
&         signs,dum_svectout,tim_nonlop,cwave0,gvnl,iatom_only=i2pert)
         gh0(:,:) = gh0(:,:) + gvnl(:,:)
       end if

!      Compute the dft contribution to the Lagrange multiplier
!      cwavef3 and cwave0 have been transferred to gh1 and gh0
!      these vectors will be used to store the wavefunctions of band iband
!      cg1 and gh0 contain the wavefunctions of band jband

       lagr = zero ; lagi = zero
       do jband = 1, nband_k

         ii = (jband - 1)*npw_k*dtset%nspinor + icg0
         jj = (iband - 1)*npw_k*dtset%nspinor + icg0

!        dot1r and dot1i contain < u_mk | v^(1) | u_nk >
!        dot2r and dot2i contain < u_nk^(1) | u_mk^(1) >
!        m -> jband and n -> iband

         dot1r = zero ; dot1i = zero
         dot2r = zero ; dot2i = zero
         do ipw = 1, npw_k
           ii = ii + 1 ; jj = jj + 1
           dot1r = dot1r + cg(1,ii)*gh0(1,ipw) + cg(2,ii)*gh0(2,ipw)
           dot1i = dot1i + cg(1,ii)*gh0(2,ipw) - cg(2,ii)*gh0(1,ipw)
           dot2r = dot2r + cg1(1,jj)*cg3(1,ii) + &
&           cg1(2,jj)*cg3(2,ii)
           dot2i = dot2i + cg1(1,jj)*cg3(2,ii) - &
&           cg1(2,jj)*cg3(1,ii)
         end do  !  ipw

         lagr = lagr + dot1r*dot2r - dot1i*dot2i
         lagi = lagi + dot1r*dot2i + dot1i*dot2r

       end do    ! jband

       sumr = sumr + &
&       dtset%wtk(ikpt)*occ(bantot+iband)*(dotr-lagr)
       sumi = sumi + &
&       dtset%wtk(ikpt)*occ(bantot+iband)*(doti-lagi)

     end do   ! end loop over bands

     bantot = bantot + nband_k
     icg0 = icg0 + npw_k*dtset%nspinor*nband_k
     ikg = ikg + npw_k

     ABI_DEALLOCATE(cwave0)
     ABI_DEALLOCATE(cwavef3)
     ABI_DEALLOCATE(gh0)
     ABI_DEALLOCATE(gh1)
     ABI_DEALLOCATE(gvnl)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ffnlk)
     ABI_DEALLOCATE(kpg_k)

   end do   ! end loop over k-points

 end do   ! end loop over spins

 if (xmpi_paral == 1) then
   buffer(1) = sumr ; buffer(2) = sumi
   call xmpi_sum(buffer,spaceComm,ierr)
   sumr = buffer(1) ; sumi = buffer(2)
 end if


 d3lo(1,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sumr
!d3lo(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = sumi

!In some cases, the imaginary part is /= 0 because of the
!use of time reversal symmetry
 d3lo(2,i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = zero

 call gs_hamk%free()

 ABI_DEALLOCATE(vlocal1)
 ABI_DEALLOCATE(wfraug)

end subroutine pead_nl_resp
!!***

!!****f* ABINIT/pead_nl_mv
!! NAME
!! pead_nl_mv
!!
!! FUNCTION
!! Compute the finite difference expression of the k-point derivative
!! using the PEAD formulation of the third-order energy
!! (see Nunes and Gonze PRB 63, 155107 (2001) [[cite:Nunes2001]] Eq. 102)
!! and the finite difference formula of Marzari and Vanderbilt
!! (see Marzari and Vanderbilt, PRB 56, 12847 (1997) [[cite:Marzari1997]], Appendix B)
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol) = array for planewave coefficients of wavefunctions
!!  cgindex(nkpt2,nsppol) = for each k-point, cgindex stores the location of the WF in the cg array
!!  cg1 = first-order wavefunction relative to the perturbations i1pert
!!  cg3 = first-order wavefunction relative to the perturbations i3pert
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  i1pert,i3pert = type of perturbation that has to be computed
!!  i1dir,i3dir=directions of the corresponding perturbations
!!  kneigh(30,nkpt2) = index of the neighbours of each k-point
!!  kg_neigh(30,nkpt2,3) = necessary to construct the vector joining a k-point
!!                         to its nearest neighbour in case of a single k-point,
!!                         a line of k-points or a plane of k-points.
!!                         See getshell.F90 for details
!!  kptindex(2,nkpt3)= index of the k-points in the reduced BZ
!!                     related to a k-point in the full BZ
!!  kpt3(3,nkpt3) = reduced coordinates of k-points in the full BZ
!!  mband = maximum number of bands
!!  mkmem = maximum number of k points which can fit in core memory
!!  mkmem_max = maximal number of k-points on each processor (MPI //)
!!  mk1mem = maximum number of k points for first-order WF which can fit in core memory
!!  mpi_enreg=MPI-parallelisation information
!!  mpw   = maximum number of planewaves in basis sphere (large number)
!!  mvwtk(30,nkpt) = weights to compute the finite difference ddk
!!  natom = number of atoms in unit cell
!!  nkpt2 = number of k-points in the reduced part of the BZ
!!          nkpt2 = nkpt/2 in case of time-reversal symmetry (kptopt = 2)
!!  nkpt3 = number of k-points in the full BZ
!!  nneigh = total number of neighbours required to evaluate the finite difference formula
!!  npwarr(nkpt) = array holding npw for each k point
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = number of channels for spin-polarization (1 or 2)
!!  pwind(mpw,nneigh,mkmem) = array used to compute the overlap matrix smat between k-points
!!
!! OUTPUT
!!  d3_berry(2,3) = Berry-phase part of the third-order energy
!!
!! SIDE EFFECTS
!!  mpi_enreg=MPI-parallelisation information
!!
!! NOTES
!! For a given set of values of i1pert,i3pert,i1dir and
!! i3dir, the routine computes the k-point derivatives for
!! 12dir = 1,2,3
!!
!! PARENTS
!!      pead_nl_loop
!!
!! CHILDREN
!!      dzgedi,dzgefa,mpi_recv,mpi_send,status,wrtout,xmpi_sum
!!
!! SOURCE

subroutine pead_nl_mv(cg,cgindex,cg1,cg3,dtset,dtfil,d3_berry,gmet,&
&                   i1pert,i3pert,i1dir,i3dir,kneigh,kg_neigh,kptindex,&
&                   kpt3,mband,mkmem,mkmem_max,mk1mem,mpi_enreg,&
&                   mpw,mvwtk,natom,nkpt2,nkpt3,nneigh,npwarr,nspinor,&
&                   nsppol,pwind)

 use m_hide_lapack, only : dzgedi, dzgefa

!Arguments ------------------------------------
!
!---  Arguments : integer scalars
 integer, intent(in) :: i1dir,i1pert,i3dir,i3pert,mband,mk1mem
 integer, intent(in) :: mkmem,mkmem_max,mpw,natom
 integer, intent(in) :: nkpt2,nkpt3,nneigh,nspinor,nsppol
!
!---  Arguments : integer arrays
 integer, intent(in) :: cgindex(nkpt2,nsppol)
 integer, intent(in) :: kneigh(30,nkpt2),kg_neigh(30,nkpt2,3),kptindex(2,nkpt3)
 integer, intent(in) :: npwarr(nkpt2),pwind(mpw,nneigh,mkmem)
!
!---  Arguments : real(dp) scalars
!
!---  Arguments : real(dp) arrays
 real(dp), intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp), intent(in) :: cg1(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp), intent(in) :: cg3(2,mpw*nspinor*mband*mk1mem*nsppol)
 real(dp), intent(in) :: gmet(3,3),kpt3(3,nkpt3)
 real(dp), intent(in) :: mvwtk(30,nkpt2)
 real(dp), intent(out) :: d3_berry(2,3)
!
!---  Arguments : structured datatypes
 type(MPI_type), intent(in) :: mpi_enreg
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(in) :: dtset

!Local variables-------------------------------
!
!---- Local variables : integer scalars
 integer :: count,counter,count1,iband,icg
 integer :: ierr,ii,ikpt,ikpt_loc,ikpt2
 integer :: ikpt_rbz,ineigh,info,ipw,isppol,jband,jcg,jj,jkpt,job,jpw, jkpt2, jkpt_rbz
 integer :: lband,lpband,nband_occ,npw_k,npw_k1,my_source,his_source,dest,tag
 integer :: spaceComm
 integer,parameter :: level=52
 integer :: bdtot_index
!
!---- Local variables : integer arrays
 integer,allocatable :: ipvt(:)
 integer, allocatable :: bd_index(:,:)
!
!---- Local variables : real(dp) scalars
 real(dp) :: dotnegi,dotnegr,dotposi,dotposr
! real(dp) :: c1,c2 ! appear commented out below
!
!---- Local variables : real(dp) arrays
 real(dp) :: d3_aux(2,3),det(2,2),dk(3),dk_(3)
 real(dp) :: z1(2),z2(2)
 real(dp),allocatable :: buffer(:,:),cgq(:,:),cg1q(:,:),cg3q(:,:)
 real(dp),allocatable :: qmat(:,:,:),s13mat(:,:,:),s1mat(:,:,:),s3mat(:,:,:)
 real(dp),allocatable :: smat(:,:,:),zgwork(:,:)
!
!---- Local variables : character variables
 character(len=500) :: message
!
!---- Local variables : structured datatypes


#if defined HAVE_MPI
integer :: status1(MPI_STATUS_SIZE)
spaceComm=mpi_enreg%comm_cell
#endif

 ABI_UNUSED(dtfil%ireadwf)

! ***********************************************************************

 write(message,'(8a)') ch10,&
& ' pead_nl_mv : finite difference expression of the k-point derivative',ch10,&
& '           is performed using the PEAD formulation of ',&
& 'the third-order energy',ch10,&
& '           (see Nunes and Gonze PRB 63, 155107 (2001) [[cite:Nunes2001]] Eq. 102)',ch10
!call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')


!fab: I think that the following restriction must be eliminated:
!isppol = 1

 ikpt_loc = 0
 d3_aux(:,:) = 0_dp

 ABI_ALLOCATE(s13mat,(2,mband,mband))
 ABI_ALLOCATE(smat,(2,mband,mband))
 ABI_ALLOCATE(s1mat,(2,mband,mband))
 ABI_ALLOCATE(qmat,(2,mband,mband))
 ABI_ALLOCATE(ipvt,(mband))
 ABI_ALLOCATE(s3mat,(2,mband,mband))
 ABI_ALLOCATE(zgwork,(2,mband))
 ABI_ALLOCATE(bd_index, (nkpt2, nsppol))

 bdtot_index = 0
 do isppol = 1, nsppol
   do ikpt_rbz = 1, nkpt2
     bd_index(ikpt_rbz,isppol) = bdtot_index
     bdtot_index = bdtot_index + dtset%nband(ikpt_rbz+nkpt2*(isppol-1))
   end do
 end do

!fab: I think here I have to add the loop over spin

 do isppol = 1, nsppol

!  Loop over k-points
!  COMMENT: Every processor has to make mkmem_max iterations
!  even if mkmem < mkemem_max. This is due to the fact
!  that it still has to communicate its wavefunctions
!  to other processors even if it has no more overlap
!  matrices to compute.

   ikpt_loc = 0 ; ikpt = 0

   do while (ikpt_loc < mkmem_max)

     if (ikpt_loc < mkmem) ikpt = ikpt + 1

     if (xmpi_paral == 1) then
!      if ((minval(abs(mpi_enreg%proc_distrb(ikpt,1:mband,1:dtset%nsppol) &
!      &       - mpi_enreg%me)) /= 0).and.(ikpt_loc < mkmem)) cycle
       if(ikpt>nkpt2)then
         ikpt_loc=mkmem_max
         cycle
       end if
       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,mband,-1,mpi_enreg%me)) then
         if(ikpt==nkpt2) ikpt_loc=mkmem_max
         cycle
       end if
     end if

     ikpt_loc = ikpt_loc + 1
     npw_k = npwarr(ikpt)
     counter = 100*ikpt

     ii = cgindex(ikpt,isppol)

!    Loop on the  neighbours

     do ineigh = 1,nneigh

       s13mat(:,:,:) = zero
       smat(:,:,:) = zero
       s1mat(:,:,:) = zero
       s3mat(:,:,:) = zero
       qmat(:,:,:) = zero

       ikpt2  = kneigh(ineigh,ikpt)
       ikpt_rbz = kptindex(1,ikpt2)   ! index of the k-point in the reduced BZ
       jj = cgindex(ikpt_rbz,isppol)
       ! previous fixed value for nband_k now called nband_occ:
       !nband_occ = dtset%nband(ikpt_rbz+nkpt2*(isppol-1))
       ! TODO: check if all these bands are occupied in nsppol = 2 case
       nband_occ = 0
       do iband = 1, dtset%nband(ikpt_rbz+nkpt2*(isppol-1))
         !Note, only one image is allowed here (or occ_orig should be the same or all images)
         if (dtset%occ_orig(bd_index(ikpt_rbz,isppol) + iband,1) > tol10) nband_occ = nband_occ + 1
       end do
       npw_k1 = npwarr(ikpt_rbz)
       dk_(:) = kpt3(:,ikpt2) - dtset%kptns(:,ikpt)
       dk(:)  = dk_(:) - nint(dk_(:)) + real(kg_neigh(ineigh,ikpt,:),dp)

       count = nspinor*mband*npw_k1
       ABI_ALLOCATE(cgq,(2,count))
       ABI_ALLOCATE(cg1q,(2,count))
       ABI_ALLOCATE(cg3q,(2,count))

#if defined HAVE_MPI

       my_source = mpi_enreg%proc_distrb(ikpt_rbz,1,1)

!      do dest = 0, mpi_enreg%nproc-1
       do dest = 0, maxval(mpi_enreg%proc_distrb(1:nkpt2,1:mband,1:dtset%nsppol))

         if ((dest==mpi_enreg%me).and.(ikpt_loc <= mkmem)) then
!          I am dest and have something to do

           if (my_source == mpi_enreg%me) then
!            I am destination and source
             jcg = cgindex(ikpt_rbz,isppol)

             cgq(:,1:count)  = cg(:,jcg+1:jcg+count)
             cg1q(:,1:count) = cg1(:,jcg+1:jcg+count)
             cg3q(:,1:count) = cg3(:,jcg+1:jcg+count)

           else
!            I am the destination but not the source -> receive

             tag = ikpt_rbz

             ABI_ALLOCATE(buffer,(2,3*count))

             call MPI_RECV(buffer,2*3*count,MPI_DOUBLE_PRECISION,my_source,tag,spaceComm,status1,ierr)

             cgq(:,1:count)  = buffer(:,1:count)
             cg1q(:,1:count) = buffer(:,count+1:2*count)
             cg3q(:,1:count) = buffer(:,2*count+1:3*count)
             ABI_DEALLOCATE(buffer)

           end if

         else if (ikpt_loc <= mpi_enreg%mkmem(dest)) then  ! dest != me and the dest has a k-point to treat

           jkpt=mpi_enreg%kpt_loc2ibz_sp(dest, ikpt_loc,1)
           jkpt2  = kneigh(ineigh,jkpt)
           jkpt_rbz = kptindex(1,jkpt2)   ! index of the k-point in the reduced BZ

           his_source = mpi_enreg%proc_distrb(jkpt_rbz,1,1)

           if (his_source == mpi_enreg%me) then

             jcg = cgindex(jkpt_rbz,isppol)

             tag = jkpt_rbz
             count1 = npwarr(jkpt_rbz)*mband*nspinor
             ABI_ALLOCATE(buffer,(2,3*count1))
             buffer(:,1:count1)            = cg(:,jcg+1:jcg+count1)
             buffer(:,count1+1:2*count1)   = cg1(:,jcg+1:jcg+count1)
             buffer(:,2*count1+1:3*count1) = cg3(:,jcg+1:jcg+count1)

             call MPI_SEND(buffer,2*3*count1,MPI_DOUBLE_PRECISION,dest,tag,spaceComm,ierr)

             ABI_DEALLOCATE(buffer)

           end if

         end if

       end do
!
!      do jkpt = 1, nkpt2
!
!      if ((jkpt == ikpt_rbz).and.(source /= mpi_enreg%me).and.&
!      &         (ikpt_loc <= mkmem)) then
!
!      tag = jkpt
!
!      allocate(buffer(2,3*count))
!      call MPI_RECV(buffer,2*3*count,MPI_DOUBLE_PRECISION,&
!      source,tag,spaceComm,status1,ierr)
!
!      cgq(:,1:count)  = buffer(:,1:count)
!      cg1q(:,1:count) = buffer(:,count+1:2*count)
!      cg3q(:,1:count) = buffer(:,2*count+1:3*count)
!      deallocate(buffer)
!
!      end if
!
!      !        ----------------------------------------------------------------------------
!      !        --------------- Here: send the WF to all the cpus that need it -------------
!      !        ----------------------------------------------------------------------------
!
!      do dest = 1, mpi_enreg%nproc
!
!      if ((minval(abs(mpi_enreg%proc_distrb(jkpt,1:mband,1:dtset%nsppol) &
!      &           - mpi_enreg%me)) == 0).and.&
!      &           (mpi_enreg%kptdstrb(dest,ineigh,ikpt_loc) == jkpt)) then
!
!
!
!      jcg = cgindex(jkpt,isppol)
!
!      if (((dest-1) == mpi_enreg%me)) then
!
!      cgq(:,1:count)  = cg(:,jcg+1:jcg+count)
!      cg1q(:,1:count) = cg1(:,jcg+1:jcg+count)
!      cg3q(:,1:count) = cg3(:,jcg+1:jcg+count)
!
!      else
!
!      tag = jkpt
!      count1 = npwarr(jkpt)*mband*nspinor
!      allocate(buffer(2,3*count1))
!      buffer(:,1:count1)            = cg(:,jcg+1:jcg+count1)
!      buffer(:,count1+1:2*count1)   = cg1(:,jcg+1:jcg+count1)
!      buffer(:,2*count1+1:3*count1) = cg3(:,jcg+1:jcg+count1)
!
!      call MPI_SEND(buffer,2*3*count1,MPI_DOUBLE_PRECISION,(dest-1),tag,spaceComm,ierr)
!
!      deallocate(buffer)
!
!      end if
!
!      end if
!
!      end do          ! loop over dest
!
!      end do          ! loop over jkpt

       if (ikpt_loc > mkmem) then
         ABI_DEALLOCATE(cgq)
         ABI_DEALLOCATE(cg1q)
         ABI_DEALLOCATE(cg3q)
         cycle
       end if

#else
!      no // over k-points

       cgq(:,1:count)  = cg(:,jj+1:jj+count)
       cg1q(:,1:count) = cg1(:,jj+1:jj+count)
       cg3q(:,1:count) = cg3(:,jj+1:jj+count)

#endif

!      Compute overlap matrices

       if (kptindex(2,ikpt2) == 0) then  ! no time-reversal symmetry

         do ipw = 1, npw_k

           jpw = pwind(ipw,ineigh,ikpt_loc)
           if (jpw /= 0) then

             do iband = 1, nband_occ
               do jband = 1, nband_occ

                 icg = ii + (iband-1)*npw_k + ipw
                 jcg = (jband-1)*npw_k1 + jpw

                 smat(1,iband,jband) = smat(1,iband,jband) + &
&                 cg(1,icg)*cgq(1,jcg) + cg(2,icg)*cgq(2,jcg)
                 smat(2,iband,jband) = smat(2,iband,jband) + &
&                 cg(1,icg)*cgq(2,jcg) - cg(2,icg)*cgq(1,jcg)

                 s13mat(1,iband,jband) = s13mat(1,iband,jband) + &
&                 cg1(1,icg)*cg3q(1,jcg) + cg1(2,icg)*cg3q(2,jcg)
                 s13mat(2,iband,jband) = s13mat(2,iband,jband) + &
&                 cg1(1,icg)*cg3q(2,jcg) - cg1(2,icg)*cg3q(1,jcg)

                 s1mat(1,iband,jband) = s1mat(1,iband,jband) + &
&                 cg1(1,icg)*cgq(1,jcg) + cg1(2,icg)*cgq(2,jcg) + &
&                 cg(1,icg)*cg1q(1,jcg) + cg(2,icg)*cg1q(2,jcg)
                 s1mat(2,iband,jband) = s1mat(2,iband,jband) + &
&                 cg1(1,icg)*cgq(2,jcg) - cg1(2,icg)*cgq(1,jcg) + &
&                 cg(1,icg)*cg1q(2,jcg) - cg(2,icg)*cg1q(1,jcg)

                 s3mat(1,iband,jband) = s3mat(1,iband,jband) + &
&                 cg3(1,icg)*cgq(1,jcg) + cg3(2,icg)*cgq(2,jcg) + &
&                 cg(1,icg)*cg3q(1,jcg) + cg(2,icg)*cg3q(2,jcg)
                 s3mat(2,iband,jband) = s3mat(2,iband,jband) + &
&                 cg3(1,icg)*cgq(2,jcg) - cg3(2,icg)*cgq(1,jcg) + &
&                 cg(1,icg)*cg3q(2,jcg) - cg(2,icg)*cg3q(1,jcg)

               end do
             end do

           end if

         end do   ! ipw

       else                              ! use time-reversal symmetry

         do ipw = 1,npw_k

           jpw = pwind(ipw,ineigh,ikpt_loc)
           if (jpw /= 0) then

             do iband = 1, nband_occ
               do jband = 1, nband_occ

                 icg = ii + (iband-1)*npw_k + ipw
                 jcg = (jband-1)*npw_k1 + jpw

                 smat(1,iband,jband) = smat(1,iband,jband) + &
&                 cg(1,icg)*cgq(1,jcg) - cg(2,icg)*cgq(2,jcg)
                 smat(2,iband,jband) = smat(2,iband,jband) - &
&                 cg(1,icg)*cgq(2,jcg) - cg(2,icg)*cgq(1,jcg)

                 s13mat(1,iband,jband) = s13mat(1,iband,jband) + &
&                 cg1(1,icg)*cg3q(1,jcg) - cg1(2,icg)*cg3q(2,jcg)
                 s13mat(2,iband,jband) = s13mat(2,iband,jband) - &
&                 cg1(1,icg)*cg3q(2,jcg) - cg1(2,icg)*cg3q(1,jcg)

                 s1mat(1,iband,jband) = s1mat(1,iband,jband) + &
&                 cg1(1,icg)*cgq(1,jcg) - cg1(2,icg)*cgq(2,jcg) + &
&                 cg(1,icg)*cg1q(1,jcg) - cg(2,icg)*cg1q(2,jcg)
                 s1mat(2,iband,jband) = s1mat(2,iband,jband) - &
&                 cg1(1,icg)*cgq(2,jcg) - cg1(2,icg)*cgq(1,jcg) - &
&                 cg(1,icg)*cg1q(2,jcg) - cg(2,icg)*cg1q(1,jcg)

                 s3mat(1,iband,jband) = s3mat(1,iband,jband) + &
&                 cg3(1,icg)*cgq(1,jcg) - cg3(2,icg)*cgq(2,jcg) + &
&                 cg(1,icg)*cg3q(1,jcg) - cg(2,icg)*cg3q(2,jcg)
                 s3mat(2,iband,jband) = s3mat(2,iband,jband) - &
&                 cg3(1,icg)*cgq(2,jcg) - cg3(2,icg)*cgq(1,jcg) - &
&                 cg(1,icg)*cg3q(2,jcg) - cg(2,icg)*cg3q(1,jcg)

               end do
             end do

           end if

         end do   ! ipw

       end if

       ABI_DEALLOCATE(cgq)
       ABI_DEALLOCATE(cg1q)
       ABI_DEALLOCATE(cg3q)

!      Compute qmat, the inverse of smat

       job = 1  ! compute inverse only
       qmat(:,:,:) = smat(:,:,:)

       call dzgefa(qmat,mband,nband_occ,ipvt,info)
       call dzgedi(qmat,mband,nband_occ,ipvt,det,zgwork,job)

!      DEBUG
!      write(100,*)
!      write(100,*)'ikpt = ',ikpt,'ineigh = ',ineigh
!      do iband = 1,nband_occ
!      do jband = 1,nband_occ
!      c1 = 0_dp ; c2 = 0_dp
!      do lband = 1,nband_occ
!      c1 = c1 + smat(1,iband,lband)*qmat(1,lband,jband) - &
!      &           smat(2,iband,lband)*qmat(2,lband,jband)
!      c2 = c2 + smat(1,iband,lband)*qmat(2,lband,jband) + &
!      &           smat(2,iband,lband)*qmat(1,lband,jband)
!      end do
!      write(100,'(2(2x,i2),2(2x,f16.9))')iband,jband,&
!      & c1,c2
!      end do
!      end do
!      ENDDEBUG



!      Accumulate sum over bands

       dotposr = 0_dp ; dotposi = 0_dp
       dotnegr = 0_dp ; dotnegi = 0_dp
       do iband = 1, nband_occ
         do jband = 1, nband_occ

           dotposr = dotposr + &
&           s13mat(1,iband,jband)*qmat(1,jband,iband) - &
&           s13mat(2,iband,jband)*qmat(2,jband,iband)
           dotposi = dotposi + &
&           s13mat(1,iband,jband)*qmat(2,jband,iband) + &
&           s13mat(2,iband,jband)*qmat(1,jband,iband)


           do lband = 1, nband_occ
             do lpband= 1, nband_occ

               z1(1) = s1mat(1,iband,jband)*qmat(1,jband,lband) - &
&               s1mat(2,iband,jband)*qmat(2,jband,lband)
               z1(2) = s1mat(1,iband,jband)*qmat(2,jband,lband) + &
&               s1mat(2,iband,jband)*qmat(1,jband,lband)

               z2(1) = s3mat(1,lband,lpband)*qmat(1,lpband,iband) - &
&               s3mat(2,lband,lpband)*qmat(2,lpband,iband)
               z2(2) = s3mat(1,lband,lpband)*qmat(2,lpband,iband) + &
&               s3mat(2,lband,lpband)*qmat(1,lpband,iband)

               dotnegr = dotnegr + &
&               z1(1)*z2(1) - z1(2)*z2(2)
               dotnegi = dotnegi + &
&               z1(1)*z2(2) + z1(2)*z2(1)

             end do   ! lpband
           end do   ! lband

         end do   ! jband
       end do   ! iband

       d3_aux(1,:) = d3_aux(1,:) + &
&       dk(:)*mvwtk(ineigh,ikpt)*dtset%wtk(ikpt)*(2_dp*dotposr-dotnegr)
       d3_aux(2,:) = d3_aux(2,:) + &
&       dk(:)*mvwtk(ineigh,ikpt)*dtset%wtk(ikpt)*(2_dp*dotposi-dotnegi)

     end do        ! End loop over neighbours


   end do      ! End loop over k-points

 end do  ! fab: end loop over spin




 call xmpi_sum(d3_aux,spaceComm,ierr)


 ABI_DEALLOCATE(s13mat)
 ABI_DEALLOCATE(smat)
 ABI_DEALLOCATE(s1mat)
 ABI_DEALLOCATE(qmat)
 ABI_DEALLOCATE(ipvt)
 ABI_DEALLOCATE(s3mat)
 ABI_DEALLOCATE(zgwork)
 ABI_DEALLOCATE(bd_index)


!fab: I think that in the following we have to make a distinction:
!for the spin unpolarized case we leave the PEAD expression as it is, while
!in the spin polarized case we have simply to divide by 2
!(see eq.19 di PRB 63,155107 [[cite:Nunes2001]], eq. 7 di PRB 71,125107 [[cite:Veithen2005]]
! and eq 13 di PRB 71, 125107 [[cite:Veithen2005]] ...
!in this latter equation the 2 must be simply replaced by the sum over the spin components...
!and indeed we have inserted the loop over the spin,
!but there was a factor 2 already present in the routine due to spin degenracy that had to be removed)


 if (nsppol==1) then

!  Take minus the imaginary part

   d3_berry(1,:) = -1_dp*d3_aux(2,:)
   d3_berry(2,:) = d3_aux(1,:)

   d3_berry(2,:) = 0_dp

 else

   d3_berry(1,:) = -1_dp*d3_aux(2,:)/2._dp
   d3_berry(2,:) = d3_aux(1,:)/2._dp

   d3_berry(2,:) = 0_dp/2._dp

 end if

!DEBUG
!write(100,*)'pead_nl_mv.f : d3_berry'
!write(100,*)'Perturbation',i1dir,i3dir
!write(100,*)
!write(100,*)'before transformation'
!write(100,*)'real part'
!write(100,'(3(2x,f20.9))')d3_berry(1,:)
!write(100,*)
!write(100,*)'imaginary part'
!write(100,'(3(2x,f20.9))')d3_berry(2,:)
!write(100,*)
!write(100,*)'after transformation'
!ENDDEBUG

!Compute the projection on the basis vectors of
!reciprocal space

 d3_aux(:,:) = 0_dp
 do ii = 1,3
   do jj = 1,3
     d3_aux(:,ii) = d3_aux(:,ii) + gmet(ii,jj)*d3_berry(:,jj)
   end do
 end do
 d3_berry(:,:) = d3_aux(:,:)

!Write out the berryphase part of the third order energy

 if (mpi_enreg%me == 0) then

   write(message,'(a,a,a)')ch10,' Berryphase part of the third-order energy:',ch10
   call wrtout(std_out,  message,'COLL')

   if (i1pert < natom + 1) then
     write(message,'(a,i3,a,i3)')&
&     '            j1: Displacement of atom ',i1pert,&
&     ' along direction ',i1dir
   else if (i1pert == natom + 2) then
     write(message,'(a,i3)')&
&     '            j1: homogenous electric field along direction ',i1dir
   end if
   call wrtout(std_out,  message,'COLL')

   write(message,'(a)')&
&   '            j2: k-point derivative along direction i2dir '
   call wrtout(std_out,  message,'COLL')

   if (i3pert < natom + 1) then
     write(message,'(a,i3,a,i3,a)')&
&     '            j3: Displacement of atom ',i3pert,&
&     ' along direction ',i3dir,ch10
   else if (i3pert == natom + 2) then
     write(message,'(a,i3,a)')&
&     '            j3: homogenous electric field along direction ',i3dir,ch10
   end if
   call wrtout(std_out,  message,'COLL')

!  write(ab_out,'(5x,a5,8x,a9,5x,a14)')'i2dir','real part','imaginary part'
   write(std_out,'(5x,a5,8x,a9,5x,a14)')'i2dir','real part','imaginary part'
   do ii = 1,3
     write(std_out,'(7x,i1,3x,f16.9,3x,f16.9)')ii,&
&     d3_berry(1,ii),d3_berry(2,ii)
     write(std_out,'(7x,i1,3x,f16.9,3x,f16.9)')ii,&
&     d3_berry(1,ii),d3_berry(2,ii)
   end do

 end if    ! mpi_enreg%me == 0

!DEBUG
!write(100,*)'real part'
!write(100,'(3(2x,f20.9))')d3_berry(1,:)
!write(100,*)
!write(100,*)'imaginary part'
!write(100,'(3(2x,f20.9))')d3_berry(2,:)
!ENDDEBUG

end subroutine pead_nl_mv
!!***

end module m_pead_nl_loop
!!***
