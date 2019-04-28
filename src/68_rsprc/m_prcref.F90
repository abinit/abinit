!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_prcref
!! NAME
!!  m_prcref
!!
!! FUNCTION
!!  Routines to precondition residual potential (or density) and forces.
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, MT, PMA)
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

module m_prcref

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_errors
 use m_abicore
 use m_xmpi
 use m_xcdata

 use m_frskerker1
 use m_frskerker2
 use mod_prc_memory

 use m_time,     only : timab
 use m_numeric_tools, only : dotproduct
 use m_geometry, only : xcart2xred, metric
 use m_cgtools,  only : dotprod_vn, mean_fftr
 use m_mpinfo,   only : ptabs_fourdp, destroy_mpi_enreg, initmpi_seq
 use m_pawtab,   only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type
 use m_fftcore,  only : kgindex
 use m_fft,      only : zerosym, indirect_parallel_fourier, fourdp
 use m_kg,       only : getph
 use m_dtset,    only : testsusmat
 use m_spacepar, only : hartre, laplacian
 use m_distribfft, only : init_distribfft_seq
 use m_forces,     only : fresid
 use m_atm2fft,    only : atm2fft
 use m_rhotoxc,    only : rhotoxc
 use m_mklocl,     only : mklocl
 use m_mkcore,     only : mkcore

 implicit none

 private
!!***

 public :: prcref
 public :: prcref_PMA
 public :: moddiel           ! Precondition the residual, using a model dielectric function.
!!***

contains
!!***

!!****f* ABINIT/prcref
!!
!! NAME
!! prcref
!!
!! FUNCTION
!! Compute preconditioned residual potential (or density) and forces.
!! iprcel, densfor_pred and iprcfc govern the choice of the preconditioner.
!! Three tasks are done:
!! 1) Preconditioning of the forces (residual has already been included)
!!     using the approximate force constant matrix. Get proposed
!!     change of atomic positions.
!! 2) Precondition the residual, get first part of proposed trial
!!     potential change.
!! 3) PAW only: precondition the rhoij residuals (simple preconditionning)
!! 4) Take into account the proposed change of atomic positions to
!!     modify the proposed trial potential change.
!!
!! NOTE
!! This routine is almost similar to prcref_PMA.F90 which is employed in
!! case of potential mixing. Yet it has undergone strong changes simultaneously
!! from two different sources at the same time which resulted in a splitting.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  dielstrt=number of the step at which the dielectric preconditioning begins.
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | intxc=control xc quadrature
!!   | densfor_pred= not yet used here
!!   | iprcel= governs the preconditioning of the potential residual
!!   |    0 => simple model dielectric matrix, described by the
!!   |              parameters dielng, diemac, diemix and diemixmag contained in dielar.
!!   |    between 21 and 39 => until istep=dielstart, same as iprcel=0, then uses
!!   |              the RPA dielectric matrix (routine dielmt)
!!   |    between 41 and 49 => uses the RPA dielectric matrix (routine dielmt).
!!   |    between 51 and 59 => uses the RPA dielectric matrix (routine dieltcel).
!!   |    between 61 and 69 => uses the electronic dielectric matr (routine dieltcel).
!!   |    between 71 and 79 => uses the real-space preconditioner based on Kerker prc (prcrskerkerN)
!!   |    between 81 and 99 => reserved for futur version of the real-space preconditioner
!!   |    between 141 and 169 -> same as between 41 and 69 but with a different periodicity: modulo(iprcel modulo (10))
!!   | iprcfc= governs the preconditioning of the forces
!!   |         0 => hessian is the identity matrix
!!   |         1 => hessian is 0.5 times the identity matrix
!!   |         2 => hessian is 0.25 times the identity matrix
!!   | ixc=exchange-correlation choice parameter.
!!   | natom=number of atoms
!!   | nspden=number of spin-density components
!!   | occopt=option for occupancies
!!   | prtvol=control print volume and debugging
!!   | typat(natom)=integer type for each atom in cell
!!  etotal=total ennergy
!!  fcart(3,natom)=cartesian forces (hartree/bohr)
!!  ffttomix(nfft*(1-nfftprc/nfft))=Index of the points of the FFT (fine) grid on the grid used for mixing (coarse)
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  istep= number of the step in the SCF cycle
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  mgfft=maximum size of 1D FFTs
!!  moved_atm_inside= if 1, then the preconditioned forces
!!    as well as the preconditioned potential residual must be computed;
!!    otherwise, compute only the preconditioned potential residual.
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=number of fft grid points
!!  nfftprc=size of FFT grid on which the potential residual will be preconditionned
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ngfftprc(18)=contain all needed information about 3D FFT for the grid corresponding to nfftprc
!!  nkxc=second dimension of the array kxc, see rhotoxc.F90 for a description
!!  npawmix=-PAW only- number of spherical part elements to be mixed
!!  npwdiel=number of planewaves for dielectric matrix
!!  ntypat=number of types of atoms in cell.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  optreal=1 if residual potential is is REAL space, 2 if it is in RECIPROCAL SPACE
!!  optres=0: the array vresid contains a potential residual
!!         1: the array vresid contains a density residual
!!  pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!                                    Use here rhoij residuals (and gradients)
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=array for electron density in reciprocal space
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!  vresid(optreal*nfftprc,nspden)=residual potential
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree)
!!  vhartr(nfft)=array for holding Hartree potential
!!  vlspl(mqgrid,2,ntypat)=q^2 v(q) spline for each type of atom.
!!  vpsp(nfft)=array for holding local psp
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  dtn_pc(3,natom)=preconditioned change of atomic position,
!!                                          in reduced coordinates
!!  vrespc(optreal*nfftprc,nspden)=preconditioned residual of the potential
!!  ==== if psps%usepaw==1
!!    rhoijrespc(npawmix)= preconditionned rhoij residuals at output
!!
!! SIDE EFFECT
!!  dielinv(2,npwdiel,nspden,npwdiel,nspden)=
!!                              inverse of the dielectric matrix in rec. space
!!  kxc(nfft,nkxc)=exchange-correlation kernel,
!!       needed if the electronic dielectric matrix is computed
!!  ===== if densfor_pred==3 .and. moved_atm_inside==1 =====
!!    ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phases
!!
!! PARENTS
!!      newrho
!!
!! CHILDREN
!!      atm2fft,dielmt,dieltcel,fourdp,fresid,getph,hartre
!!      indirect_parallel_fourier,kgindex,mean_fftr,metric,mkcore,mklocl
!!      moddiel,prcrskerker1,prcrskerker2,rhotoxc,testsusmat,xcart2xred
!!      xcdata_init,xmpi_sum,zerosym
!!
!! SOURCE

subroutine prcref(atindx,dielar,dielinv,&
&  dielstrt,dtn_pc,dtset,etotal,fcart,ffttomix,gmet,gsqcut,&
&  istep,kg_diel,kxc,&
&  mgfft,moved_atm_inside,mpi_enreg,my_natom,&
&  nattyp,nfft,nfftprc,ngfft,ngfftprc,nkxc,npawmix,npwdiel,ntypat,n1xccc,&
&  optreal,optres,pawrhoij,pawtab,ph1d,psps,rhog,rhoijrespc,rhor,rprimd,&
&  susmat,vhartr,vpsp,vresid,vrespc,vxc,wvl,wvl_den,xred)

!Arguments-------------------------------
!scalars
 integer,intent(in) :: dielstrt,istep,my_natom,mgfft,moved_atm_inside,n1xccc
 integer,intent(in) :: nfft,nfftprc,nkxc,npawmix,npwdiel,ntypat,optreal,optres
 real(dp),intent(in) :: etotal,gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type), intent(in) :: wvl
 type(wvl_denspot_type), intent(inout) :: wvl_den
!arrays
 integer,intent(in) :: atindx(dtset%natom),ffttomix(nfft*(1-nfftprc/nfft))
 integer,intent(in) :: kg_diel(3,npwdiel),nattyp(ntypat),ngfft(18),ngfftprc(18)
 real(dp),intent(in) :: dielar(7),fcart(3,dtset%natom),rhog(2,nfft)
 real(dp),intent(in) :: rhor(nfft,dtset%nspden),rprimd(3,3)
 real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(in) :: vhartr(nfft),vresid(nfftprc*optreal,dtset%nspden)
 real(dp),intent(in) :: vxc(nfft,dtset%nspden)
 real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(inout) :: gmet(3,3),kxc(nfft,nkxc)
 real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom),vpsp(nfft)
 real(dp),intent(inout) :: xred(3,dtset%natom)
 real(dp),intent(out) :: dtn_pc(3,dtset%natom),rhoijrespc(npawmix)
 real(dp),intent(out) :: vrespc(nfftprc*optreal,dtset%nspden)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: coredens_method,cplex,dielop,iatom,ier,ifft,ii,index,ipw1
 integer :: ipw2,iq,iq0,ispden,klmn,kmix,n1,n2,n3,n3xccc,nfftot,nk3xc,optatm
 integer :: optdyfr,opteltfr,optgr,option,optn,optn2,optstr,optv,vloc_method
 real(dp) :: ai,ar,diemix,diemixmag,eei,enxc
 real(dp) :: mixfac
 real(dp) :: mixfac_eff,mixfacmag,ucvol,vxcavg
 logical :: computediel,non_magnetic_xc
 character(len=500) :: message
 type(xcdata_type) :: xcdata
!arrays
 integer :: qprtrb(3)
 integer,allocatable :: indpw_prc(:)
 real(dp) :: dummy6(6),dummy7(6),gprimd(3,3),qphon(3),rmet(3,3),strsxc(6)
 real(dp) :: vmean(dtset%nspden),vprtrb(2)
 real(dp),allocatable :: dummy(:),dummy1(:),dummy2(:),dummy3(:),dummy4(:),dummy5(:),dummy8(:),dummy9(:)
 real(dp),allocatable :: dyfrlo_indx(:,:,:),dyfrx2(:,:,:)
 real(dp),allocatable :: fcart_pc(:,:),gresid(:,:),grtn_indx(:,:)
 real(dp),allocatable :: grxc(:,:),grxc_indx(:,:),rhog_wk(:,:),rhor_new(:,:)
 real(dp),allocatable :: rhor_wk(:,:),rhor_wk0(:,:),vhartr_wk(:),vpsp_wk(:)
 real(dp),allocatable :: vres_diel(:,:),vxc_wk(:,:),work(:),work1(:,:),work2(:)
 real(dp),allocatable :: work3(:,:),xccc3d(:),xred_wk(:,:)
 logical,allocatable :: mask(:)

! *************************************************************************

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!1) Eventually take care of the forces

 if(moved_atm_inside==1)then
   ABI_ALLOCATE(fcart_pc,(3,dtset%natom))

   if(dtset%iprcfc==0)then
     fcart_pc(:,:)=fcart(:,:)
   else
     fcart_pc(:,:)= (two**dtset%iprcfc) * fcart(:,:)
   end if

!  Compute preconditioned delta xred from preconditioned fcart and rprimd
   call xcart2xred(dtset%natom,rprimd,fcart_pc,dtn_pc)

   ABI_DEALLOCATE(fcart_pc)
 end if

!#######################################################################

!2) Take care of the potential residual

!Compute the residuals corresponding to the solution
!of an approximate realspace dielectric function according
!to X. Gonze PRB vol54 nb7 p4383 (1996) [[cite:Gonze1996]]
 if(dtset%iprcel>=71.and.dtset%iprcel<=79) then
   if (nfft==nfftprc) then
     if (dtset%iprcel<=78) then
       call prcrskerker1(dtset,mpi_enreg,nfft,dtset%nspden,ngfft,dielar,etotal, &
&       gprimd,vresid,vrespc,rhor(:,1))
     else
       call prcrskerker2(dtset,nfft,dtset%nspden,ngfft,dielar,gprimd,rprimd, &
&       vresid,vrespc,dtset%natom,xred,mpi_enreg,ucvol)
     end if
   else
!    If preconditionning has to be done on a coarse grid,
!    has to transfer several arrays
     ABI_ALLOCATE(work1,(nfftprc,dtset%nspden))
     ABI_ALLOCATE(work3,(nfftprc,dtset%nspden))
     ABI_ALLOCATE(work,(2*nfftprc))
     do ispden=1,dtset%nspden
       work(:)=vresid(:,ispden)
       call fourdp(1,work,work1(:,ispden),+1,mpi_enreg,nfftprc,1,ngfftprc,0)
     end do
     ABI_DEALLOCATE(work)
     if (dtset%iprcel<=78) then
       ABI_ALLOCATE(rhog_wk,(2,nfftprc))
       rhog_wk(:,:)=zero
       if (mpi_enreg%nproc_fft>1.and. mpi_enreg%paral_kgb==1) then
         nfftot=PRODUCT(ngfft(1:3))
         call indirect_parallel_Fourier(ffttomix,rhog_wk,mpi_enreg,ngfftprc,&
&         ngfft,nfftprc,nfft,dtset%paral_kgb,rhog,nfftot)
       else
         do ii=1,nfft
           if (ffttomix(ii)>0) rhog_wk(:,ffttomix(ii))=rhog(:,ii)
         end do
       end if
       call zerosym(rhog_wk,2,ngfftprc(1),ngfftprc(2),ngfftprc(3),&
&       comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
       ABI_ALLOCATE(work,(nfftprc))
       call fourdp(1,rhog_wk,work,+1,mpi_enreg,nfftprc,1,ngfftprc,0)
       call prcrskerker1(dtset,mpi_enreg,nfftprc,dtset%nspden,ngfftprc,dielar,etotal, &
&       gprimd,work1,work3,work)
       ABI_DEALLOCATE(work)
     else
       call prcrskerker2(dtset,nfftprc,dtset%nspden,ngfftprc,dielar,gprimd,rprimd, &
&       work1,work3,dtset%natom,xred,mpi_enreg,ucvol)
     end if
     do ispden=1,dtset%nspden
       call fourdp(1,vrespc(:,ispden),work3(:,ispden),-1,mpi_enreg,nfftprc,1,ngfftprc,0)
     end do
     ABI_DEALLOCATE(work1)
     ABI_DEALLOCATE(work3)
   end if

 else

   if(dtset%iprcel==0 .or. (dtset%iprcel<40.and.istep<dielstrt) )then
     cplex=optreal
     qphon(:)=zero
!    Simple scalar multiplication, or model dielectric function
     call moddiel(cplex,dielar,mpi_enreg,nfftprc,ngfftprc,dtset%nspden,optreal,optres,qphon,rprimd,vresid,vrespc)

!    Use the inverse dielectric matrix in a small G sphere
   else if( (istep>=dielstrt .and. dtset%iprcel>=21) .or. modulo(dtset%iprcel,100)>=41 )then

!    With dielop=1, the matrices will be computed when istep=dielstrt
!    With dielop=2, the matrices will be computed when istep=dielstrt and 1
     dielop=1
     if(modulo(dtset%iprcel,100)>=41)dielop=2
     call testsusmat(computediel,dielop,dielstrt,dtset,istep) !test if the matrix is to be computed
     if(computediel) then
!      Compute the inverse dielectric matrix from the susceptibility matrix
!      There are two routines for the RPA matrix, while for the electronic
!      dielectric matrix, only dieltcel will do the work
       if(modulo(dtset%iprcel,100)<=49)then
         call dielmt(dielinv,gmet,kg_diel,&
&         npwdiel,dtset%nspden,dtset%occopt,dtset%prtvol,susmat)
       else
         option=1
         if(modulo(dtset%iprcel,100)>=61)option=2
         call dieltcel(dielinv,gmet,kg_diel,kxc,&
&         nfft,ngfft,nkxc,npwdiel,dtset%nspden,dtset%occopt,option,dtset%prtvol,susmat)
       end if
     end if

     ABI_ALLOCATE(work1,(2,nfftprc))
     ABI_ALLOCATE(work2,(optreal*nfftprc))

!    Presently, one uses the inverse of the RPA dielectric matrix,
!    for which spin must be averaged.

!    Do fft from real space (work2) to G space (work1)
     if (optreal==1) then
       work2(:)=vresid(:,1)
!      Must average over spins in the case of a potential residual
       if(dtset%nspden/=1.and.optres==0)work2(:)=(work2(:)+vresid(:,2))*half
       call fourdp(1,work1,work2,-1,mpi_enreg,nfftprc,1,ngfftprc,0)
     else
       work1(:,:)=reshape(vresid(:,1),(/2,nfftprc/))
       if(dtset%nspden/=1.and.optres==0)work1(:,:)=(work1(:,:)+reshape(vresid(:,2),(/2,nfftprc/)))*half
     end if

!    Multiply by restricted inverse of dielectric matrix.
!    Must first copy relevant elements of work1 to a npwdiel-dimensioned array,
!    then zero work1, operate with the dielinv matrix, and store in work1.

     ABI_ALLOCATE(vres_diel,(2,npwdiel))
     ABI_ALLOCATE(indpw_prc,(npwdiel))
     ABI_ALLOCATE(mask,(npwdiel))
     mask(:)=.true.
     call kgindex(indpw_prc,kg_diel,mask,mpi_enreg,ngfftprc,npwdiel)

     do ipw1=1,npwdiel
       if(mask(ipw1)) then
         vres_diel(1,ipw1)=work1(1,indpw_prc(ipw1))
         vres_diel(2,ipw1)=work1(2,indpw_prc(ipw1))
       end if
     end do
     work1(:,:)=zero
     do ipw1=1,npwdiel
       ar=zero ; ai=zero

!      Use inverse of dielectric matrix (potential mixing)
       if (optres==0) then
         do ipw2=1,npwdiel
           if(mask(ipw2))then
             ar=ar+dielinv(1,ipw1,1,ipw2,1)*vres_diel(1,ipw2) &
&             -dielinv(2,ipw1,1,ipw2,1)*vres_diel(2,ipw2)
             ai=ai+dielinv(2,ipw1,1,ipw2,1)*vres_diel(1,ipw2) &
&             +dielinv(1,ipw1,1,ipw2,1)*vres_diel(2,ipw2)
           end if
         end do
       else
!        Use symetric of inverse of dielectric matrix (density mixing)
         do ipw2=1,npwdiel
           if(mask(ipw2))then
             ar=ar+dielinv(1,ipw2,1,ipw1,1)*vres_diel(1,ipw2) &
&             +dielinv(2,ipw2,1,ipw1,1)*vres_diel(2,ipw2)
             ai=ai-dielinv(2,ipw2,1,ipw1,1)*vres_diel(1,ipw2) &
&             +dielinv(1,ipw2,1,ipw1,1)*vres_diel(2,ipw2)
           end if
         end do
       end if
!      Must be careful not to count the diagonal 1 twice : it is added later,
!      so must be subtracted now.
       call xmpi_sum(ar,mpi_enreg%comm_fft,ier)
       call xmpi_sum(ai,mpi_enreg%comm_fft,ier)
       if(mask(ipw1)) then
         work1(1,indpw_prc(ipw1))=ar-vres_diel(1,ipw1)
         work1(2,indpw_prc(ipw1))=ai-vres_diel(2,ipw1)
       end if !mask(ipw1)
     end do ! ipw1
     ABI_DEALLOCATE(vres_diel)
     ABI_DEALLOCATE(indpw_prc)
     ABI_DEALLOCATE(mask)
!    Fourier transform
     if (optreal==1) then
       call fourdp(1,work1,work2,1,mpi_enreg,nfftprc,1,ngfftprc,0)
     else
       work2(:)=reshape(work1(:,:),(/nfftprc*2/))
     end if

!    Add to get the preconditioned vresid, must be careful about spins.
     if(dtset%iprcel>=30)then
       diemix=dielar(4);diemixmag=abs(dielar(7))
       vrespc(:,1)=diemix*(vresid(:,1)+work2(:))
       if(dtset%nspden/=1.and.optres==0)vrespc(:,2)=diemixmag*(vresid(:,2)+work2(:))
       if(dtset%nspden==4.and.optres==0)vrespc(:,3:4)=diemixmag*vresid(:,3:4)
       if(dtset%nspden/=1.and.optres==1)vrespc(:,2:dtset%nspden)=diemixmag*vresid(:,2:dtset%nspden)
     else
       vrespc(:,1)=vresid(:,1)+work2(:)
       if(dtset%nspden/=1.and.optres==0)vrespc(:,2)=vresid(:,2)+work2(:)
       if(dtset%nspden==4.and.optres==0)vrespc(:,3:4)=vresid(:,3:4)
       if(dtset%nspden/=1.and.optres==1)vrespc(:,2:dtset%nspden)=vresid(:,2:dtset%nspden)
     end if

     ABI_DEALLOCATE(work1)
     ABI_DEALLOCATE(work2)

!    Other choice ?
   else
     write(message, '(a,i3,a,a,a,a)' )&
&     'From the calling routine, iprcel=',dtset%iprcel,ch10,&
&     'The only allowed values are 0 or larger than 20.',ch10,&
&     'Action: correct your input file.'
     MSG_ERROR(message)
   end if
 end if
!#######################################################################

!3) PAW only : precondition the rhoij quantities (augmentation
!occupancies) residuals. Use a simple preconditionning
!with the same mixing factor as the model dielectric function.

 if (psps%usepaw==1.and.my_natom>0) then
   if (istep>=dielstrt.and.dtset%iprcel>=21.and.dtset%iprcel<30) then
     mixfac=one;mixfacmag=one
   else
     mixfac=dielar(4);mixfacmag=abs(dielar(7))
   end if
   if (pawrhoij(1)%cplex_rhoij==1) then
     index=0
     do iatom=1,my_natom
       do iq=1,pawrhoij(iatom)%qphase
         iq0=merge(0,pawrhoij(iatom)%lmn2_size,iq==1)
         do ispden=1,pawrhoij(iatom)%nspden
           mixfac_eff=mixfac;if (ispden>1) mixfac_eff=mixfacmag
           do kmix=1,pawrhoij(iatom)%lmnmix_sz
             index=index+1;klmn=iq0+pawrhoij(iatom)%kpawmix(kmix)
             rhoijrespc(index)=mixfac_eff*pawrhoij(iatom)%rhoijres(klmn,ispden)
           end do
         end do
       end do
     end do
   else
     index=-1
     do iatom=1,my_natom
       do iq=1,pawrhoij(iatom)%qphase
         iq0=merge(0,2*pawrhoij(iatom)%lmn2_size,iq==1)
         do ispden=1,pawrhoij(iatom)%nspden
           mixfac_eff=mixfac;if (ispden>1) mixfac_eff=mixfacmag
           do kmix=1,pawrhoij(iatom)%lmnmix_sz
             index=index+2;klmn=iq0+2*pawrhoij(iatom)%kpawmix(kmix)-1
             rhoijrespc(index:index+1)=mixfac_eff*pawrhoij(iatom)%rhoijres(klmn:klmn+1,ispden)
           end do
         end do
       end do
     end do
   end if
 end if

!#######################################################################

!4) Take care of the change of atomic positions
!Note : this part is very demanding on memory...
!however, since this algorithm is still in development,
!it was NOT included in the estimation provided by memory.f
 if(abs(dtset%densfor_pred)==3 .and. moved_atm_inside==1)then

!  Not yet compatible with resid given in reciprocal space
   if (optreal/=1) then
     write(message, '(5a)' )&
&     'From the calling routine, densfor_pred=3',ch10,&
&     'You cannot use residuals in reciprocal space.',ch10,&
&     'Action: correct your input file.'
     MSG_ERROR(message)
   end if
!  Not compatible with non-collinear magnetism
   if(dtset%nspden==4)then
     MSG_ERROR('densfor_pred=3 does not work for nspden=4 !')
   end if

   n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
   nfftot=PRODUCT(ngfft(1:3))

   if (optres==0) then  ! Array vresid contains a potential residual
!    -----------------------------------------------------------------

!    First subtract the current local, hartree and exchange correlation potentials
     do ispden=1,min(dtset%nspden,2)
       vrespc(:,ispden)=vrespc(:,ispden)-vpsp(:)-vhartr(:)-vxc(:,ispden)
     end do
     if (dtset%nspden==4) then
       do ispden=3,4
         vrespc(:,ispden)=vrespc(:,ispden)-vxc(:,ispden)
       end do
     end if

!    Compute the modified density, in rhor_wk
     option=2
     ABI_ALLOCATE(gresid,(3,dtset%natom))
     ABI_ALLOCATE(grxc,(3,dtset%natom))
     ABI_ALLOCATE(rhor_wk,(nfft,dtset%nspden))
     ABI_ALLOCATE(rhor_wk0,(nfft,dtset%nspden))
     ABI_ALLOCATE(xred_wk,(3,dtset%natom))
     xred_wk(:,:)=xred(:,:)+dtn_pc(:,:)
     call fresid(dtset,gresid,mpi_enreg,nfft,ngfft,&
&     ntypat,option,pawtab,rhor,rprimd,&
&     ucvol,rhor_wk,xred_wk,xred,psps%znuclpsp)

!    Compute up+down rhog_wk(G) by fft
     ABI_ALLOCATE(work,(nfft))
     ABI_ALLOCATE(rhog_wk,(2,nfft))
     work(:)=rhor_wk(:,1)
     call fourdp(1,rhog_wk,work,-1,mpi_enreg,nfft,1,ngfft,0)
     ABI_DEALLOCATE(work)

!    Compute structure factor phases for new atomic pos:
     call getph(atindx,dtset%natom,n1,n2,n3,ph1d,xred_wk)

!    Compute local ionic pseudopotential vpsp:
!    and core electron density xccc3d, if needed.
     n3xccc=0;if (n1xccc/=0) n3xccc=nfft
     ABI_ALLOCATE(xccc3d,(n3xccc))
     ABI_ALLOCATE(vpsp_wk,(nfft))
     vprtrb(1:2)=zero

!    Determine by which method the local ionic potential and/or
!    the pseudo core charge density contributions have to be computed
!    Local ionic potential:
!     Method 1: PAW
!     Method 2: Norm-conserving PP, icoulomb>0, wavelets
     vloc_method=1;if (psps%usepaw==0) vloc_method=2
     if (dtset%icoulomb>0) vloc_method=2
     if (psps%usewvl==1) vloc_method=2
!    Pseudo core charge density:
!     Method 1: PAW, nc_xccc_gspace
!     Method 2: Norm-conserving PP, wavelets
     coredens_method=1;if (psps%usepaw==0) coredens_method=2
     if (psps%nc_xccc_gspace==1) coredens_method=1
     if (psps%nc_xccc_gspace==0) coredens_method=2
     if (psps%usewvl==1) coredens_method=2

!    Local ionic potential and/or pseudo core charge by method 1
     if (vloc_method==1.or.coredens_method==1) then
       optv=0;if (vloc_method==1) optv=1
       optn=0;if (coredens_method==1) optn=n3xccc/nfft
       optatm=1;optdyfr=0;optgr=0;optstr=0;optn2=1;opteltfr=0
!      Note: atindx1 should be passed to atm2fft (instead of atindx) but it is unused...
       call atm2fft(atindx,xccc3d,vpsp,dummy,dummy2,dummy9,dummy1,gmet,gprimd,dummy3,dummy4,gsqcut,&
&       mgfft,psps%mqgrid_vl,dtset%natom,nattyp,nfft,ngfft,ntypat,optatm,optdyfr,opteltfr,optgr,optn,optn2,&
&       optstr,optv,psps,pawtab,ph1d,psps%qgrid_vl,qprtrb,dummy5,dummy6,dummy7,&
&       ucvol,psps%usepaw,dummy8,dummy8,dummy8,vprtrb,psps%vlspl,&
&       comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,&
&       paral_kgb=mpi_enreg%paral_kgb,distribfft=mpi_enreg%distribfft)
     end if

!    Local ionic potential by method 2
     if (vloc_method==2) then
       option=1
       ABI_ALLOCATE(dyfrlo_indx,(3,3,dtset%natom))
       ABI_ALLOCATE(grtn_indx,(3,dtset%natom))
       call mklocl(dtset,dyfrlo_indx,eei,gmet,gprimd,grtn_indx,gsqcut,dummy6,&
&       mgfft,mpi_enreg,dtset%natom,nattyp,nfft,ngfft,dtset%nspden,&
&       ntypat,option,pawtab,ph1d,psps,qprtrb,rhog_wk,rhor_wk,rprimd,&
&       ucvol,vprtrb,vpsp_wk,wvl,wvl_den,xred)
       ABI_DEALLOCATE(dyfrlo_indx)
       ABI_DEALLOCATE(grtn_indx)
     end if

!    Pseudo core electron density by method 2
     if (coredens_method==2.and.n1xccc/=0) then
       option=1
       ABI_ALLOCATE(dyfrx2,(3,3,dtset%natom))
       ABI_ALLOCATE(grxc_indx,(3,dtset%natom))
       call mkcore(dummy6,dyfrx2,grxc_indx,mpi_enreg,dtset%natom,nfft,dtset%nspden,ntypat,&
&       n1,n1xccc,n2,n3,option,rprimd,dtset%typat,ucvol,vxc,psps%xcccrc,&
&       psps%xccc1d,xccc3d,xred_wk)
       ABI_DEALLOCATE(dyfrx2)
       ABI_DEALLOCATE(grxc_indx)
     end if

!    Compute Hartree+xc potentials
     ABI_ALLOCATE(vxc_wk,(nfft,dtset%nspden))
     ABI_ALLOCATE(vhartr_wk,(nfft))
     option=1

     call hartre(1,gsqcut,psps%usepaw,mpi_enreg,nfft,ngfft,rhog_wk,rprimd,vhartr_wk)

!    Prepare the call to rhotoxc
     call xcdata_init(xcdata,dtset=dtset)
     nk3xc=1 ; non_magnetic_xc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)
     call rhotoxc(enxc,kxc,mpi_enreg,nfft,ngfft,&
&     work,0,work,0,nkxc,nk3xc,non_magnetic_xc,n3xccc,option,rhor_wk,rprimd,strsxc,1,&
&     vxc_wk,vxcavg,xccc3d,xcdata,vhartr=vhartr_wk)
     ABI_DEALLOCATE(xccc3d)

!    Sum all contributions
     do ispden=1,min(dtset%nspden,2)
       do ifft=1,nfft
         vrespc(ifft,ispden)=vrespc(ifft,ispden)+vpsp_wk(ifft)+vhartr_wk(ifft)+vxc_wk(ifft,ispden)
       end do
     end do
     if (dtset%nspden==4) then
       do ispden=3,4
         do ifft=1,nfft
           vrespc(ifft,ispden)=vrespc(ifft,ispden)+vxc_wk(ifft,ispden)
         end do
       end do
     end if
     call mean_fftr(vrespc,vmean,nfft,nfftot,dtset%nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
     if(dtset%nspden==2) then
       vmean(1)=half*(vmean(1)+vmean(2))
       vmean(2)=vmean(1)
     end if
     do ispden=1,dtset%nspden
       vrespc(:,ispden)=vrespc(:,ispden)-vmean(ispden)
     end do
     ABI_DEALLOCATE(gresid)
     ABI_DEALLOCATE(grxc)
     ABI_DEALLOCATE(rhog_wk)
     ABI_DEALLOCATE(rhor_wk)
     ABI_DEALLOCATE(rhor_wk0)
     ABI_DEALLOCATE(xred_wk)
     ABI_DEALLOCATE(vhartr_wk)
     ABI_DEALLOCATE(vpsp_wk)
     ABI_DEALLOCATE(vxc_wk)

   else                 ! Array vresid contains a density residual
!    -----------------------------------------------------------------

!    Only have to compute the modified preconditionned density residual
     option=2
     ABI_ALLOCATE(gresid,(3,dtset%natom))
     ABI_ALLOCATE(grxc,(3,dtset%natom))
     ABI_ALLOCATE(rhor_new,(nfft,dtset%nspden))
     ABI_ALLOCATE(rhor_wk,(nfft,dtset%nspden))
     ABI_ALLOCATE(rhor_wk0,(nfft,dtset%nspden))
     ABI_ALLOCATE(xred_wk,(3,dtset%natom))
     xred_wk(:,:)=xred(:,:)+dtn_pc(:,:)
     rhor_new(:,1)=rhor(:,1)+vrespc(:,1)
     if (dtset%nspden==2) then
       rhor_new(:,1)=rhor_new(:,1)+vrespc(:,2)
       rhor_new(:,2)=rhor(:,2)+vrespc(:,1)
     end if
     call fresid(dtset,gresid,mpi_enreg,nfft,ngfft,&
&     ntypat,option,pawtab,rhor,rprimd,&
&     ucvol,rhor_wk0,xred_wk,xred,psps%znuclpsp)
     call fresid(dtset,gresid,mpi_enreg,nfft,ngfft,&
&     ntypat,option,pawtab,rhor_new,rprimd,&
&     ucvol,rhor_wk,xred_wk,xred,psps%znuclpsp)
     vrespc(:,1)=rhor_wk(:,dtset%nspden)-rhor_wk0(:,dtset%nspden)
     if (dtset%nspden==2) vrespc(:,2)=rhor_wk(:,1)-rhor_wk0(:,1)-vrespc(:,1)
     ABI_DEALLOCATE(gresid)
     ABI_DEALLOCATE(grxc)
     ABI_DEALLOCATE(rhor_new)
     ABI_DEALLOCATE(rhor_wk)
     ABI_DEALLOCATE(rhor_wk0)
     ABI_DEALLOCATE(xred_wk)
   end if

 end if

end subroutine prcref
!!***

!!****f* ABINIT/prcref_PMA
!!
!! NAME
!! prcref_PMA
!!
!! FUNCTION
!! Compute preconditioned residual potential (or density) and forces.
!! iprcel, densfor_pred and iprcfc govern the choice of the preconditioner.
!! Three tasks are done:
!! 1) Preconditioning of the forces (residual has already been included)
!!     using the approximate force constant matrix. Get proposed
!!     change of atomic positions.
!! 2) Precondition the residual, get first part of proposed trial
!!     potential change.
!! 3) PAW only: precondition the rhoij residuals (simple preconditionning)
!! 4) Take into account the proposed change of atomic positions to
!!     modify the proposed trial potential change.
!!
!! NOTE
!! This routine is almost similar to prcref.F90 which is employed in
!! case of density mixing. Yet it has undergone strong changes simultaneously
!! from two different sources at the same time which resulted in a splitting.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  dielstrt=number of the step at which the dielectric preconditioning begins.
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | intxc=control xc quadrature
!!   | densfor_pred= not yet used here
!!   | iprcel= governs the preconditioning of the potential residual
!!   |    0 => simple model dielectric matrix, described by the
!!   |              parameters dielng, diemac, diemix and diemixmag contained in dielar.
!!   |    between 21 and 39 => until istep=dielstart, same as iprcel=0, then uses
!!   |              the RPA dielectric matrix (routine dielmt)
!!   |    between 41 and 49 => uses the RPA dielectric matrix (routine dielmt).
!!   |    between 51 and 59 => uses the RPA dielectric matrix (routine dieltcel).
!!   |    between 61 and 69 => uses the electronic dielectric matr (routine dieltcel).
!!   |    between 71 and 79 => uses the real-space preconditioner based on Kerker prc (prcrskerkerN)
!!   |    between 81 and 99 => reserved for futur version of the real-space preconditioner
!!   |    between 141 and 169 -> same as between 41 and 69 but with a different periodicity: modulo(iprcel modulo (10))
!!   | iprcfc= governs the preconditioning of the forces
!!   |         0 => hessian is the identity matrix
!!   |         1 => hessian is 0.5 times the identity matrix
!!   |         2 => hessian is 0.25 times the identity matrix
!!   | ixc=exchange-correlation choice parameter.
!!   | natom=number of atoms
!!   | nspden=number of spin-density components
!!   | occopt=option for occupancies
!!   | prtvol=control print volume and debugging
!!   | typat(natom)=integer type for each atom in cell
!!  fcart(3,natom)=cartesian forces (hartree/bohr)
!!  ffttomix(nfft*(1-nfftprc/nfft))=Index of the points of the FFT (fine) grid on the grid used for mixing (coarse)
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  istep= number of the step in the SCF cycle
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  mgfft=maximum size of 1D FFTs
!!  moved_atm_inside= if 1, then the preconditioned forces
!!    as well as the preconditioned potential residual must be computed;
!!    otherwise, compute only the preconditioned potential residual.
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=number of fft grid points
!!  nfftprc=size of FFT grid on which the potential residual will be preconditionned
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ngfftprc(18)=contain all needed information about 3D FFT for the grid corresponding to nfftprc
!!  nkxc=second dimension of the array kxc, see rhotoxc.f for a description
!!  npawmix=-PAW only- number of spherical part elements to be mixed
!!  npwdiel=number of planewaves for dielectric matrix
!!  ntypat=number of types of atoms in cell.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  optreal=1 if residual potential is is REAL space, 2 if it is in RECIPROCAL SPACE
!!  optres=0: the array vresid contains a potential residual
!!         1: the array vresid contains a density residual
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!                                    Use here rhoij residuals (and gradients)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=array for electron density in reciprocal space
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!  vresid(optreal*nfftprc,nspden)=residual potential
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree)
!!  vhartr(nfft)=array for holding Hartree potential
!!  vlspl(mqgrid,2,ntypat)=q^2 v(q) spline for each type of atom.
!!  vpsp(nfft)=array for holding local psp
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!!  etotal
!!  pawtab
!!
!! OUTPUT
!!  dtn_pc(3,natom)=preconditioned change of atomic position,
!!                                          in reduced coordinates
!!  vrespc(optreal*nfftprc,nspden)=preconditioned residual of the potential
!!  ==== if psps%usepaw==1
!!    rhoijrespc(npawmix)= preconditionned rhoij residuals at output
!!
!! SIDE EFFECT
!!  dielinv(2,npwdiel,nspden,npwdiel,nspden)=
!!                              inverse of the dielectric matrix in rec. space
!!  kxc(nfft,nkxc)=exchange-correlation kernel,
!!       needed if the electronic dielectric matrix is computed
!!  ===== if densfor_pred==3 .and. moved_atm_inside==1 =====
!!    ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phases
!!
!! PARENTS
!!      newvtr
!!
!! CHILDREN
!!      atm2fft,dielmt,dieltcel,fourdp,fresid,getph,hartre
!!      indirect_parallel_fourier,kgindex,mean_fftr,metric,mkcore,mklocl
!!      moddiel,prcrskerker1,prcrskerker2,rhotoxc,testsusmat,xcart2xred
!!      xcdata_init,xmpi_sum,zerosym
!!
!! SOURCE

  subroutine prcref_PMA(atindx,dielar,dielinv,&
&  dielstrt,dtn_pc,dtset,fcart,ffttomix,gmet,gsqcut,&
&  istep,kg_diel,kxc,&
&  mgfft,moved_atm_inside,mpi_enreg,my_natom,&
&  nattyp,nfft,nfftprc,ngfft,ngfftprc,nkxc,npawmix,npwdiel,ntypat,n1xccc,&
&  optreal,optres,pawrhoij,ph1d,psps,rhog, rhoijrespc,rhor,rprimd,&
&  susmat,vhartr,vpsp,vresid,vrespc,vxc,xred,&
&  etotal,pawtab,wvl)

!Arguments-------------------------------
!variables used for tfvw
!scalars
 integer,intent(in) :: dielstrt,istep,mgfft,moved_atm_inside,my_natom,n1xccc
 integer,intent(in) :: nfft,nfftprc,nkxc,npawmix,npwdiel,ntypat
 integer,intent(in) :: optreal,optres
 real(dp),intent(in) :: etotal,gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_data), intent(inout) :: wvl
!arrays
 integer,intent(in) :: atindx(dtset%natom),ffttomix(nfft*(1-nfftprc/nfft))
 integer,intent(in) :: kg_diel(3,npwdiel),nattyp(ntypat),ngfft(18),ngfftprc(18)
 real(dp),intent(in) :: dielar(7),fcart(3,dtset%natom)
 real(dp),intent(in) :: rhog(2,nfft)
 real(dp),intent(in) :: rhor(nfft,dtset%nspden),rprimd(3,3)
 real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(in) :: vhartr(nfft),vresid(nfftprc*optreal,dtset%nspden)
 real(dp),intent(in) :: vxc(nfft,dtset%nspden)
 real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(inout) :: gmet(3,3),kxc(nfft,nkxc)
 real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom),vpsp(nfft)
 real(dp),intent(inout) :: xred(3,dtset%natom)
 real(dp),intent(out) :: dtn_pc(3,dtset%natom),rhoijrespc(npawmix)
 real(dp),intent(out) :: vrespc(nfftprc*optreal,dtset%nspden)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: coredens_method,cplex,dielop,iatom,ier,ifft,ii,index,ipw1
 integer :: ipw2,ispden,klmn,kmix,n1,n2,n3,n3xccc,nfftot,nk3xc,optatm
 integer :: optdyfr,opteltfr,optgr,option,optn,optn2,optstr,optv,vloc_method
 real(dp) :: ai,ar,diemix,diemixmag,eei,enxc
 real(dp) :: mixfac
 real(dp) :: mixfac_eff,mixfacmag,ucvol,vxcavg
 logical :: computediel,non_magnetic_xc
 character(len=500) :: message
 type(xcdata_type) :: xcdata
!arrays
 integer :: qprtrb(3)
 integer,allocatable :: indpw_prc(:)
 real(dp) :: dummy6(6),gprimd(3,3),qphon(3),rmet(3,3),strsxc(6)
 real(dp) :: vmean(dtset%nspden),vprtrb(2)
 real(dp),allocatable :: dummy_in(:)
 real(dp) :: dummy_out1(0),dummy_out2(0),dummy_out3(0),dummy_out4(0),dummy_out5(0),dummy_out6(0),dummy_out7(0)
 real(dp),allocatable :: dyfrlo_indx(:,:,:),dyfrx2(:,:,:)
 real(dp),allocatable :: fcart_pc(:,:),gresid(:,:),grtn_indx(:,:)
 real(dp),allocatable :: grxc(:,:),grxc_indx(:,:),rhog_wk(:,:)
 real(dp),allocatable :: rhor_wk(:,:),rhor_wk0(:,:),vhartr_wk(:),vpsp_wk(:)
 real(dp),allocatable :: vres_diel(:,:),vxc_wk(:,:),work(:),work1(:,:),work2(:)
 real(dp),allocatable :: work3(:,:),xccc3d(:),xred_wk(:,:)
 logical,allocatable :: mask(:)

! *************************************************************************

 if(optres==1)then
   MSG_ERROR('density mixing (optres=1) not admitted!')
 end if

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!1) Eventually take care of the forces

 if(moved_atm_inside==1)then
   ABI_ALLOCATE(fcart_pc,(3,dtset%natom))

   if(dtset%iprcfc==0)then
     fcart_pc(:,:)=fcart(:,:)
   else
     fcart_pc(:,:)= (two**dtset%iprcfc) * fcart(:,:)
   end if

!  Compute preconditioned delta xred from preconditioned fcart and rprimd
   call xcart2xred(dtset%natom,rprimd,fcart_pc,dtn_pc)

   ABI_DEALLOCATE(fcart_pc)
 end if

!#######################################################################

!2) Take care of the potential residual

!Compute the residuals corresponding to the solution
!of an approximate realspace dielectric function according
!to X. Gonze PRB vol54 nb7 p4383 (1996) [[cite:Gonze1996]]
 if(dtset%iprcel>=71.and.dtset%iprcel<=79) then
   if (nfft==nfftprc) then
     if (dtset%iprcel<=78) then
       call prcrskerker1(dtset,mpi_enreg,nfft,dtset%nspden,ngfft,dielar,etotal, &
&       gprimd,vresid,vrespc,rhor(:,1))
     else
       call prcrskerker2(dtset,nfft,dtset%nspden,ngfft,dielar,gprimd,rprimd, &
&       vresid,vrespc,dtset%natom,xred,mpi_enreg,ucvol)
     end if
   else
!    If preconditionning has to be done on a coarse grid,
!    has to transfer several arrays
     ABI_ALLOCATE(work1,(nfftprc,dtset%nspden))
     ABI_ALLOCATE(work3,(nfftprc,dtset%nspden))
     ABI_ALLOCATE(work,(2*nfftprc))
     do ispden=1,dtset%nspden
       work(:)=vresid(:,ispden)
       call fourdp(1,work,work1(:,ispden),+1,mpi_enreg,nfftprc,1,ngfftprc,0)
     end do
     ABI_DEALLOCATE(work)
     if (dtset%iprcel<=78) then
       ABI_ALLOCATE(rhog_wk,(2,nfftprc))
       rhog_wk(:,:)=zero
       if (mpi_enreg%nproc_fft>1.and. mpi_enreg%paral_kgb==1) then
         nfftot=PRODUCT(ngfft(1:3))
         call indirect_parallel_Fourier(ffttomix,rhog_wk,mpi_enreg,ngfftprc,&
&         ngfft,nfftprc,nfft,dtset%paral_kgb,rhog,nfftot)
       else
         do ii=1,nfft
           if (ffttomix(ii)>0) rhog_wk(:,ffttomix(ii))=rhog(:,ii)
         end do
       end if
       call zerosym(rhog_wk,2,ngfftprc(1),ngfftprc(2),ngfftprc(3),&
&       comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
       ABI_ALLOCATE(work,(nfftprc))
       call fourdp(1,rhog_wk,work,+1,mpi_enreg,nfftprc,1,ngfftprc,0)
       call prcrskerker1(dtset,mpi_enreg,nfftprc,dtset%nspden,ngfftprc,dielar,etotal, &
&       gprimd,work1,work3,work)
       ABI_DEALLOCATE(work)
     else
       call prcrskerker2(dtset,nfftprc,dtset%nspden,ngfftprc,dielar,gprimd,rprimd, &
&       work1,work3,dtset%natom,xred,mpi_enreg,ucvol)
     end if
     do ispden=1,dtset%nspden
       call fourdp(1,vrespc(:,ispden),work3(:,ispden),-1,mpi_enreg,nfftprc,1,ngfftprc,0)
     end do
     ABI_DEALLOCATE(work1)
     ABI_DEALLOCATE(work3)
   end if

 else

   if(dtset%iprcel==0 .or. (dtset%iprcel<40.and.istep<dielstrt) )then
     cplex=optreal
     qphon(:)=zero
!    Simple scalar multiplication, or model dielectric function
     call moddiel(cplex,dielar,mpi_enreg,nfftprc,ngfftprc,dtset%nspden,optreal,optres,qphon,rprimd,vresid,vrespc)

!    Use the inverse dielectric matrix in a small G sphere
   else if( (istep>=dielstrt .and. dtset%iprcel>=21) .or. modulo(dtset%iprcel,100)>=41 )then

!    Wnith dielop=1, the matrices will be computed when istep=dielstrt
!    With dielop=2, the matrices will be computed when istep=dielstrt and 1
     dielop=1
     if(modulo(dtset%iprcel,100)>=41)dielop=2
     call testsusmat(computediel,dielop,dielstrt,dtset,istep) !test if the matrix is to be computed
     if(computediel) then
!      Compute the inverse dielectric matrix from the susceptibility matrix
!      There are two routines for the RPA matrix, while for the electronic
!      dielectric matrix, only dieltcel will do the work
       if(modulo(dtset%iprcel,100)<=49)then
         call dielmt(dielinv,gmet,kg_diel,&
&         npwdiel,dtset%nspden,dtset%occopt,dtset%prtvol,susmat)
       else
         option=1
         if(modulo(dtset%iprcel,100)>=61)option=2
         call dieltcel(dielinv,gmet,kg_diel,kxc,&
&         nfft,ngfft,nkxc,npwdiel,dtset%nspden,dtset%occopt,option,dtset%prtvol,susmat)
       end if
     end if

     ABI_ALLOCATE(work1,(2,nfftprc))
     ABI_ALLOCATE(work2,(optreal*nfftprc))

!    Presently, one uses the inverse of the RPA dielectric matrix,
!    for which spin must be averaged.

!    Do fft from real space (work2) to G space (work1)
     if (optreal==1) then
       work2(:)=vresid(:,1)
!      Must average over spins if needed.
       if(dtset%nspden/=1)work2(:)=(work2(:)+vresid(:,2))*half
       call fourdp(1,work1,work2,-1,mpi_enreg,nfftprc,1,ngfftprc,0)
     else
       work1(:,:)=reshape(vresid(:,1),(/2,nfftprc/))
       if (dtset%nspden/=1) work1(:,:)=(work1(:,:)+reshape(vresid(:,2),(/2,nfftprc/)))*half
     end if

!    Multiply by restricted inverse of dielectric matrix.
!    Must first copy relevant elements of work1 to a npwdiel-dimensioned array,
!    then zero work1, operate with the dielinv matrix, and store in work1.

     ABI_ALLOCATE(vres_diel,(2,npwdiel))
     ABI_ALLOCATE(indpw_prc,(npwdiel))
     ABI_ALLOCATE(mask,(npwdiel))
     mask(:)=.true.
     call kgindex(indpw_prc,kg_diel,mask,mpi_enreg,ngfftprc,npwdiel)
     do ipw1=1,npwdiel
       if(mask(ipw1)) then
         vres_diel(1,ipw1)=work1(1,indpw_prc(ipw1))
         vres_diel(2,ipw1)=work1(2,indpw_prc(ipw1))
       end if
     end do

     work1(:,:)=zero
     do ipw1=1,npwdiel
       ar=zero ; ai=zero

!      Use inverse of dielectric matrix (potential mixing)
       do ipw2=1,npwdiel
         if(mask(ipw2))then
           ar=ar+dielinv(1,ipw1,1,ipw2,1)*vres_diel(1,ipw2) &
&           -dielinv(2,ipw1,1,ipw2,1)*vres_diel(2,ipw2)
           ai=ai+dielinv(2,ipw1,1,ipw2,1)*vres_diel(1,ipw2) &
&           +dielinv(1,ipw1,1,ipw2,1)*vres_diel(2,ipw2)
         end if
       end do
!      Must be careful not to count the diagonal 1 twice : it is added later,
!      so must be subtracted now.
       call xmpi_sum(ar,mpi_enreg%comm_fft,ier)
       call xmpi_sum(ai,mpi_enreg%comm_fft,ier)
       if(mask(ipw1)) then
         work1(1,indpw_prc(ipw1))=ar-vres_diel(1,ipw1)
         work1(2,indpw_prc(ipw1))=ai-vres_diel(2,ipw1)
       end if !mask(ipw1)
     end do ! ipw1
     ABI_DEALLOCATE(vres_diel)
     ABI_DEALLOCATE(indpw_prc)
     ABI_DEALLOCATE(mask)

!    Fourier transform
     if (optreal==1) then
       call fourdp(1,work1,work2,1,mpi_enreg,nfftprc,1,ngfftprc,0)
     else
       work2(:)=reshape(work1(:,:),(/nfftprc*2/))
     end if

!    Add to get the preconditioned vresid, must be careful about spins.
     if(dtset%iprcel>=30)then
       diemix=dielar(4);diemixmag=abs(dielar(7))
       vrespc(:,1)=diemix*(vresid(:,1)+work2(:))
       if(dtset%nspden/=1)vrespc(:,2)=diemixmag*(vresid(:,2)+work2(:))
       if(dtset%nspden==4)vrespc(:,3:4)=diemixmag*vresid(:,3:4)
     else
       vrespc(:,1)=vresid(:,1)+work2(:)
       if(dtset%nspden/=1)vrespc(:,2)=vresid(:,2)+work2(:)
       if(dtset%nspden==4)vrespc(:,3:4)=vresid(:,3:4)
     end if

     ABI_DEALLOCATE(work1)
     ABI_DEALLOCATE(work2)

!    Other choice ?
   else
     write(message, '(a,i0,a,a,a,a)' )&
&     'From the calling routine, iprcel= ',dtset%iprcel,ch10,&
&     'The only allowed values are 0 or larger than 20.',ch10,&
&     'Action: correct your input file.'
     MSG_ERROR(message)
   end if
 end if
!#######################################################################

!3) PAW only : precondition the rhoij quantities (augmentation
!occupancies) residuals. Use a simple preconditionning
!with the same mixing factor as the model dielectric function.

 if (psps%usepaw==1.and.my_natom>0) then
   if (istep>=dielstrt.and.dtset%iprcel>=21.and.dtset%iprcel<30) then
     mixfac=one;mixfacmag=one
   else
     mixfac=dielar(4);mixfacmag=abs(dielar(7))
   end if
   if (pawrhoij(1)%cplex_rhoij==1) then
     index=0
     do iatom=1,my_natom
       do ispden=1,pawrhoij(iatom)%nspden
         mixfac_eff=mixfac;if (ispden>1) mixfac_eff=mixfacmag
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           index=index+1;klmn=pawrhoij(iatom)%kpawmix(kmix)
           rhoijrespc(index)=mixfac_eff*pawrhoij(iatom)%rhoijres(klmn,ispden)
         end do
       end do
     end do
   else
     index=-1
     do iatom=1,my_natom
       do ispden=1,pawrhoij(iatom)%nspden
         mixfac_eff=mixfac;if (ispden>1) mixfac_eff=mixfacmag
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           index=index+2;klmn=2*pawrhoij(iatom)%kpawmix(kmix)-1
           rhoijrespc(index:index+1)=mixfac_eff*pawrhoij(iatom)%rhoijres(klmn:klmn+1,ispden)
         end do
       end do
     end do
   end if
 end if
!#######################################################################

!4) Take care of the change of atomic positions
!Note : this part is very demanding on memory...
!however, since this algorithm is still in development,
!it was NOT included in the estimation provided by memory.f
 if(abs(dtset%densfor_pred)==3 .and. moved_atm_inside==1)then

!  Not yet compatible with resid given in reciprocal space
   if (optreal/=1) then
     write(message, '(5a)' )&
&     'From the calling routine, densfor_pred=3',ch10,&
&     'You cannot use residuals in reciprocal space.',ch10,&
&     'Action: correct your input file.'
     MSG_ERROR(message)
   end if

!  Not compatible with non-collinear magnetism
   if(dtset%nspden==4)then
     MSG_ERROR('densfor_pred=3 does not work for nspden=4!')
   end if

   n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
   nfftot=PRODUCT(ngfft(1:3))

!  First subtract the current local, hartree and exchange correlation potentials
   do ispden=1,min(dtset%nspden,2)
     vrespc(:,ispden)=vrespc(:,ispden)-vpsp(:)-vhartr(:)-vxc(:,ispden)
   end do
   if (dtset%nspden==4) then
     do ispden=3,4
       vrespc(:,ispden)=vrespc(:,ispden)-vxc(:,ispden)
     end do
   end if

!  Compute the modified density, in rhor_wk
   option=2
   ABI_ALLOCATE(gresid,(3,dtset%natom))
   ABI_ALLOCATE(grxc,(3,dtset%natom))
   ABI_ALLOCATE(rhor_wk,(nfft,dtset%nspden))
   ABI_ALLOCATE(rhor_wk0,(nfft,dtset%nspden))
   ABI_ALLOCATE(xred_wk,(3,dtset%natom))
   xred_wk(:,:)=xred(:,:)+dtn_pc(:,:)

   call fresid(dtset,gresid,mpi_enreg,nfft,ngfft,&
&   ntypat,option,pawtab,rhor,rprimd,&
&   ucvol,rhor_wk,xred_wk,xred,psps%znuclpsp)

!  Compute up+down rhog_wk(G) by fft
   ABI_ALLOCATE(work,(nfft))
   ABI_ALLOCATE(rhog_wk,(2,nfft))
   work(:)=rhor_wk(:,1)
   call fourdp(1,rhog_wk,work,-1,mpi_enreg,nfft,1,ngfft,0)
   ABI_DEALLOCATE(work)

!  Compute structure factor phases for new atomic pos:
   call getph(atindx,dtset%natom,n1,n2,n3,ph1d,xred_wk)

!  Compute local ionic pseudopotential vpsp:
!  and core electron density xccc3d, if needed.
   n3xccc=0;if (n1xccc/=0) n3xccc=nfft
   ABI_ALLOCATE(xccc3d,(n3xccc))
   ABI_ALLOCATE(vpsp_wk,(nfft))
   vprtrb(1:2)=zero

!  Determine by which method the local ionic potential and/or
!  the pseudo core charge density have to be computed
!  Local ionic potential:
!   Method 1: PAW
!   Method 2: Norm-conserving PP, icoulomb>0, wavelets
   vloc_method=1;if (psps%usepaw==0) vloc_method=2
   if (dtset%icoulomb>0) vloc_method=2
   if (psps%usewvl==1) vloc_method=2
!  Pseudo core charge density:
!   Method 1: PAW, nc_xccc_gspace
!   Method 2: Norm-conserving PP, wavelets
   coredens_method=1;if (psps%usepaw==0) coredens_method=2
   if (psps%nc_xccc_gspace==1) coredens_method=1
   if (psps%nc_xccc_gspace==0) coredens_method=2
   if (psps%usewvl==1) coredens_method=2

!  Local ionic potential and/or pseudo core charge by method 1
   if (vloc_method==1.or.coredens_method==1) then
     optv=0;if (vloc_method==1) optv=1
     optn=0;if (coredens_method==1) optn=n3xccc/nfft
     optatm=1;optdyfr=0;opteltfr=0;optgr=0;optstr=0;optn2=1
!    Note: atindx1 should be passed to atm2fft (instead of atindx) but it is unused...
     call atm2fft(atindx,xccc3d,vpsp,dummy_out1,dummy_out2,dummy_out3,dummy_in,gmet,&
&     gprimd,dummy_out4,dummy_out5,gsqcut,mgfft,psps%mqgrid_vl,dtset%natom,nattyp,&
&     nfft,ngfft,ntypat,optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,&
&     psps,pawtab,ph1d,psps%qgrid_vl,qprtrb,dummy_in,dummy_out6,dummy_out7,ucvol,&
&     psps%usepaw,dummy_in,dummy_in,dummy_in,vprtrb,psps%vlspl,&
&     comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,&
&     paral_kgb=mpi_enreg%paral_kgb,distribfft=mpi_enreg%distribfft)
   end if

!  Local ionic potential by method 2
   if (vloc_method==2) then
     option=1
     ABI_ALLOCATE(dyfrlo_indx,(3,3,dtset%natom))
     ABI_ALLOCATE(grtn_indx,(3,dtset%natom))
     call mklocl(dtset,dyfrlo_indx,eei,gmet,gprimd,grtn_indx,gsqcut,dummy6,&
&     mgfft,mpi_enreg,dtset%natom,nattyp,nfft,ngfft,dtset%nspden,&
&     ntypat,option,pawtab,ph1d,psps,qprtrb,rhog_wk,rhor_wk,rprimd,&
&     ucvol,vprtrb,vpsp_wk,wvl%descr,wvl%den,xred)
     ABI_DEALLOCATE(dyfrlo_indx)
     ABI_DEALLOCATE(grtn_indx)
   end if

!  Pseudo core electron density by method 2
   if (coredens_method==2.and.n1xccc/=0) then
     option=1
     ABI_ALLOCATE(dyfrx2,(3,3,dtset%natom))
     ABI_ALLOCATE(grxc_indx,(3,dtset%natom))
     call mkcore(dummy6,dyfrx2,grxc_indx,mpi_enreg,dtset%natom,nfft,dtset%nspden,ntypat,&
&     n1,n1xccc,n2,n3,option,rprimd,dtset%typat,ucvol,vxc,psps%xcccrc,&
&     psps%xccc1d,xccc3d,xred_wk)
     ABI_DEALLOCATE(dyfrx2)
     ABI_DEALLOCATE(grxc_indx)
   end if

!  Compute Hartree+xc potentials
   ABI_ALLOCATE(vxc_wk,(nfft,dtset%nspden))
   ABI_ALLOCATE(vhartr_wk,(nfft))
   option=1

   call hartre(1,gsqcut,psps%usepaw,mpi_enreg,nfft,ngfft,rhog_wk,rprimd,vhartr_wk)

!  Prepare the call to rhotoxc
   call xcdata_init(xcdata,dtset=dtset)
   nk3xc=1 ; non_magnetic_xc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)
   ABI_ALLOCATE(work,(0))
   call rhotoxc(enxc,kxc,mpi_enreg,nfft,ngfft,&
&   work,0,work,0,nkxc,nk3xc,non_magnetic_xc,n3xccc,option,rhor_wk,rprimd,strsxc,1,&
&   vxc_wk,vxcavg,xccc3d,xcdata,vhartr=vhartr_wk)
   ABI_DEALLOCATE(work)
   ABI_DEALLOCATE(xccc3d)

!  Sum all contributions
   do ispden=1,min(dtset%nspden,2)
     do ifft=1,nfft
       vrespc(ifft,ispden)=vrespc(ifft,ispden)+vpsp_wk(ifft)+vhartr_wk(ifft)+vxc_wk(ifft,ispden)
     end do
   end do
   if (dtset%nspden==4) then
     do ispden=3,4
       do ifft=1,nfft
         vrespc(ifft,ispden)=vrespc(ifft,ispden)+vxc_wk(ifft,ispden)
       end do
     end do
   end if
   call mean_fftr(vrespc,vmean,nfft,nfftot,dtset%nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
   if(dtset%nspden==2) then
     vmean(1)=half*(vmean(1)+vmean(2))
     vmean(2)=vmean(1)
   end if
   do ispden=1,dtset%nspden
     vrespc(:,ispden)=vrespc(:,ispden)-vmean(ispden)
   end do
   ABI_DEALLOCATE(gresid)
   ABI_DEALLOCATE(grxc)
   ABI_DEALLOCATE(rhog_wk)
   ABI_DEALLOCATE(rhor_wk)
   ABI_DEALLOCATE(rhor_wk0)
   ABI_DEALLOCATE(xred_wk)
   ABI_DEALLOCATE(vhartr_wk)
   ABI_DEALLOCATE(vpsp_wk)
   ABI_DEALLOCATE(vxc_wk)

 end if

end subroutine prcref_PMA
!!***


!!****f* ABINIT/moddiel
!! NAME
!! moddiel
!!
!! FUNCTION
!! Precondition the residual, using a model dielectric function.
!! When cplex=1, assume q=(0 0 0), and vresid and vrespc will be REAL
!! When cplex=2, q must be taken into account, and vresid and vrespc will be COMPLEX
!!
!! INPUTS
!!  cplex= if 1, vhartr is REAL, if 2, vhartr is COMPLEX
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  optreal=1 if residual potential is in REAL space, 2 if it is in RECIPROCAL SPACE
!!  optres=0: the array vresid contains a potential residual
!!         1: the array vresid contains a density residual
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  vresid(cplex*nfft,nspden)=residual density/potential in REAL space       (if optreal==1)
!!                            residual density/potential in RECIPROCAL space (if optreal==2)
!!
!! OUTPUT
!!  vrespc(cplex*nfft,nspden)=preconditioned residual of the density/potential
!!                            in REAL space       if optreal==1
!!                            in RECIPROCAL space if optreal==2
!!
!! SIDE EFFECTS
!!
!! NOTES
!! optreal==2 is not compatible with cplex==1
!!
!! PARENTS
!!      dfpt_newvtr,prcref,prcref_PMA
!!
!! CHILDREN
!!      fourdp,metric,ptabs_fourdp
!!
!! SOURCE

subroutine moddiel(cplex,dielar,mpi_enreg,nfft,ngfft,nspden,optreal,optres,qphon,rprimd,vresid,vrespc)

!Arguments-------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,optreal,optres
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: dielar(7),qphon(3),rprimd(3,3)
 real(dp),intent(in) :: vresid(cplex*nfft,nspden)
 real(dp),intent(out) :: vrespc(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i23,i3,ifft,ig,ii,ii1,ing,ispden,me_fft,mg,n1,n2,n3,nproc_fft
 integer :: nspden_eff,qeq0
 logical :: magn_precon
 real(dp) :: dielng,diemac,diemac_inv,diemix,diemixmag,diemix_eff,factor,gqg2p3,gqgm12,gqgm13
 real(dp) :: gqgm23,gs,gs2,gs3,l2g2,length2,ucvol
 character(len=500) :: message
!arrays
 integer :: id(3)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: gmet(3,3),gprimd(3,3),potg0(4),rmet(3,3)
 real(dp),allocatable :: gq(:,:),work1(:,:),work2(:)

! *************************************************************************

!Check that cplex has an allowed value
 if(cplex/=1 .and. cplex/=2)then
   write(message,'(a,i0,a,a)')&
&   '  From the calling routine, cplex=',cplex,ch10,&
&   '  but the only value allowed are 1 and 2.'
   MSG_BUG(message)
 end if

 if(cplex==1.and.optreal==2)then
   MSG_BUG('When optreal=2, cplex must be 2.')
 end if

!This is to allow q=0
 qeq0=0
 if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15) qeq0=1

!If cplex=1 then qphon should be 0 0 0
 if (cplex==1.and. qeq0/=1) then
   write(message,'(a,3e12.4,a)' )' cplex=1 but qphon=',qphon,' qphon should be 0 0 0.'
   MSG_BUG(message)
 end if

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 me_fft=ngfft(11)
 nproc_fft=ngfft(10)

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 dielng=dielar(2) ; diemac=dielar(3) ; diemix=dielar(4) ; diemixmag=dielar(7)

 magn_precon=(diemixmag>=zero) ! Set to true if magnetization has to be preconditionned
 diemixmag=abs(diemixmag)

!write(std_out,*)' moddiel : diemac, diemix, diemixmag =',diemac,diemix,diemixmag

 if(abs(diemac-1.0_dp)<1.0d-6)then

!  Here, simple mixing is required, through macroscopic
!  dielectric constant set to 1.0_dp .
   vrespc(:,1)=diemix*vresid(:,1)
   if (nspden/=1) vrespc(:,2:nspden)=diemixmag*vresid(:,2:nspden)
 else

!  Magnetization is not preconditionned
   if (optres==1.and.nspden>1.and.(.not.magn_precon)) vrespc(:,2:nspden)=diemixmag*vresid(:,2:nspden)

!  Here, model dielectric function (G-diagonal operator)

   length2=(two_pi*dielng)**2
   diemac_inv=1.0_dp/diemac
   ABI_ALLOCATE(work1,(2,nfft))
   if (optreal==1) then
     ABI_ALLOCATE(work2,(cplex*nfft))
   end if

!  In order to speed the routine, precompute the components of g
   mg=maxval(ngfft)
   ABI_ALLOCATE(gq,(3,mg))
   do ii=1,3
     id(ii)=ngfft(ii)/2+2
     do ing=1,ngfft(ii)
       ig=ing-(ing/id(ii))*ngfft(ii)-1
       gq(ii,ing)=ig+qphon(ii)
     end do
   end do

!  Do-loop on spins
!  Note XG 010922 : I doubt the preconditioner is OK for the magnetization
   nspden_eff=nspden;if (optres==1.and.(.not.magn_precon)) nspden_eff=1
   do ispden=1,nspden_eff

     diemix_eff=diemix;if (ispden>1) diemix_eff=diemixmag

!    Do fft from real space (work2) to G space (work1)
     if (optreal==1) then
       work2(:)=vresid(:,ispden)
       call fourdp(cplex,work1,work2,-1,mpi_enreg,nfft,1,ngfft,0)
     else
!      work1(:,:)=reshape(vresid(:,ispden),(/2,nfft/))
!      Reshape function does not work with big arrays for some compilers
       do ifft=1,nfft
         work1(1,ifft)=vresid(2*ifft-1,ispden)
         work1(2,ifft)=vresid(2*ifft  ,ispden)
       end do
     end if

!    Store G=0 value
     potg0(ispden)=work1(re,1)

!    Triple loop, for the three dimensions
     do i3=1,n3
!      Precompute some products that do not depend on i2 and i1
       gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
       gqgm23=gq(3,i3)*gmet(2,3)*2
       gqgm13=gq(3,i3)*gmet(1,3)*2
       do i2=1,n2
         if (fftn2_distrib(i2)==me_fft) then
           gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
           gqgm12=gq(2,i2)*gmet(1,2)*2
           gqg2p3=gqgm13+gqgm12
           i23=n1*( ffti2_local(i2)-1+(n2/nproc_fft)*(i3-1))

!          Do the test that eliminates the Gamma point outside
!          of the inner loop
           ii1=1
           if(i2 == 1 .and. i3 == 1 .and. qeq0==1)then
!            if(i23==0 .and. qeq0==1)then: this changes with the number of fft procs...
!            and seems to be wrong.Pls check
             ii1=2
           end if

!          Here, unlike in hartre.f, the G=0 term is not eliminated, albeit
!          not changed.
           do i1=ii1,n1

!            One obtains the square of the norm of q+G (defined by i1,i2,i3)
             gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
             ifft=i1+i23

             l2g2=length2*gs
!            The model dielectric function is now computed
             factor = (l2g2+diemac_inv)/(l2g2+1.0_dp) * diemix_eff
             work1(re,ifft)=work1(re,ifft)*factor
             work1(im,ifft)=work1(im,ifft)*factor

           end do
         end if
       end do
     end do

!    Might get rid of the G=0 term
!    if(qeq0==1)then
!    work1(re,1)=0.0_dp
!    work1(im,1)=0.0_dp
!    end if

!    Fourier transform
     if (optreal==1) then
       call fourdp(cplex,work1,work2,1,mpi_enreg,nfft,1,ngfft,0)
       vrespc(:,ispden)=work2(:)
     else
!      vrespc(:,ispden)=reshape(work1(:,:),(/nfft*2/))
!      Reshape function does not work with big arrays for some compilers
       do ifft=1,nfft
         vrespc(2*ifft-1,ispden)=work1(1,ifft)
         vrespc(2*ifft  ,ispden)=work1(2,ifft)
       end do
     end if

!    End of loop on spin polarizations
   end do

   ABI_DEALLOCATE(gq)
   ABI_DEALLOCATE(work1)
   if (optreal==1) then
     ABI_DEALLOCATE(work2)
   end if

!  End condition diemac/=1.0
 end if

end subroutine moddiel
!!***

!!****f* ABINIT/dielmt
!! NAME
!! dielmt
!!
!! FUNCTION
!! Compute dielectric matrix from susceptibility matrix
!! Diagonalize it, then invert it.
!!
!! INPUTS
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  npwdiel=size of the dielinv and susmat arrays.
!!  nspden=number of spin-density components
!!  occopt=option for occupancies
!!  prtvol=control print volume and debugging output
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!
!! OUTPUT
!!  dielinv(2,npwdiel,(nspden+4)/3,npwdiel,(nspden+4)/3)=inverse of the (non-hermitian)
!!      TC dielectric matrix in reciprocal space.
!!
!! NOTES
!! Warning : will not work in the spin-polarized, metallic case.
!! Output (not cleaned)
!! !!! Spin behaviour is not obvious !!!
!!
!! TODO
!! Write equation below (hermitian matrix)
!!
!! PARENTS
!!      prcref,prcref_PMA
!!
!! CHILDREN
!!      timab,wrtout,zhpev
!!
!! SOURCE

subroutine dielmt(dielinv,gmet,kg_diel,npwdiel,nspden,occopt,prtvol,susmat)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwdiel,nspden,occopt,prtvol
!arrays
 integer,intent(in) :: kg_diel(3,npwdiel)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 real(dp),intent(out) :: dielinv(2,npwdiel,nspden,npwdiel,nspden)

!Local variables-------------------------------
!scalars
 integer :: ieig,ier,ii,index,ipw,ipw1,ipw2,isp,jj,npwsp
 real(dp) :: ai1,ai2,ar1,ar2,eiginv,gfact,gfactinv,gred1,gred2,gred3,gsquar
 real(dp) :: tpisq
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: dielh(:),dielmat(:,:,:,:,:),dielvec(:,:,:)
 real(dp),allocatable :: eig_diel(:),zhpev1(:,:),zhpev2(:)
!no_abirules
!integer :: ipw3
!real(dp) :: elementi,elementr

! *************************************************************************

!DEBUG
!write(std_out,*)' dielmt : enter '
!if(.true.)stop
!ENDDEBUG

!tpisq is (2 Pi) **2:
 tpisq=(two_pi)**2

 call timab(90,1,tsec)

!-Compute now the hermitian dielectric matrix------------------------------
!Following remarks are only valid within RPA approximation (Kxc=0):

!for the spin-unpolarized case, 1 - 4pi (1/G) chi0(G,Gp) (1/Gp)

!for the spin-polarized case,
!( 1  0 ) - 4pi ( 1/G  1/G )   ( chi0 upup  chi0 updn )   ( 1/Gp 1/Gp )
!( 0  1 )       ( 1/G  1/G )   ( chi0 dnup  chi0 dndn )   ( 1/Gp 1/Gp )
!which is equal to
!( 1  0 ) - 4pi (1/G  0 ) (chi0 upup+dndn+updn+dnup  chi0 upup+dndn+updn+dnup) (1/Gp 0  )
!( 0  1 )       ( 0  1/G) (chi0 upup+dndn+updn+dnup  chi0 upup+dndn+updn+dnup) ( 0  1/Gp)
!So, if spin-polarized, sum all spin contributions
!Note: chi0 updn = chi0 dnup = zero for non-metallic systems

!In the case of non-collinear magnetism, within RPA, this is the same because:
!chi0_(s1,s2),(s3,s4) = delta_s1,s2 * delta_s3,s4 * chi0_(s1,s1),(s3,s3)
!Only chi_upup,upup, chi_dndn,dndn, chi_upup,dndn and chi_dndn,upup
!have to be taken into account (stored, susmat(:,ipw1,1:2,ipw2,1:2)

 ABI_ALLOCATE(dielmat,(2,npwdiel,min(nspden,2),npwdiel,min(nspden,2)))

 if(nspden/=1)then
   if (occopt<3) then
     do ipw2=1,npwdiel
       do ipw1=1,npwdiel
         dielmat(1,ipw1,1,ipw2,1)=susmat(1,ipw1,1,ipw2,1)+susmat(1,ipw1,2,ipw2,2)
         dielmat(2,ipw1,1,ipw2,1)=susmat(2,ipw1,1,ipw2,1)+susmat(2,ipw1,2,ipw2,2)
       end do
     end do
   else
     do ipw2=1,npwdiel
       do ipw1=1,npwdiel
         dielmat(1,ipw1,1,ipw2,1)=susmat(1,ipw1,1,ipw2,1)+susmat(1,ipw1,2,ipw2,2)+susmat(1,ipw1,1,ipw2,2)+susmat(1,ipw1,2,ipw2,1)
         dielmat(2,ipw1,1,ipw2,1)=susmat(2,ipw1,1,ipw2,1)+susmat(2,ipw1,2,ipw2,2)+susmat(2,ipw1,1,ipw2,2)+susmat(2,ipw1,2,ipw2,1)
       end do
     end do
   end if
 else
   do ipw2=1,npwdiel
     do ipw1=1,npwdiel
       dielmat(1,ipw1,1,ipw2,1)=susmat(1,ipw1,1,ipw2,1)
       dielmat(2,ipw1,1,ipw2,1)=susmat(2,ipw1,1,ipw2,1)
     end do
   end do
 end if
!Compute 1/G factors and include them in the dielectric matrix
 do ipw1=1,npwdiel
   gred1=dble(kg_diel(1,ipw1))
   gred2=dble(kg_diel(2,ipw1))
   gred3=dble(kg_diel(3,ipw1))
   gsquar=tpisq*(gmet(1,1)*gred1**2+gmet(2,2)*gred2**2+gmet(3,3)*gred3**2 &
&   +two*( (gmet(1,2)*gred2+gmet(1,3)*gred3)* gred1 +      &
&   gmet(2,3)*gred2*gred3)                        )
!  Distinguish G=0 from other elements
   if(gsquar>tol12)then
!    !$ gfact=\sqrt (4.0_dp \pi/gsquar/dble(nspden))$
     gfact=sqrt(four_pi/gsquar)
     do ipw2=1,npwdiel
!      Must multiply both rows and columns, and also changes the sign
       dielmat(1,ipw2,1,ipw1,1)=-dielmat(1,ipw2,1,ipw1,1)*gfact
       dielmat(2,ipw2,1,ipw1,1)=-dielmat(2,ipw2,1,ipw1,1)*gfact
       dielmat(1,ipw1,1,ipw2,1)= dielmat(1,ipw1,1,ipw2,1)*gfact
       dielmat(2,ipw1,1,ipw2,1)= dielmat(2,ipw1,1,ipw2,1)*gfact
     end do
   else
!    Zero the G=0 elements, head and wings
     do ipw2=1,npwdiel
       dielmat(1,ipw2,1,ipw1,1)=zero
       dielmat(2,ipw2,1,ipw1,1)=zero
       dielmat(1,ipw1,1,ipw2,1)=zero
       dielmat(2,ipw1,1,ipw2,1)=zero
     end do
   end if
 end do

!Complete the matrix in the spin-polarized case
!should this be nspden==2??
 if(nspden/=1)then
   do ipw1=1,npwdiel
     do ipw2=1,npwdiel
       dielmat(1,ipw1,1,ipw2,2)=dielmat(1,ipw1,1,ipw2,1)
       dielmat(2,ipw1,1,ipw2,2)=dielmat(2,ipw1,1,ipw2,1)
       dielmat(1,ipw1,2,ipw2,1)=dielmat(1,ipw1,1,ipw2,1)
       dielmat(2,ipw1,2,ipw2,1)=dielmat(2,ipw1,1,ipw2,1)
       dielmat(1,ipw1,2,ipw2,2)=dielmat(1,ipw1,1,ipw2,1)
       dielmat(2,ipw1,2,ipw2,2)=dielmat(2,ipw1,1,ipw2,1)
     end do
   end do
 end if

!DEBUG
!write(std_out,*)' dielmt : make dielmat equal to identity matrix '
!do ipw1=1,npwdiel
!do ipw2=1,npwdiel
!dielmat(1,ipw1,1,ipw2,1)=0.0_dp
!dielmat(2,ipw1,1,ipw2,1)=0.0_dp
!end do
!end do
!ENDDEBUG

!Add the diagonal part
 do isp=1,min(nspden,2)
   do ipw=1,npwdiel
     dielmat(1,ipw,isp,ipw,isp)=one+dielmat(1,ipw,isp,ipw,isp)
   end do
 end do

!-The hermitian dielectric matrix is computed ------------------------------
!-Now, diagonalize it ------------------------------------------------------

!In RPA, everything is projected on the spin-symmetrized
!space. This was coded here (for the time being).

!Diagonalize the hermitian dielectric matrix

!npwsp=npwdiel*nspden
 npwsp=npwdiel

 ABI_ALLOCATE(dielh,(npwsp*(npwsp+1)))
 ABI_ALLOCATE(dielvec,(2,npwsp,npwsp))
 ABI_ALLOCATE(eig_diel,(npwsp))
 ABI_ALLOCATE(zhpev1,(2,2*npwsp-1))
 ABI_ALLOCATE(zhpev2,(3*npwsp-2))
 ier=0
!Store the dielectric matrix in proper mode before calling zhpev
 index=1
 do ii=1,npwdiel
   do jj=1,ii
     dielh(index  )=dielmat(1,jj,1,ii,1)
     dielh(index+1)=dielmat(2,jj,1,ii,1)
     index=index+2
   end do
 end do
!If spin-polarized and non RPA, need to store other parts of the matrix
!if(nspden/=1)then
!do ii=1,npwdiel
!Here, spin-flip contribution
!do jj=1,npwdiel
!dielh(index  )=dielmat(1,jj,1,ii,2)
!dielh(index+1)=dielmat(2,jj,1,ii,2)
!index=index+2
!end do
!Here spin down-spin down upper matrix
!do jj=1,ii
!dielh(index  )=dielmat(1,jj,2,ii,2)
!dielh(index+1)=dielmat(2,jj,2,ii,2)
!index=index+2
!end do
!end do
!end if

 call ZHPEV ('V','U',npwsp,dielh,eig_diel,dielvec,npwdiel,zhpev1,&
& zhpev2,ier)
 ABI_DEALLOCATE(zhpev1)
 ABI_DEALLOCATE(zhpev2)

 if(prtvol>=10)then
   write(message, '(a,a,a,5es12.4)' )ch10,&
&   ' Five largest eigenvalues of the hermitian RPA dielectric matrix:',&
&   ch10,eig_diel(npwdiel:npwdiel-4:-1)
   call wrtout(ab_out,message,'COLL')
 end if

 write(message, '(a,a)' )ch10,&
& ' dielmt : 15 largest eigenvalues of the hermitian RPA dielectric matrix'
 call wrtout(std_out,message,'COLL')
 write(message, '(a,5es12.5)' )'  1-5  :',eig_diel(npwdiel:npwdiel-4:-1)
 call wrtout(std_out,message,'COLL')
 write(message, '(a,5es12.5)' )'  6-10 :',eig_diel(npwdiel-5:npwdiel-9:-1)
 call wrtout(std_out,message,'COLL')
 write(message, '(a,5es12.5)' )'  11-15:',eig_diel(npwdiel-10:npwdiel-14:-1)
 call wrtout(std_out,message,'COLL')
 write(message, '(a,a)' )ch10,&
& ' dielmt : 5 smallest eigenvalues of the hermitian RPA dielectric matrix'
 call wrtout(std_out,message,'COLL')
 write(message, '(a,5es12.5)' )'  1-5  :',eig_diel(1:5)
 call wrtout(std_out,message,'COLL')

!Invert the hermitian dielectric matrix,
!Should use a BLAS call !
 do ipw2=1,npwdiel
   do ipw1=ipw2,npwdiel
     dielinv(1,ipw1,1,ipw2,1)=zero
     dielinv(2,ipw1,1,ipw2,1)=zero
   end do
 end do
 do ieig=1,npwdiel
   eiginv=one/eig_diel(ieig)
   do ipw2=1,npwdiel
     do ipw1=ipw2,npwdiel
       ar1=dielvec(1,ipw1,ieig)
       ai1=dielvec(2,ipw1,ieig)
       ar2=dielvec(1,ipw2,ieig)
       ai2=dielvec(2,ipw2,ieig)
       dielinv(1,ipw1,1,ipw2,1)=dielinv(1,ipw1,1,ipw2,1)+&
&       (ar1*ar2+ai1*ai2)*eiginv
       dielinv(2,ipw1,1,ipw2,1)=dielinv(2,ipw1,1,ipw2,1)+&
&       (ai1*ar2-ar1*ai2)*eiginv
     end do
   end do
 end do
 do ipw2=1,npwdiel-1
   do ipw1=ipw2+1,npwdiel
     dielinv(1,ipw2,1,ipw1,1)= dielinv(1,ipw1,1,ipw2,1)
     dielinv(2,ipw2,1,ipw1,1)=-dielinv(2,ipw1,1,ipw2,1)
   end do
 end do

 ABI_DEALLOCATE(dielh)
 ABI_DEALLOCATE(dielvec)
 ABI_DEALLOCATE(eig_diel)

!DEBUG
!Checks whether the inverse of the hermitian dielectric matrix
!has been correctly generated
!do ipw1=1,npwdiel
!do ipw2=1,npwdiel
!elementr=0.0_dp
!elementi=0.0_dp
!do ipw3=1,npwdiel
!elementr=elementr+dielinv(1,ipw1,1,ipw3,1)*dielmat(1,ipw3,1,ipw2,1)&
!&                    -dielinv(2,ipw1,1,ipw3,1)*dielmat(2,ipw3,1,ipw2,1)
!elementi=elementi+dielinv(1,ipw1,1,ipw3,1)*dielmat(2,ipw3,1,ipw2,1)&
!&                    +dielinv(2,ipw1,1,ipw3,1)*dielmat(1,ipw3,1,ipw2,1)
!end do
!if(elementr**2+elementi**2 > 1.0d-12)then
!if( ipw1 /= ipw2 .or. &
!&        ( abs(elementr-1.0_dp)>1.0d-6 .or. abs(elementi)>1.0d-6 ))then
!write(std_out,*)' dielmt : the inversion procedure is not correct '
!write(std_out,*)' ipw1, ipw2 =',ipw1,ipw2
!write(std_out,*)' elementr,elementi=',elementr,elementi
!stop
!end if
!end if
!end do
!end do
!write(std_out,*)'dielmt : matrix has been inverted successfully '
!stop
!ENDDEBUG

!Then get the inverse of the asymmetric
!dielectric matrix, as required for the preconditioning.

!Inverse of the dielectric matrix : ( 1 - 4pi (1/G^2) chi0(G,Gp) )^(-1)
!In dielinv there is now (1 - 4pi (1/G) chi0(G,Gp) (1/Gp) )^(-1)
!So, evaluate dielinv_after(G,Gp) =
!(4pi/G^2)^(1/2) dielinv_before(G,Gp) (4pi/Gp^2)^(-1/2)
!In RPA, can focus on the spin-averaged quantities
 do ipw1=1,npwdiel
   gred1=dble(kg_diel(1,ipw1))
   gred2=dble(kg_diel(2,ipw1))
   gred3=dble(kg_diel(3,ipw1))
   gsquar=tpisq*(gmet(1,1)*gred1**2+gmet(2,2)*gred2**2+gmet(3,3)*gred3**2 &
&   +two*( (gmet(1,2)*gred2+gmet(1,3)*gred3)* gred1 +      &
&   gmet(2,3)*gred2*gred3)                        )
!  Distinguish G=0 from other elements
   if(gsquar>tol12)then
     gfact=sqrt(four_pi/gsquar)
     gfactinv=one/gfact
     do ipw2=1,npwdiel
!      Must multiply both rows and columns
       dielinv(1,ipw2,1,ipw1,1)=dielinv(1,ipw2,1,ipw1,1)*gfactinv
       dielinv(2,ipw2,1,ipw1,1)=dielinv(2,ipw2,1,ipw1,1)*gfactinv
       dielinv(1,ipw1,1,ipw2,1)=dielinv(1,ipw1,1,ipw2,1)*gfact
       dielinv(2,ipw1,1,ipw2,1)=dielinv(2,ipw1,1,ipw2,1)*gfact
     end do
   else
!    Zero the G=0 elements, head
     do ipw2=1,npwdiel
       if (ipw2/=ipw1) dielinv(1:2,ipw1,1,ipw2,1)=zero
     end do
   end if
 end do

 ABI_DEALLOCATE(dielmat)

 call timab(90,2,tsec)

end subroutine dielmt
!!***


!!****f* ABINIT/dieltcel
!! NAME
!! dieltcel
!!
!! FUNCTION
!! Compute either test charge or electronic dielectric matrices
!! from susceptibility matrix
!! Diagonalize it, then invert it.
!!
!! INPUTS
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kxc(nfft,nkxc)=exchange-correlation kernel,
!!       needed if the electronic dielectric matrix is computed
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkxc=second dimension of the array kxc, see rhotoxc.f for a description
!!  npwdiel=size of the dielinv and susmat arrays.
!!  nspden=number of spin-density components
!!  occopt=option for occupancies
!!  option=1 for Test Charge dielectric matrix, 2 for electronic dielectric matrix
!!  prtvol=control print volume and debugging output
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!
!! OUTPUT
!!  dielinv(2,npwdiel,nspden,npwdiel,nspden)=inverse of the (non-hermitian)
!!      TC dielectric matrix in reciprocal space.
!!
!! NOTES
!! Output (not cleaned)
!! !!! Spin behaviour is not obvious !!!
!! Will not work in the spin-polarized, metallic case.
!!
!! PARENTS
!!      prcref,prcref_PMA
!!
!! CHILDREN
!!      destroy_mpi_enreg,fourdp,init_distribfft_seq,initmpi_seq,timab,wrtout
!!      zhpev
!!
!! SOURCE

subroutine dieltcel(dielinv,gmet,kg_diel,kxc,nfft,ngfft,nkxc,npwdiel,nspden,occopt,option,prtvol,susmat)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nkxc,npwdiel,nspden,occopt,option
 integer,intent(in) :: prtvol
!arrays
 integer,intent(in) :: kg_diel(3,npwdiel),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),kxc(nfft,nkxc)
 real(dp),intent(in) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 real(dp),intent(out) :: dielinv(2,npwdiel,nspden,npwdiel,nspden)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,ieig,ier,ifft,ii,index,ipw0,ipw1,ipw2,ispden,j1
 integer :: j2,j3,jj,k1,k2,k3,n1,n2,n3
 real(dp) :: ai,ai2,ar,ar2,eiginv,gred1,gred2,gred3,gsquar,si
 real(dp) :: sr,tpisq
 character(len=500) :: message
 type(MPI_type) :: mpi_enreg_seq
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: eig_msusinvsqr(:),eig_msussqr(:)
 real(dp),allocatable :: eig_sus(:),eig_sym(:),invsqrsus(:,:,:,:,:)
 real(dp),allocatable :: khxc(:,:,:,:,:),kxcg(:,:),sqrsus(:,:,:,:,:),sush(:)
 real(dp),allocatable :: susvec(:,:,:),symdielmat(:,:,:,:,:),symh(:)
 real(dp),allocatable :: symvec(:,:,:,:,:),wkxc(:),work(:,:,:,:,:)
 real(dp),allocatable :: work2(:,:,:,:,:),zhpev1(:,:),zhpev2(:)
!no_abirules
!integer :: ipw3
!real(dp) :: elementi,elementr
!DEBUG
!Used to moderate divergence effect near rho=0 of the Kxc
!this limit value is truly empirical (exprmt on small Sr cell).
!real(dp) :: kxc_min=-200.0
!ENDDEBUG

! *************************************************************************

 call timab(96,1,tsec)

!tpisq is (2 Pi) **2:
 tpisq=(two_pi)**2

 if(nspden/=1 .and. (occopt>=3 .and. occopt<=8) )then
   write(message, '(a,a,a)' )&
&   'In the present version of the code, one cannot produce',ch10,&
&   'the dielectric matrix in the metallic, spin-polarized case.'
   MSG_BUG(message)
 end if

 if(nspden==4)then
   write(message,'(a,a,a)')&
&   'In the present version of the code, one cannot produce',ch10,&
&   'the dielectric matrix in the non-collinear spin-polarized case.'
   MSG_ERROR(message)
 end if


!-Diagonalize the susceptibility matrix

 ABI_ALLOCATE(sush,(npwdiel*(npwdiel+1)))
 ABI_ALLOCATE(susvec,(2,npwdiel,npwdiel))
 ABI_ALLOCATE(eig_msusinvsqr,(npwdiel))
 ABI_ALLOCATE(eig_msussqr,(npwdiel))
 ABI_ALLOCATE(eig_sus,(npwdiel))
 ABI_ALLOCATE(zhpev1,(2,2*npwdiel-1))
 ABI_ALLOCATE(zhpev2,(3*npwdiel-2))
 ABI_ALLOCATE(work,(2,npwdiel,nspden,npwdiel,nspden))
 ABI_ALLOCATE(work2,(2,npwdiel,nspden,npwdiel,nspden))
 ABI_ALLOCATE(sqrsus,(2,npwdiel,nspden,npwdiel,nspden))
 ABI_ALLOCATE(invsqrsus,(2,npwdiel,nspden,npwdiel,nspden))

!At some time, should take care of different spin channels
 do ispden=1,nspden

   if(nspden/=1)then
     MSG_ERROR('dieltcel : stop, nspden/=1')
   end if

!  Store the susceptibility matrix in proper mode before calling zhpev
   index=1
   do ii=1,npwdiel
     do jj=1,ii
       sush(index  )=susmat(1,jj,1,ii,1)
       sush(index+1)=susmat(2,jj,1,ii,1)
       index=index+2
     end do
   end do

   ier=0
   call ZHPEV ('V','U',npwdiel,sush,eig_sus,susvec,npwdiel,zhpev1,zhpev2,ier)

!  DEBUG
!  write(std_out,*)' dieltcel : print eigenvalues of the susceptibility matrix'
!  do ii=1,npwdiel
!  write(std_out,'(i5,es16.6)' )ii,eig_sus(ii)
!  end do
!  ENDDEBUG

   do ii=1,npwdiel
     if(-eig_sus(ii)>1.d-12)then
       eig_msussqr(ii)=sqrt(-eig_sus(ii))
       eig_msusinvsqr(ii)=1._dp/eig_msussqr(ii)
     else if(-eig_sus(ii)< -1.d-12)then
       message = "Found positive eigenvalue of susceptibility matrix."
       MSG_BUG(message)
     else
!      Set the eigenvalue corresponding to a constant potential change to 1,
!      while it will be set to zero in Khx.
       eig_msussqr(ii)=1._dp
       eig_msusinvsqr(ii)=1._dp
     end if
   end do

!  Compute square root of minus susceptibility matrix
!  and inverse square root of minus susceptibility matrix
   do ii=1,npwdiel
     work(:,:,1,ii,1)=susvec(:,:,ii)*eig_msussqr(ii)
     work2(:,:,1,ii,1)=susvec(:,:,ii)*eig_msusinvsqr(ii)
   end do
   do ipw2=1,npwdiel
     do ipw1=ipw2,npwdiel
       ar=0._dp ; ai=0._dp ; ar2=0._dp ; ai2=0._dp
       do ii=1,npwdiel
         sr=susvec(1,ipw2,ii) ; si=susvec(2,ipw2,ii)
         ar =ar  +work(1,ipw1,1,ii,1)*sr  +work(2,ipw1,1,ii,1)*si
         ai =ai  +work(2,ipw1,1,ii,1)*sr  -work(1,ipw1,1,ii,1)*si
         ar2=ar2 +work2(1,ipw1,1,ii,1)*sr +work2(2,ipw1,1,ii,1)*si
         ai2=ai2 +work2(2,ipw1,1,ii,1)*sr -work2(1,ipw1,1,ii,1)*si
       end do
       sqrsus(1,ipw1,1,ipw2,1)=ar
       sqrsus(2,ipw1,1,ipw2,1)=ai
       invsqrsus(1,ipw1,1,ipw2,1)=ar2
       invsqrsus(2,ipw1,1,ipw2,1)=ai2
       if(ipw1/=ipw2)then
         sqrsus(1,ipw2,1,ipw1,1)=ar
         sqrsus(2,ipw2,1,ipw1,1)=-ai
         invsqrsus(1,ipw2,1,ipw1,1)=ar2
         invsqrsus(2,ipw2,1,ipw1,1)=-ai2
       end if
     end do
   end do

!  DEBUG
!  Checks whether sqrsus and invsqrsus are inverse of each other.
!  do ipw1=1,npwdiel
!  do ipw2=1,npwdiel
!  elementr=0.0_dp
!  elementi=0.0_dp
!  do ipw3=1,npwdiel
!  elementr=elementr+sqrsus(1,ipw1,1,ipw3,1)*invsqrsus(1,ipw3,1,ipw2,1)&
!  &                    -sqrsus(2,ipw1,1,ipw3,1)*invsqrsus(2,ipw3,1,ipw2,1)
!  elementi=elementi+sqrsus(1,ipw1,1,ipw3,1)*invsqrsus(2,ipw3,1,ipw2,1)&
!  &                    +sqrsus(2,ipw1,1,ipw3,1)*invsqrsus(1,ipw3,1,ipw2,1)
!  end do
!  if(elementr**2+elementi**2 > 1.0d-12)then
!  if( ipw1 /= ipw2 .or. &
!  &        ( abs(elementr-1.0_dp)>1.0d-6 .or. abs(elementi)>1.0d-6 ))then
!  write(std_out,*)' dieltcel : sqrsus and invsqrsus are not (pseudo)',&
!  &        'inverse of each other'
!  write(std_out,*)' ipw1, ipw2 =',ipw1,ipw2
!  write(std_out,*)' elementr,elementi=',elementr,elementi
!  stop
!  end if
!  end if
!  end do
!  end do
!  ENDDEBUG

!  End loop over spins
 end do

 ABI_DEALLOCATE(eig_msusinvsqr)
 ABI_DEALLOCATE(eig_msussqr)
 ABI_DEALLOCATE(eig_sus)
 ABI_DEALLOCATE(sush)
 ABI_DEALLOCATE(susvec)

!-Compute the Hxc kernel

 ABI_ALLOCATE(khxc,(2,npwdiel,nspden,npwdiel,nspden))
 ABI_ALLOCATE(symdielmat,(2,npwdiel,nspden,npwdiel,nspden))

 khxc(:,:,:,:,:)=0.0_dp

!Compute Hartree kernel
 do ipw1=1,npwdiel
   gred1=dble(kg_diel(1,ipw1))
   gred2=dble(kg_diel(2,ipw1))
   gred3=dble(kg_diel(3,ipw1))
   gsquar=tpisq*(gmet(1,1)*gred1**2+gmet(2,2)*gred2**2+gmet(3,3)*gred3**2 &
&   +2.0_dp*( (gmet(1,2)*gred2+gmet(1,3)*gred3)* gred1 +      &
&   gmet(2,3)*gred2*gred3)                        )
!  Distinguish G=0 from other elements
   if(gsquar>1.0d-12)then
     khxc(1,ipw1,1,ipw1,1)= 4.0_dp*pi/gsquar
   else
!    G=0
     ipw0=ipw1
   end if
 end do

!Eventually add the xc part
 if(option>=2)then

   ABI_ALLOCATE(wkxc,(nfft))
   ABI_ALLOCATE(kxcg,(2,nfft))
   wkxc(:)=kxc(:,1)
!  DEBUG
!  Used to moderate divergenc effect near rho=0 of the Kxc (see above).
!  wkxc(:)=merge(kxc(:,1), kxc_min, kxc(:,1) > kxc_min)
!  ENDDEBUG
   call initmpi_seq(mpi_enreg_seq)
   call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfft(2),ngfft(3),'all')
   call fourdp(1,kxcg,wkxc,-1,mpi_enreg_seq,nfft,1,ngfft,0) ! trsfrm R to G
   call destroy_mpi_enreg(mpi_enreg_seq)

!  Compute difference in G vectors
   n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
   do ipw2=1,npwdiel
     if(ipw2/=ipw0)then

       j1=kg_diel(1,ipw2) ; j2=kg_diel(2,ipw2) ; j3=kg_diel(3,ipw2)
!      Fills diagonal
       khxc(1,ipw2,1,ipw2,1)=khxc(1,ipw2,1,ipw2,1)+kxcg(1,1)
       khxc(2,ipw2,1,ipw2,1)=khxc(2,ipw2,1,ipw2,1)+kxcg(2,1)

       if(ipw2/=npwdiel)then
!        Fills off-diagonal part of the matrix, except G=0
         do ipw1=ipw2+1,npwdiel
           if(ipw1/=ipw0)then
             i1=kg_diel(1,ipw1) ; i2=kg_diel(2,ipw1) ; i3=kg_diel(3,ipw1)
!            Use of two mod calls handles both i1-j1>=ndiel1 AND i1-j1<0
             k1=mod(n1+mod(i1-j1,n1),n1)
             k2=mod(n2+mod(i2-j2,n2),n2)
             k3=mod(n3+mod(i3-j3,n3),n3)
             ifft=k1+1+n1*(k2+n2*k3)
!            The signs of imaginary contributions have been checked
             khxc(1,ipw1,1,ipw2,1)=kxcg(1,ifft)
             khxc(2,ipw1,1,ipw2,1)=kxcg(2,ifft)
             khxc(1,ipw2,1,ipw1,1)=kxcg(1,ifft)
             khxc(2,ipw2,1,ipw1,1)=-kxcg(2,ifft)
           end if
         end do
       end if

     end if
   end do

   ABI_DEALLOCATE(wkxc)
   ABI_DEALLOCATE(kxcg)

!  Endif option 2
 end if

!Now, get the symmetric dielectric matrix
!Premultiplication by square root of minus susceptibility matrix
 do ipw2=1,npwdiel
   do ipw1=1,npwdiel
     ar=0._dp ; ai=0._dp
     do ii=1,npwdiel
       ar=ar+sqrsus(1,ipw1,1,ii,1)*khxc(1,ii,1,ipw2,1) &
&       -sqrsus(2,ipw1,1,ii,1)*khxc(2,ii,1,ipw2,1)
       ai=ai+sqrsus(2,ipw1,1,ii,1)*khxc(1,ii,1,ipw2,1) &
&       +sqrsus(1,ipw1,1,ii,1)*khxc(2,ii,1,ipw2,1)
     end do
     work(1,ipw1,1,ipw2,1)=ar
     work(2,ipw1,1,ipw2,1)=ai
   end do
 end do
!Postmultiplication by square root of minus susceptibility matrix
 do ipw2=1,npwdiel
!  do ipw1=ipw2,npwdiel
   do ipw1=1,npwdiel
     ar=0._dp ; ai=0._dp
     do ii=1,npwdiel
       ar=ar+work(1,ipw1,1,ii,1)*sqrsus(1,ii,1,ipw2,1) &
&       -work(2,ipw1,1,ii,1)*sqrsus(2,ii,1,ipw2,1)
       ai=ai+work(2,ipw1,1,ii,1)*sqrsus(1,ii,1,ipw2,1) &
&       +work(1,ipw1,1,ii,1)*sqrsus(2,ii,1,ipw2,1)
     end do
     symdielmat(1,ipw1,1,ipw2,1)=ar
     symdielmat(2,ipw1,1,ipw2,1)=ai
!    if(ipw1/=ipw2)then
!    symdielmat(1,ipw2,1,ipw1,1)=ar
!    symdielmat(2,ipw2,1,ipw1,1)=-ai
!    end if
   end do
!  Add the unity matrix
   symdielmat(1,ipw2,1,ipw2,1)=1._dp+symdielmat(1,ipw2,1,ipw2,1)
 end do

 ABI_DEALLOCATE(khxc)

 ABI_ALLOCATE(symh,(npwdiel*(npwdiel+1)))
 ABI_ALLOCATE(symvec,(2,npwdiel,nspden,npwdiel,nspden))
 ABI_ALLOCATE(eig_sym,(npwdiel))

!Store the symmetrized dielectric matrix in proper mode before calling zhpev
 index=1
 do ii=1,npwdiel
   do jj=1,ii
     symh(index  )=symdielmat(1,jj,1,ii,1)
     symh(index+1)=symdielmat(2,jj,1,ii,1)
     index=index+2
   end do
 end do

 ier=0
 call ZHPEV ('V','U',npwdiel,symh,eig_sym,symvec,npwdiel,zhpev1,&
& zhpev2,ier)

 if(prtvol>=10)then
   write(message, '(a,a,a,5es12.4)' )ch10,&
&   ' Five largest eigenvalues of the symmetrized dielectric matrix:',&
&   ch10,eig_sym(npwdiel:npwdiel-4:-1)
   call wrtout(ab_out,message,'COLL')
 end if

 write(message,'(2a)')ch10,' dieltcel : 15 largest eigenvalues of the symmetrized dielectric matrix'
 call wrtout(std_out,message,'COLL')
 write(message, '(a,5es12.5)' )'  1-5  :',eig_sym(npwdiel:npwdiel-4:-1)
 call wrtout(std_out,message,'COLL')
 write(message, '(a,5es12.5)' )'  6-10 :',eig_sym(npwdiel-5:npwdiel-9:-1)
 call wrtout(std_out,message,'COLL')
 write(message, '(a,5es12.5)' )'  11-15:',eig_sym(npwdiel-10:npwdiel-14:-1)
 call wrtout(std_out,message,'COLL')
 write(message, '(2a)' )ch10,' dieltcel : 5 smallest eigenvalues of the symmetrized dielectric matrix'
 call wrtout(std_out,message,'COLL')
 write(message, '(a,5es12.5)' )'  1-5  :',eig_sym(1:5)
 call wrtout(std_out,message,'COLL')

!Invert the hermitian dielectric matrix,
 work(:,:,:,:,:)=0.0_dp
 do ieig=1,npwdiel
   eiginv=1.0_dp/eig_sym(ieig)
   do ipw2=1,npwdiel
!    do ipw1=ipw2,npwdiel
     do ipw1=1,npwdiel
       work(1,ipw1,1,ipw2,1)=work(1,ipw1,1,ipw2,1)+&
&       (symvec(1,ipw1,1,ieig,1)*symvec(1,ipw2,1,ieig,1)+ &
&       symvec(2,ipw1,1,ieig,1)*symvec(2,ipw2,1,ieig,1) ) * eiginv
       work(2,ipw1,1,ipw2,1)=work(2,ipw1,1,ipw2,1)+&
&       (symvec(2,ipw1,1,ieig,1)*symvec(1,ipw2,1,ieig,1)- &
&       symvec(1,ipw1,1,ieig,1)*symvec(2,ipw2,1,ieig,1) ) * eiginv
     end do
   end do
 end do
!if(npwdiel>1)then
!do ipw2=2,npwdiel
!do ipw1=1,ipw2-1
!work(1,ipw1,1,ipw2,1)= work(1,ipw2,1,ipw1,1)
!work(2,ipw1,1,ipw2,1)=-work(2,ipw2,1,ipw1,1)
!end do
!end do
!end if

 ABI_DEALLOCATE(eig_sym)
 ABI_DEALLOCATE(symh)
 ABI_DEALLOCATE(symvec)

!DEBUG
!Checks whether the inverse of the symmetric dielectric matrix
!has been correctly generated
!do ipw1=1,npwdiel
!do ipw2=1,npwdiel
!elementr=0.0_dp
!elementi=0.0_dp
!do ipw3=1,npwdiel
!elementr=elementr+work(1,ipw1,1,ipw3,1)*symdielmat(1,ipw3,1,ipw2,1)&
!&                    -work(2,ipw1,1,ipw3,1)*symdielmat(2,ipw3,1,ipw2,1)
!elementi=elementi+work(1,ipw1,1,ipw3,1)*symdielmat(2,ipw3,1,ipw2,1)&
!&                    +work(2,ipw1,1,ipw3,1)*symdielmat(1,ipw3,1,ipw2,1)
!end do
!if(elementr**2+elementi**2 > 1.0d-12)then
!if( ipw1 /= ipw2 .or. &
!&        ( abs(elementr-1.0_dp)>1.0d-6 .or. abs(elementi)>1.0d-6 ))then
!write(std_out,*)' dieltcel : the inversion procedure is not correct '
!write(std_out,*)' ipw1, ipw2 =',ipw1,ipw2
!write(std_out,*)' elementr,elementi=',elementr,elementi
!stop
!end if
!end if
!end do
!end do
!write(std_out,*)'dieltcel : matrix has been inverted successfully '
!ENDDEBUG

 ABI_DEALLOCATE(symdielmat)

!Then get the inverse of the asymmetric
!dielectric matrix, as required for the preconditioning.
!Premultiplication by square root of minus susceptibility matrix
 do ipw2=1,npwdiel
   do ipw1=1,npwdiel
     ar=0._dp ; ai=0._dp
     do ii=1,npwdiel
       ar=ar+invsqrsus(1,ipw1,1,ii,1)*work(1,ii,1,ipw2,1) &
&       -invsqrsus(2,ipw1,1,ii,1)*work(2,ii,1,ipw2,1)
       ai=ai+invsqrsus(2,ipw1,1,ii,1)*work(1,ii,1,ipw2,1) &
&       +invsqrsus(1,ipw1,1,ii,1)*work(2,ii,1,ipw2,1)
     end do
     work2(1,ipw1,1,ipw2,1)=ar
     work2(2,ipw1,1,ipw2,1)=ai
   end do
 end do
!Postmultiplication by square root of minus susceptibility matrix
 do ipw2=1,npwdiel
   do ipw1=1,npwdiel
     ar=0._dp ; ai=0._dp
     do ii=1,npwdiel
       ar=ar+work2(1,ipw1,1,ii,1)*sqrsus(1,ii,1,ipw2,1) &
&       -work2(2,ipw1,1,ii,1)*sqrsus(2,ii,1,ipw2,1)
       ai=ai+work2(2,ipw1,1,ii,1)*sqrsus(1,ii,1,ipw2,1) &
&       +work2(1,ipw1,1,ii,1)*sqrsus(2,ii,1,ipw2,1)
     end do
     dielinv(1,ipw1,1,ipw2,1)=ar
     dielinv(2,ipw1,1,ipw2,1)=ai
   end do
 end do

 ABI_DEALLOCATE(invsqrsus)
 ABI_DEALLOCATE(sqrsus)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(work2)
 ABI_DEALLOCATE(zhpev1)
 ABI_DEALLOCATE(zhpev2)

 call timab(96,2,tsec)

end subroutine dieltcel
!!***

!!****f* ABINIT/prcrskerker1
!! NAME
!! prcrskerker1
!!
!! FUNCTION
!! preconditionning by a real-space conjugate gradient on residual
!! using a model dielectric function in real space
!!
!! INPUTS
!!  nfft=number of fft grid points
!!  nspden=number of spin-density components
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  gprimd(3,3)=dimensional primitive translations in fourier space (bohr**-1)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  vresid(nfft,nspden)=residual potential
!!  base(nfft) = real space function used as a basis to guess a fine dielectric funtion
!!  see the calling routine to know the content
!!
!! OUTPUT
!!  vrespc(nfft,nspden)=preconditioned residual of the potential
!!
!! WARNINGS
!! This is experimental code : input, ouptput, results and any other feature may vary greatly.
!!
!! NOTES
!!  needs severe cleaning and this is abuse of modules as common blocks...
!!
!! PARENTS
!!      prcref,prcref_PMA
!!
!! CHILDREN
!!      cgpr,frskerker1__end,frskerker1__init,laplacian,prc_mem_init
!!
!! SOURCE

subroutine prcrskerker1(dtset,mpi_enreg,nfft,nspden,ngfft,dielar,etotal,gprimd,vresid,vrespc,base)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nspden
 real(dp) :: etotal
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: base(nfft),dielar(7),gprimd(3,3)
 real(dp),intent(in) :: vresid(nfft,nspden)
 real(dp),intent(out) :: vrespc(nfft,nspden)

!Local variables-------------------------------
!scalars
 integer ::  ifft,ispden,n1,n2,n3
 real(dp) :: base_delta,base_max,base_min,dielng,diemac,diemix
 real(dp) :: diemixmag
 real(dp) :: rdummy1,rdummy2
 logical ::  new_prc_func
!arrays
 real(dp) :: deltaW(nfft,nspden)
 real(dp) :: g2cart(nfft)
 real(dp) :: mat(nfft,nspden)

! *************************************************************************

!DEBUG
!write(std_out,*)' prckerker1 : enter '
!ENDDEBUG
!if(cycle==0) then
 call prc_mem_init(nfft)

 if(cycle==0) then
   new_prc_func=.TRUE.
   energy_min=etotal
 else if(etotal < energy_min) then
   new_prc_func=.TRUE.
   energy_min=etotal
 else
   new_prc_func=.FALSE.
 end if


 dielng=dielar(2)
 diemac=dielar(3)
 diemix=dielar(4)
 diemixmag=dielar(7)
!******************************************************************
!compute the diemac(r)                                          **
!******************************************************************
!this task will be devoted to a general function later
 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)
!base_cp=base
 if(new_prc_func) then
   base_min=base(1)
   base_max=base(1)
   do ifft=1,nfft
     base_min = min(base_min,base(ifft))
     base_max = max(base_max,base(ifft))
   end do
   base_delta = base_max - base_min
!  if(cycle.lt.2) then
   rdiemac(:) = (((base(:)-base_min) / (base_delta) ) *(diemac-one) + one)
!  else
!  rdiemac(:) = rdiemac(:)*0.5_dp+0.5_dp*(((base(:)-base_min) / (base_delta) ) *(diemac-one) + one)
!  end if
!  if(cycle==0) rdiemac(:) = (((base(:)-base_min) / (base_delta) ) *(diemac-one) + one)
!  rdiemac(:) = exp(((base(:)-base_min) / (base_delta) *log(diemac)))
 end if
 cycle=cycle+1
!if(cycle==5) cycle=0
!end if
!******************************************************************
!compute deltaW                                                 **
!******************************************************************
 vrespc=vresid !starting point
 ! put the laplacian of the residuals into deltaW
 call laplacian(gprimd,mpi_enreg,nfft,nspden,ngfft,rdfuncr=vrespc,laplacerdfuncr=deltaW,g2cart_out=g2cart) 

!call laplacian(vrespc,buffer,ngfft,gprimd) ! put the laplacian of the residuals into deltaW
!do ifft=1,nfft
!if (buffer(ifft,1)/=deltaW(ifft,1)) then
!stop
!end if
!end do
 deltaW(:,1)= diemix*(((one/rdiemac(:))*vresid(:,1))-(((dielng)**2)*deltaW(:,1)))
 if (nspden>1.and.(diemixmag>=zero)) then
   do ispden=2,nspden
     deltaW(:,ispden)= abs(diemixmag)*(((one/rdiemac(:))*vresid(:,ispden))-(((dielng)**2)*deltaW(:,ispden)))
   end do
 end if
!call random_number(deltaW)
!call random_number(vrespc)
!******************************************************************
!Finding the preconditionned residuals which minimizes          **
!half*(vrespc*(1-dielng2/4pi2 nabla2) vrespc) - vrespc * deltaW **
!***********************************************************************
 vrespc(:,1)=diemix*vrespc(:,1)
 if (nspden>1) vrespc(:,2:nspden)=abs(diemixmag)*vrespc(:,2:nspden)
!buffer=vrespc


!==============================================================================
!==============================================================================
!! Original loop
!==============================================================================
!==============================================================================

 call frskerker1__init(dtset,mpi_enreg,nfft,ngfft,nspden,dielng,deltaW,gprimd,mat,g2cart)

!call cgpr(pf_rscgres,dpf_rscgres,newvres,real(1e-40,dp),700,vrespc,rdummy1,rdummy2)
!rdummy1 = pf_rscgres(nfft,nspden,vrespc)
 call cgpr(nfft,nspden,frskerker1__pf,frskerker1__dpf,frskerker1__newvres,&
& real(1e-10,dp),700,vrespc,rdummy1,rdummy2)
 call frskerker1__end()

!==============================================================================
!==============================================================================
!! Original loop end
!==============================================================================
!==============================================================================


!cplex=1
!qphon(:)=zero
!call moddiel(cplex,dielar,nfft,ngfft,nspden,1,0,qphon,rprimd,vresid,buffer)
!c1=0
!do ifft=1,nfft,1
!if((abs(buffer(ifft,1)-vrespc(ifft,1))/(abs(buffer(ifft,1)+vrespc(ifft,1))*half)) > 5e-3) then
!c1=c1+1
!end if
!end do
!call laplacian(vrespc,buffer,ngfft,gprimd)
!buffer=vrespc(:,:)-buffer(:,:)*dielng**2
!c2=0
!do ifft=1,nfft,1
!if((abs(buffer(ifft,1)-deltaW(ifft,1))/(abs(buffer(ifft,1)+deltaW(ifft,1))*half)) > 5e-3) then
!c2=c2+1
!end if
!end do
!!!  !stop
!call laplacian(gprimd,mpi_enreg,nfft,nspden,ngfft,&
!& g2cart_out=g2cart)

!vrespc=vresid
!do ispden=1,nspden
!call fourdp(1, gvrespc(:,:,ispden), vrespc(:,ispden),-1,mpi_enreg,nfft,ngfft,0)
!end do
!filtering
!do ispden=1,nspden
!do ifft=1,nfft
!!    gvrespc(:,ifft,ispden)=(one-exp(-g2cart(ifft)*15.0_dp))*gvrespc(:,ifft,ispden)
!!      gvrespc(:,ifft,ispden)=(exp(-g2cart(ifft)*10.0_dp))*gvrespc(:,ifft,ispden)
!!      gvrespc(:,ifft,ispden)=(one-one/(exp(-0.002_dp/g2cart(ifft)**2)+one))*gvrespc(:,ifft,ispden)
!gvrespc(:,ifft,ispden)=(two-2_dp/(exp(-0.008_dp/(g2cart(ifft)+0.0012_dp))+one))*gvrespc(:,ifft,ispden)
!gvrespc(:,ifft,ispden)=min(one,(sqrt(g2cart(ifft)/0.006_dp))**(one))*gvrespc(:,ifft,ispden)
!end do
!end do
!change resulting potential to real space
!do ispden=1,nspden
!call fourdp(1,gvrespc(:,:,ispden),vrespc(:,ispden),1,mpi_enreg,nfft,ngfft,0)
!end do
!vrespc=vrespc*diemix
!maxg2=g2cart(1)
!ming2=g2cart(5)
!do ifft=1,nfft
!maxg2=max(g2cart(ifft),maxg2)
!if(g2cart(ifft) .gt. zero) ming2=min(g2cart(ifft),ming2)
!end do
!stop

!DEBUG
!write(std_out,*)' prckerker1 : exit '
!ENDDEBUG

end subroutine prcrskerker1
!!***

!!****f* ABINIT/prcrskerker2
!! NAME
!! prcrskerker2
!!
!! FUNCTION
!! preconditionning by a real-space conjugate gradient on residual
!! using a model dielectric function in real space
!! differing from prcrskerker1 by the
!! use of a linear response approach
!!
!! INPUTS
!!  nfft=number of fft grid points
!!  nspden=number of spin-density components
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  gprimd(3,3)=dimensional primitive translations in fourier space (bohr**-1)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  vresid(nfft,nspden)=residual potential
!!
!! OUTPUT
!!  vrespc(nfft,nspden)=preconditioned residual of the potential
!!
!! WARNINGS
!! This is experimental code : input, ouptput, results and any other feature may vary greatly.
!!
!! NOTES
!!
!! PARENTS
!!      prcref,prcref_PMA
!!
!! CHILDREN
!!      cgpr,dotprod_vn,frskerker2__end,frskerker2__init,laplacian,ptabs_fourdp
!!
!! SOURCE

subroutine prcrskerker2(dtset,nfft,nspden,ngfft,dielar,gprimd,rprimd,vresid,vrespc,natom,xred,mpi_enreg,ucvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: dielar(7),gprimd(3,3),rprimd(3,3),vresid(nfft,nspden)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(out) :: vrespc(nfft,nspden)

!Local variables-------------------------------
  !logical,save ::ok=.FALSE.
!scalars
 integer :: cplex,i1,i2,i3,iatom,iatom27,ifft,ispden,n1,n2,n3,natom27,nfftotf
 integer :: option
 real(dp),save :: lastp1=one,lastp2=one
 real(dp) :: C1,C2,DE,core,dielng,diemac,diemix,diemixmag,doti,dr,l1,l2,l3,l4,r
 real(dp) :: rdummy1,rdummy2,rmin,xr,y,yr,zr
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: V1(nfft,nspden),V2(nfft,nspden),buffer(nfft,nspden)
 real(dp) :: deltaW(nfft,nspden)
 real(dp) :: mat(nfft,nspden)
 real(dp) :: rdielng(nfft),rdiemac(nfft),xcart(3,natom)
 real(dp) :: xcart27(3,natom*27)

! *************************************************************************

 dielng=dielar(2)
 diemac=dielar(3)
 diemix=dielar(4)
 diemixmag=dielar(7)
!******************************************************************
!compute the diemac(r)                                          **
!******************************************************************
!this task will be devoted to a general function later
 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)
 nfftotf=n1*n2*n3
!if(.not.ok) then
 xcart(1,:)=xred(1,:)*rprimd(1,1)+xred(2,:)*rprimd(1,2)+xred(3,:)*rprimd(1,3)
 xcart(2,:)=xred(1,:)*rprimd(2,1)+xred(2,:)*rprimd(2,2)+xred(3,:)*rprimd(2,3)
 xcart(3,:)=xred(1,:)*rprimd(3,1)+xred(2,:)*rprimd(3,2)+xred(3,:)*rprimd(3,3)

 iatom27=0
 do i1=-1,1
   do i2=-1,1
     do i3=-1,1
       do iatom=1,natom
         iatom27=iatom27+1
         xcart27(:,iatom27)=xcart(:,iatom)+rprimd(:,1)*i1+rprimd(:,2)*i2+rprimd(:,3)*i3
       end do
     end do
   end do
 end do

!stop
 natom27=27*natom

 l1=0.34580850339844665
!l2=0.5123510203906797 !0.41242551019533985
!l3=0.8001489796093203 !0.90007448980466009

 l2=0.41242551019533985
 l3=0.90007448980466009
 l4=0.9666914966015534


 l1=0.31387233559896449
 l2=0.35828367346355994
 l3=0.9333829932031068
 l4=0.9777943310677023

 l1=3.5
 l2=11.5
 l3=2.5
 l4=6.5
!l1=30. !cellules pleines

 rdielng=zero
 core=1. !(value of Er at the core of atoms)
 dr=2.65 ! radius of atoms=2.65165
 y=1. ! value of Er in the empty region

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 do i3=1,n3
   ifft=(i3-1)*n1*(n2/mpi_enreg%nproc_fft)
   do i2=1,n2
     if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
       do i1=1,n1
         ifft=ifft+1
!        !!!!!!!!!!!!!!!!!!!!!!!!!
!        ! calculation of the simplest part void/metal
!        !!              x=real(real(i3,dp)/real(n3,dp),dp)
!        !!              !x=i3/n3
!        !!              if(x < l1) then
!        !!                 rdiemac(ifft)=diemac
!        !!                 rdielng(ifft)=dielng
!        !!              else if(x < l2) then
!        !!                 xp=(l2-x)/(l2-l1)
!        !!                 rdiemac(ifft)=y+(diemac-y)&
!        !!                      & * (1.-(1.-xp)**4)**4
!        !!                 rdielng(ifft)=dielng*(1.-(1.-xp)**4)**4
!        !!              else if(x < l3) then
!        !!                 rdiemac(ifft)=y
!        !!                 rdielng(ifft)=zero
!        !!              else if(x < l4) then
!        !!                 xp=(l3-x)/(l3-l4)
!        !!                 rdiemac(ifft)=y+(diemac-y)&
!        !!                      & * (1.-(1.-xp)**4)**4
!        !!                 rdielng(ifft)=dielng*(1.-(1.-xp)**4)**4
!        !!              else
!        !!                 rdiemac(ifft)=diemac
!        !!                 rdielng(ifft)=dielng
!        !!              end if
!        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        !!!! calculation of atomic core dielectric
!        !!              rmin=1e16
!        !!              xr=real(real((i1-1),dp)/n1,dp)*rprimd(1,1)+real(real((i2-1),dp)/n2,dp)*rprimd(1,2)&
!        !!                   &+real((i3-1),dp)/real(n3,dp)*rprimd(1,3)
!        !!              yr=real(real((i1-1),dp)/n1,dp)*rprimd(2,1)+real(real((i2-1),dp)/n2,dp)*rprimd(2,2)&
!        !!                   &+real((i3-1),dp)/real(n3,dp)*rprimd(2,3)
!        !!              zr=real(real((i1-1),dp)/n1,dp)*rprimd(3,1)+real(real((i2-1),dp)/n2,dp)*rprimd(3,2)&
!        !!                   &+real((i3-1),dp)/real(n3,dp)*rprimd(3,3)
!        !!              do iatom=1,natom27
!        !!                 r=(xr-xcart27(1,iatom))**2+(yr-xcart27(2,iatom))**2+(zr-xcart27(3,iatom))**2
!        !!                 if (r<rmin) then
!        !!                    rmin=r
!        !!                 end if
!        !!              end do
!        !!              if(rmin < dr**2) then
!        !!                 rdiemac(ifft)=min(rdiemac(ifft),core+(diemac-core)*(1.-(1.-sqrt(rmin)/dr)**2)**2)
!        !!                 rdielng(ifft)=dielng-dielng*(1.-(1.-sqrt(rmin)/dr)**4)**4
!        !!              else
!        !!                 rdiemac(ifft)=min(rdiemac(ifft),diemac)
!        !!              end if
         rmin=1e16
         xr=real(real((i1-1),dp)/n1,dp)*rprimd(1,1)+real(real((i2-1),dp)/n2,dp)*rprimd(1,2)&
&         +real((i3-1),dp)/real(n3,dp)*rprimd(1,3)
         yr=real(real((i1-1),dp)/n1,dp)*rprimd(2,1)+real(real((i2-1),dp)/n2,dp)*rprimd(2,2)&
&         +real((i3-1),dp)/real(n3,dp)*rprimd(2,3)
         zr=real(real((i1-1),dp)/n1,dp)*rprimd(3,1)+real(real((i2-1),dp)/n2,dp)*rprimd(3,2)&
&         +real((i3-1),dp)/real(n3,dp)*rprimd(3,3)

         rdiemac(ifft)=y
         rdielng(ifft)=zero
         do iatom=1,natom27
           r=(xr-xcart27(1,iatom))**2+(yr-xcart27(2,iatom))**2+(zr-xcart27(3,iatom))**2
           if (r<rmin) then
             rmin=r
           end if
           if(r < l1) then
             rdiemac(ifft)= rdiemac(ifft) +  0.7_dp * (diemac-y)
           else if(r < l2) then
             rdiemac(ifft)= rdiemac(ifft) + 0.7_dp * (diemac-y)*(one-((sqrt(r)-l1)/(l2-l1))**2)**2
           else
             rdiemac(ifft)=rdiemac(ifft)
           end if
           if(r < l3) then
             rdielng(ifft)= rdielng(ifft) +  0.5_dp * (dielng)
           else if(r < l4) then
             rdielng(ifft)= rdielng(ifft) + 0.5_dp * (dielng)  *(one-((sqrt(r)-l3)/(l4-l3))**2)**2
           end if
         end do

         rdielng(ifft)=min(rdielng(ifft),dielng)
!        rdielng(ifft)=dielng
         rdiemac(ifft)=min(rdiemac(ifft),diemac)
         rdiemac(ifft)=diemac
       end do
     end if
   end do
 end do
!rdielng(:)=dielng

!****************************************************************************************
!****************************************************************************************
!****************************************************************************************
!****************************************************************************************
!******************************************************************
!compute V1
!******************************************************************
 V1=vresid
 call laplacian(gprimd,mpi_enreg,nfft,nspden,ngfft,rdfuncr=V1,laplacerdfuncr=deltaW)
 deltaW(:,1)= (((one/rdiemac(:))*V1(:,1))-(((rdielng(:))**2)*deltaW(:,1)))
!deltaW(:,1)= -diemix*(((rdielng(:))**2)*deltaW(:,ispden))
 if (nspden>1) then
   do ispden=2,nspden
     deltaW(:,ispden)= (((one/rdiemac(:))*V1(:,ispden))-(((rdielng(:))**2)*deltaW(:,ispden)))
!    deltaW(:,ispden)= -abs(diemixmag)*(((rdielng(:))**2)*deltaW(:,ispden))
   end do
 end if
 call frskerker2__init(dtset,mpi_enreg,nfft,ngfft,nspden,rdielng,deltaW,gprimd,mat)
 call cgpr(nfft,nspden,frskerker2__pf,frskerker2__dpf,&
& frskerker2__newvres2,lastp1*real(1e-6 ,dp),700,V1,rdummy1,rdummy2)
 lastp1=min(abs(rdummy1),1e-6_dp)
 call frskerker2__end()

!******************************************************************
!compute V2
!******************************************************************
 V2=vresid
 do ispden=1,nspden
   deltaW(:,ispden)= (rdielng(:)**2)
 end do
 call frskerker2__init(dtset,mpi_enreg,nfft,ngfft,nspden,rdielng,deltaW,gprimd,mat)
 call cgpr(nfft,nspden,frskerker2__pf,frskerker2__dpf,&
& frskerker2__newvres2,lastp2*real(1e-6,dp),700,V2,rdummy1,rdummy2)
 lastp2=min(abs(rdummy1),1e-6_dp)
 call frskerker2__end()


!******************************************************************
!compute C1, C2 & DE
!******************************************************************
 cplex=1;
 option=1;
 call dotprod_vn(cplex,& !complex density/pot
&rdielng,&          !the density
&DE,&  !resulting dorproduct integrated over r  ! here DE is used has a buffer
&doti,&          !imaginary part of the integral
&size(rdielng,1),&          !number of localy(cpu) attributed grid point
&nfftotf,&        !real total number of grid point
&nspden,&        !nspden
&option,&        !1=compute only the real part 2=compute also the imaginary part
&rdielng,&          !the potential
&ucvol,&          !cell volume
&mpi_comm_sphgrid=mpi_enreg%comm_fft)
 do ispden=1,nspden
   buffer(:,ispden)=rdielng(:)*V1(:,ispden)
 end do
 call dotprod_vn(cplex,& !complex density/pot
&rdielng,&          !the density
&C1,&  !resulting dorproduct integrated over r  ! here DE is used has a buffer
&doti,&          !imaginary part of the integral
&size(rdielng,1),&          !number of localy(cpu) attributed grid point
&nfftotf,&        !real total number of grid point
&nspden,&        !nspden
&option,&        !1=compute only the real part 2=compute also the imaginary part
&buffer,&          !the potential
&ucvol,&         !cell volume
&mpi_comm_sphgrid=mpi_enreg%comm_fft)
 do ispden=1,nspden
   buffer(:,ispden)=rdielng(:)*V2(:,ispden)
 end do
 call dotprod_vn(cplex,& !complex density/pot
&rdielng,&          !the density
&C2,&  !resulting dorproduct integrated over r  ! here DE is used has a buffer
&doti,&          !imaginary part of the integral
&size(rdielng,1),&          !number of localy(cpu) attributed grid point
&nfftotf,&        !real total number of grid point
&nspden,&        !nspden
&option,&        !1=compute only the real part 2=compute also the imaginary part
&buffer,&          !the potential
&ucvol,&         !cell volume
&mpi_comm_sphgrid=mpi_enreg%comm_fft)
 C1=C1/DE
 C2=C2/DE
 DE=C1/(one-C2)

!******************************************************************
!compute the new preconditionned residuals
!******************************************************************
 vrespc(:,1)=diemix*(V1(:,1)+DE*V2(:,1))
 if (nspden>1) vrespc(:,2:nspden)=abs(diemixmag)*(V1(:,2:nspden)+DE*V2(:,2:nspden))

end subroutine prcrskerker2
!!***

!!****f* ABINIT/cgpr
!! NAME
!! cgpr
!!
!! FUNCTION
!! perform Polak-Ribiere conjugate gradient on a function f
!! implementation based on the cg recipe of "numerical recipe"
!!
!! INPUTS
!! dp_dum_vdp: function  to be minimized (return a dp from a vector of dp)
!! vdp_dum_vdp: derivative of f
!! dtol: precision precision required for the minimization
!! itmax: number of iterations allowed (each linmin will be done with at max 10 times
!! this number
!!
!! OUTPUT
!! fmin: value of f at the minimum
!! lastdelta: absolute value of the last delta between steps
!! SIDE EFFECTS
!! v: vector on which minimization is to be performed, starting point
!! and resulting min
!!
!! PARENTS
!!      prcrskerker1,prcrskerker2
!!
!! CHILDREN
!!      linmin
!!
!! SOURCE

subroutine cgpr(nv1,nv2,dp_dum_v2dp,v2dp_dum_v2dp,sub_dum_dp_v2dp_v2dp,dtol,itmax,v,fmin,delta)

!Arguments ------------------------------------
include "dummy_functions.inc"
!scalars
 integer,intent(in) :: itmax,nv1,nv2
 real(dp),intent(in) :: dtol
 real(dp),intent(out) :: delta,fmin
!arrays
 real(dp),intent(inout) :: v(nv1,nv2)

!Local variables-------------------------------
!scalars
 integer :: iiter
 real(dp) :: fv,gam,gscal,gscal2,sto
!arrays
 real(dp) :: grad0(nv1,nv2),grad1(nv1,nv2),grad2(nv1,nv2),grad3(nv1,nv2)
!no_abirules

!************************************************************************
 fv = dp_dum_v2dp(nv1,nv2,v(:,:))
 grad0(:,:) = -v2dp_dum_v2dp(nv1,nv2,v(:,:))
 grad1(:,:) = grad0(:,:)
 grad2(:,:) = grad0(:,:)
 do iiter=1,itmax
  call linmin(nv1,nv2,dp_dum_v2dp,v2dp_dum_v2dp,sub_dum_dp_v2dp_v2dp,v,grad0,fmin)
! return if the min is reached
  sto=dtol*(abs(fmin)+abs(fv)+tol14)
  delta=abs(fv-fmin)
  delta=abs(delta)
  if((delta.lt.sto).or.(iiter==itmax)) then
!  DEBUG
!  write(std_out,*) 'cgpr (01cg) : stop cond for cgpr:',sto,'delta:',delta,'fv:',fv,'fmin:',fmin
!  ENDDEBUG
   return
  end if
! a new step
  fv=fmin
  grad0(:,:)=v2dp_dum_v2dp(nv1,nv2,v(:,:))
  gscal=dotproduct(nv1,nv2,grad1(:,:),grad1(:,:))
  grad3(:,:)=grad0(:,:)+grad1(:,:)
  gscal2=dotproduct(nv1,nv2,grad3(:,:),grad0(:,:))
  gam=gscal2/gscal
  grad1(:,:)=-grad0(:,:)
  grad2(:,:)=grad1(:,:)+gam*grad2(:,:)
  grad0(:,:)=grad2(:,:)
! DEBUG
! write(std_out,*) 'cgpr (01cg) :================================================================================='
! write(std_out,*) 'cgpr (01cg) : step',iiter,'delta:',delta ,'fv',fv,'fmin',fmin
! write(std_out,*) 'cgpr (01cg) :================================================================================='
! ENDDEBUG
 end do

end subroutine cgpr
!!***

!!****f* ABINIT/linmin
!! NAME
!! linmin
!!
!! FUNCTION
!! minimizes a function along a gradient line:
!! first bracket the minimum then perform the minimization
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!! dp_dum_vdp: function  to be minimized (return a dp from a vector of dp)
!! vdp_dum_vdp: derivative of f
!!
!! OUTPUT
!! fmin: minimun value reached for dp_dum_vdp
!!
!! SIDE EFFECTS
!! grad: the gradient line along which the minimization is performed (not changed)
!! v: the starting and then ending point of the minimization
!!
!! PARENTS
!!      cgpr
!!
!! CHILDREN
!!      bracketing
!!
!! SOURCE

subroutine linmin(nv1,nv2,dp_dum_v2dp,v2dp_dum_v2dp,sub_dum_dp_v2dp_v2dp,v,grad,fmin)

!Arguments ------------------------------------
include "dummy_functions.inc"
!scalars
 integer,intent(in) :: nv1,nv2
 real(dp),intent(out) :: fmin
!arrays
 real(dp),intent(inout) :: grad(nv1,nv2),v(nv1,nv2)

!Local variables-------------------------------
!scalars
 real(dp),parameter :: maglimit=10000.0_dp,tol=tol8*tol8*tol8
 real(dp) :: a,b,fa,fb,fx,x,xmin
!no_abirules

!************************************************************************
 a=zero
 x=ninth*real(1e-4,dp)
 call bracketing (nv1,nv2,dp_dum_v2dp,v,grad,a,x,b,fa,fx,fb)
!DEBUG
!write(std_out,*) 'linmin (01cg) : linmin after bracketing'
!write(std_out,*) 'linmin (01cg) : point',a,x,b,'value',fa,fx,fb
!ENDDEBUG
 fmin =brent(nv1,nv2,dp_dum_v2dp,v2dp_dum_v2dp,sub_dum_dp_v2dp_v2dp,6,v,grad,a,x,b,tol,xmin)

end subroutine linmin
!!***

!!****f* ABINIT/bracketing
!! NAME
!! bracketing
!!
!! FUNCTION
!! bracket a minimun of a function f
!!
!! INPUTS
!! dp_dum_vdp: the function of which the mimimum should be bracketted
!!
!!
!! OUTPUT
!! b= last member of the bracketing triplet a < x < b
!! fa,fx,fb= value of the function at dp_dum_vdp(v(:)+y*grad(:))
!!
!!
!! SIDE EFFECTS
!! v: the initial vector for the function (return unchanged)
!! grad: the direction on which the bracketting is to be performed (return unchanged)
!! a,x: two members of the bracketing triplet (see b)
!!
!! PARENTS
!!      linmin
!!
!! CHILDREN
!!
!! SOURCE

subroutine bracketing (nv1,nv2,dp_dum_v2dp,v,grad,a,x,b,fa,fx,fb)

!Arguments ------------------------------------
include "dummy_functions.inc"
!scalars
 integer,intent(in) :: nv1,nv2
 real(dp),intent(inout) :: a,x
 real(dp),intent(out) :: b,fa,fb,fx
!arrays
 real(dp),intent(inout) :: grad(nv1,nv2),v(nv1,nv2)

!Local variables-------------------------------
!scalars
 real(dp),parameter :: maglimit=10000.0_dp
 real(dp) :: c,fu,q,r,u,ulim

! *************************************************************************

 fa=dp_dum_v2dp(nv1,nv2,v(:,:)+(a*grad(:,:)))
 fx=dp_dum_v2dp(nv1,nv2,(x*grad(:,:))+v(:,:))
 if(fx > fa) then
  c=a
  a=x
  x=c
  c=fa
  fa=fx
  fx=c
 end if
 b=x+gold*(x-a)
 fb=dp_dum_v2dp(nv1,nv2,(b*grad(:,:))+v(:,:))
 do
  if (fx <= fb) return
  r=(x-a)*(fx-fb)
  q=(x-b)*(fx-fa)
  u=x-((x-b)*q-(x-a)*r)/(two*sign(max(abs(q-r),smallest_real),q-r))
  ulim=x+maglimit*(b-x)
  if((x-u)*(u-b) > zero) then
   fu=dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:))
   if(fu < fb) then
    a=x
    fa=fx
    x=u
    fx=fu
    return
   else if (fx < fu) then
    b=u
    fb=fu
    return
   end if
   u=b+gold*(b-x)
   fu=dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:))
  else if((b-u)*(u-ulim) > zero) then
   fu=dp_dum_v2dp(nv1,nv2,u*grad(:,:)+v(:,:))
   if(fu<fb) then
    x=b
    b=u
    u=b+gold*(b-x)
    fx=fb
    fb=fu
    fu=dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:))
   end if
  else if((u-ulim)*(ulim-b) >= zero) then
   u=ulim
   fu=dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:))
  else
   u=b+gold*(b-x)
   fu=dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:))
  end if
  a=x
  x=b
  b=u
  fa=fx
  fx=fb
  fb=fu
 end do

end subroutine bracketing
!!***

!!****f* ABINIT/brent
!! NAME
!! brent
!!
!! FUNCTION
!! minimizes a function along a line
!!
!! INPUTS
!! dp_dum_vdp: function  to be minimized (return a dp from a vector of dp)
!! vdp_dum_vdp: derivative of the function (return a vector of dp from a vector of dp)
!! itmax: number of iterations allowed
!! tol: tolerance on error. It depend on the precision of the numbers
!! (usualy chosen as sqrt(max precision available with your floating point reresentation))
!! ax,xx,bx: a bracketing triplet around the minimum to be find
!! OUTPUT
!! xmin: value such that dp_dum_vdp(v(:)+xmin*grad(:)) is minimum
!! brent:  dp_dum_vdp(v(:)+xmin*grad(:))
!!
!! SIDE EFFECTS
!! grad(:): direction along which the minimization is performed
!! v(:): starting and ending point of the minimization
!!
!! PARENTS
!! linmin
!!
!! CHILDREN
!! dotproduct
!!
!! SOURCE

function brent(nv1,nv2,dp_dum_v2dp,v2dp_dum_v2dp,sub_dum_dp_v2dp_v2dp,itmax,v,grad,ax,xx,bx,tol,xmin)

!Arguments ------------------------------------
include "dummy_functions.inc"
!scalars
 integer,intent(in) :: itmax,nv1,nv2
 real(dp) :: brent
 real(dp),intent(in) :: ax,bx,tol,xx
 real(dp),intent(out) :: xmin
!arrays
 real(dp),intent(inout) :: grad(nv1,nv2),v(nv1,nv2)

!Local variables-------------------------------
!scalars
 integer :: iter
 real(dp) :: a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,vv,w
 real(dp) :: x,xm,zeps
 logical :: ok1,ok2,ok3,ok4

!************************************************************************
 zeps=epsilon(ax*real(1e-2,dp))
 a=min(ax,bx)
 b=max(ax,bx)
 vv=xx
 w=xx
 x=xx
 e=zero
 fx=dp_dum_v2dp(nv1,nv2,x*grad(:,:)+v(:,:))
 fv=fx
 fw=fx
!the function sub_dum_dp_v2dp_v2dp must do the equivalent of
!v(:,:)=v(:,:)+(grad(:,:)*x)
!but for instance renormilizing the density if brent is used on a density...
!vp(:,:) = v(:,:)
!sub_dum_dp_v2dp_v2dp(x,grad(:,:),vp(:,:)
!dx=dotproduct(v2dp_dum_v2dp(vp(:,:)),grad(:,:))
 dx=dotproduct(nv1,nv2,v2dp_dum_v2dp(nv1,nv2,v(:,:)+x*grad(:,:)),grad(:,:))
 dv=dx
 dw=dx
 do iter=1,itmax
  xm=half*(a+b)
  tol1=tol*abs(x)+zeps
  tol2=two*tol1
  if(abs(x-xm) <= (tol2-half*(b-a))) then
   exit
  end if
  if(abs(e) > tol1) then
   d1=two*(b-a)
   d2=d1
   if(dw /= dx) d1=(w-x)*dx/(dx-dw)
   if(dv /= dx) d2=(vv-x)*dx/(dx-dv)
   u1=x+d1
   u2=x+d2
   ok1=((a-u1)*(u1-b)>zero).and.(dx*d1<=zero)
   ok2=((a-u2)*(u2-b)>zero).and.(dx*d2<=zero)
   olde=e
   e=d
   if(ok1.or.ok2) then
    if(ok1.and.ok2) then
     d=merge(d1,d2,abs(d1)<abs(d2))
    else
     d=merge(d1,d2,ok1)
    end if
    if(abs(d)<=abs(half*olde)) then
     u=x+d
     if(((u-a)<tol2).or.((b-u)<tol2)) d=sign(tol1,xm-x)
    else
     e=merge(a,b,dx>=zero)-x
     d=half*e
    end if
   else
    e=merge(a,b,dx>=zero)-x
    d=half*e
   end if
  else
   e=merge(a,b,dx>=zero)-x
   d=half*e
  end if

  if(abs(d) >=tol1)then
   u=x+d
   fu=dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:))
  else
   u=x+sign(tol1,d)
   fu=dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:))
   if(fu>fx) then
    exit
   end if
  end if
  du=dotproduct(nv1,nv2,v2dp_dum_v2dp(nv1,nv2,(u*grad(:,:))+v(:,:)),grad(:,:))
  if(fu<=fx)then
   if(u>=x)then
    a=x
   else
    b=x
   end if
   vv=w
   fv=fw
   dv=dw
   w=x
   fw=fx
   dw=dx
   x=u
   dx=du
   fx=fu
  else
   if(u<x) then
    a=u
   else
    b=u
   end if
   ok3=(w==x).or.(fu.le.fw)
   ok4=(vv==w).or.(vv==x).or.(fu.lt.fv)
   if(ok3) then
    vv=w
    fv=fw
    dv=dw
    w=u
    fw=fu
    dw=du
   else if( ok4 ) then
    vv=u
    fv=fu
    dv=du
   end if
  end if
 end do
 xmin=x
!the function sub_dum_dp_v2dp_v2dp must do the equivalent of
!v(:,:)=v(:,:)+(grad(:,:)*x)
!but for instance renormilizing the density if brent is used on a density...
 call sub_dum_dp_v2dp_v2dp(nv1,nv2,x,grad(:,:),v(:,:))
 brent=fx

end function brent
!!***

end module m_prcref
!!***
