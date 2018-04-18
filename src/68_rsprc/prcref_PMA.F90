!{\src2tex{textfont=tt}}
!!****f* ABINIT/prcref_PMA
!!
!! NAME
!! prcref_PMA
!!
!! FUNCTION
!! Compute preconditioned residual potential (or density) and forces.
!! iprcel, densfor_pred and iprcfc govern the choice of the preconditioner.
!! Three tasks are done :
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
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA,XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

  subroutine prcref_PMA(atindx,dielar,dielinv,&
&  dielstrt,dtn_pc,dtset,fcart,ffttomix,gmet,gsqcut,&
&  istep,kg_diel,kxc,&
&  mgfft,moved_atm_inside,mpi_enreg,my_natom,&
&  nattyp,nfft,nfftprc,ngfft,ngfftprc,nkxc,npawmix,npwdiel,ntypat,n1xccc,&
&  optreal,optres,pawrhoij,ph1d,psps,rhog, rhoijrespc,rhor,rprimd,&
&  susmat,vhartr,vpsp,vresid,vrespc,vxc,xred,&
&  etotal,pawtab,wvl)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_errors
 use m_profiling_abi
 use m_xmpi
 use m_cgtools
 use m_xcdata

 use m_geometry, only : xcart2xred, metric
 use m_pawtab,   only : pawtab_type
 use m_pawrhoij, only : pawrhoij_type
 use m_fft,      only : zerosym
 use m_fftcore,  only : kgindex
 use m_kg,       only : getph
 use m_dtset,    only : testsusmat

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prcref_PMA'
 use interfaces_52_fft_mpi_noabirule
 use interfaces_53_ffts
 use interfaces_56_xc
 use interfaces_64_psp
 use interfaces_67_common
 use interfaces_68_rsprc, except_this_one => prcref_PMA
!End of the abilint section

 implicit none

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
 logical :: computediel
 logical :: non_magnetic_xc
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

! Initialise non_magnetic_xc for rhohxc
 non_magnetic_xc=(dtset%usepawu==4).or.(dtset%usepawu==14)

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
!to X. Gonze PRB vol54 nb7 p4383 (1996)
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
       call fourdp(1,work,work1(:,ispden),+1,mpi_enreg,nfftprc,ngfftprc,dtset%paral_kgb,0)
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
       call fourdp(1,rhog_wk,work,+1,mpi_enreg,nfftprc,ngfftprc,dtset%paral_kgb,0)
       call prcrskerker1(dtset,mpi_enreg,nfftprc,dtset%nspden,ngfftprc,dielar,etotal, &
&       gprimd,work1,work3,work)
       ABI_DEALLOCATE(work)
     else
       call prcrskerker2(dtset,nfftprc,dtset%nspden,ngfftprc,dielar,gprimd,rprimd, &
&       work1,work3,dtset%natom,xred,mpi_enreg,ucvol)
     end if
     do ispden=1,dtset%nspden
       call fourdp(1,vrespc(:,ispden),work3(:,ispden),-1,mpi_enreg,nfftprc,ngfftprc,dtset%paral_kgb,0)
     end do
     ABI_DEALLOCATE(work1)
     ABI_DEALLOCATE(work3)
   end if

 else

   if(dtset%iprcel==0 .or. (dtset%iprcel<40.and.istep<dielstrt) )then
     cplex=optreal
     qphon(:)=zero
!    Simple scalar multiplication, or model dielectric function
     call moddiel(cplex,dielar,mpi_enreg,nfftprc,ngfftprc,dtset%nspden,optreal,optres,dtset%paral_kgb,qphon,rprimd,vresid,vrespc)

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
&         nfft,ngfft,nkxc,npwdiel,dtset%nspden,dtset%occopt,option,dtset%paral_kgb,dtset%prtvol,susmat)
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
       call fourdp(1,work1,work2,-1,mpi_enreg,nfftprc,ngfftprc,dtset%paral_kgb,0)
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
       call fourdp(1,work1,work2,1,mpi_enreg,nfftprc,ngfftprc,dtset%paral_kgb,0)
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
   if (pawrhoij(1)%cplex==1) then
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
   call fourdp(1,rhog_wk,work,-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
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

   call hartre(1,gsqcut,psps%usepaw,mpi_enreg,nfft,ngfft,dtset%paral_kgb,rhog_wk,rprimd,vhartr_wk)

!  Prepare the call to rhotoxc
   call xcdata_init(xcdata,dtset=dtset)
   nk3xc=1
   ABI_ALLOCATE(work,(0))
   call rhotoxc(enxc,kxc,mpi_enreg,nfft,ngfft,&
&   work,0,work,0,nkxc,nk3xc,non_magnetic_xc,n3xccc,option,dtset%paral_kgb,rhor_wk,rprimd,strsxc,1,&
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
