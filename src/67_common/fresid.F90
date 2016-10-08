!{\src2tex{textfont=tt}}
!!****f* ABINIT/fresid
!!
!! NAME
!! fresid
!!
!! FUNCTION
!! If option=1, compute the forces due to the residual of the potential
!! If option=2, generate approximate new density from old one,
!!              old atomic positions, and new atomic positions
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (XG, MM, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dtset <type(dataset_type)>=all input variables in this dataset
!!  | icoulomb=0 periodic treatment of Hartree potential, 1 use of Poisson solver
!!  | natom=number of atoms in cell.
!!  | nspden=number of spin-density components
!!  | typat(natom)=integer type for each atom in cell
!!  | usepaw= 0 for non paw calculation; =1 for paw calculation
!!  | xclevel= level of the XC functional
!! mpi_enreg=information about MPI parallelization
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! ntypat=number of types of atoms in cell.
!! option=see below
!! pawtab(ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!! rhor(nfft,nspden)=electron density in electrons/bohr**3 (slices of it if FTT parallelism).
!! rprimd(3,3)=dimensional primitive translation vectors (bohr)
!! ucvol=unit cell volume (bohr**3).
!! xred_new(3,natom)=new reduced coordinates for atoms in unit cell
!! xred_old(3,natom)=old reduced coordinates for atoms in unit cell
!! znucl(ntypat)=real(dp), atomic number of atom type
!!
!! OUTPUT
!! gresid(3,natom)=forces due to the residual of the potential
!!
!! SIDE EFFECTS
!! work(nfft,nspden)=functions on the fft grid (slices of it if FTT parallelism):
!!  if option==1, the POTENTIAL residual is input
!!  if option==2, the interpolated density is output
!!
!! NOTES
!! FFT parallelism:
!! At the beginning of this routine, the plane-waves are ditributed over all the processors.
!! In the main part, all the processors perform the same calculations over the whole FFT grid.
!! At the end, each processor gets its part of the whole FFT grid.
!! These modifications are not efficient when large FFT grids are used.
!! So they have to be considered as a first step before a comprehensive parallelization of this routine.
!!
!! PARENTS
!!      forces,prcref,prcref_PMA,scfcv
!!
!! CHILDREN
!!      atomdata_from_znucl,mean_fftr,pre_gather,pre_scatter
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine fresid(dtset,gresid,mpi_enreg,nfft,ngfft,ntypat,option,&
&                 pawtab,rhor,rprimd,ucvol,work,xred_new,xred_old,znucl)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_cgtools
 use m_atomdata

 use m_pawtab,  only : pawtab_type
 use m_mpinfo,  only : pre_gather, pre_scatter

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fresid'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,ntypat,option
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(nfft,dtset%nspden),rprimd(3,3)
 real(dp),intent(in) :: xred_new(3,dtset%natom),xred_old(3,dtset%natom)
 real(dp),intent(in) :: znucl(ntypat)
 real(dp),intent(inout) :: work(nfft,dtset%nspden)
 real(dp),intent(out) :: gresid(3,dtset%natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*dtset%usepaw)

!Local variables-------------------------------
!real(dp), parameter :: app_remain=0.001_dp
!scalars
 integer,parameter :: natnum=110
 integer :: atmove,i1,i1_new,i1m,i1p,i2,i2_new,i2m,i2p,i3,i3_new,i3m,i3p
 integer :: iatom,ifft,ifft_new,iloop,ind2m,ind2m3m,ind2p,ind2p3p,ind3m,ind3p
 integer :: index,index_new,ishift,ishift1,ishift2,ishift3,ispden,ixp,mshift,mu
 integer :: n1,n2,n3,n4,nfft_tmp,nfftot,nu,quit
 real(dp),parameter :: app_remain=0.01_dp
 real(dp) :: diff_rem1,diff_rem2,diff_rem3,difmag,difmag2
 real(dp) :: difmag2_fact,difmag2_part,drho1,drho11,drho12,drho13,drho14
 real(dp) :: drho1dn,drho1mx,drho1my,drho1mz,drho1tot,drho1up,drho2,drho21
 real(dp) :: drho22,drho23,drho24,drho2dn,drho2mx,drho2my,drho2mz,drho2tot
 real(dp) :: drho2up,drho3,drho31,drho32,drho33,drho34,drho3dn,drho3mx,drho3my
 real(dp) :: drho3mz,drho3tot,drho3up,drhox00,drhox01,drhox10,drhox11,drhoxy0
 real(dp) :: drhoxy1,drhoxyz,fact,range,range2,rcov,rcov2,rcovm1,rdiff1
 real(dp) :: rdiff2,rdiff3,vresid1,vresid2,vresid3,vresid4,xx
 type(atomdata_t) :: atom
!arrays
 integer :: diff_igrid(3),igrid(3),irange(3)
 integer,allocatable :: ii(:,:)
 real(dp) :: diff_grid(3),diff_rem(3),diff_tau(3),diff_xred(3),lencp(3)
 real(dp) :: rho_tot(4),rhosum(4),rmet(3,3),scale(3),tau(3)
 real(dp),allocatable :: approp(:),atmrho(:,:),rhor_tot(:,:),rrdiff(:,:)
 real(dp),allocatable :: work_tot(:,:)
 logical,allocatable :: my_sphere(:)

! *************************************************************************

!Compute lengths of cross products for pairs of primitive
!translation vectors (used in setting index search range below)
 lencp(1)=cross_fr(rprimd(1,2),rprimd(2,2),rprimd(3,2),&
& rprimd(1,3),rprimd(2,3),rprimd(3,3))
 lencp(2)=cross_fr(rprimd(1,3),rprimd(2,3),rprimd(3,3),&
& rprimd(1,1),rprimd(2,1),rprimd(3,1))
 lencp(3)=cross_fr(rprimd(1,1),rprimd(2,1),rprimd(3,1),&
& rprimd(1,2),rprimd(2,2),rprimd(3,2))

!Compute factor R1.(R2xR3)/|R2xR3| etc for 1, 2, 3
!(recall ucvol=R1.(R2xR3))
 scale(:)=ucvol/lencp(:)

!initialize diff_igrid, otherwise valgrind complains 
 diff_igrid=0

!Compute metric tensor in real space rmet
 do nu=1,3
   rmet(:,nu)=rprimd(1,:)*rprimd(1,nu)+rprimd(2,:)*rprimd(2,nu)+rprimd(3,:)*rprimd(3,nu)
 end do

!FFT parallelization: Starting from now, calculations are performed on the whole FFT grid
!and no more on slices. The nfft variable becomes nfft_tmp until the end
 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3)
 n4=n3/mpi_enreg%nproc_fft
 nfftot=PRODUCT(ngfft(1:3));nfft_tmp=nfftot
 if(mpi_enreg%paral_kgb==1) then
   ABI_ALLOCATE(rhor_tot,(nfftot,dtset%nspden))
   ABI_ALLOCATE(work_tot,(nfftot,dtset%nspden))
   do ispden=1,dtset%nspden
     call pre_gather(rhor(:,ispden),rhor_tot(:,ispden),n1,n2,n3,n4,mpi_enreg)
     call pre_gather(work(:,ispden),work_tot(:,ispden),n1,n2,n3,n4,mpi_enreg)
   end do
 end if

 gresid(1:3,1:dtset%natom)=0.0_dp
 quit=0

!Initialize appropriation function
 ABI_ALLOCATE(approp,(nfft_tmp))
 ABI_ALLOCATE(atmrho,(nfft_tmp,dtset%nspden))
 ABI_ALLOCATE(my_sphere,(nfft_tmp))

 approp(:)=app_remain
!First loop over atoms in unit cell : build appropriation function
!Second loop : compute forces
 do iloop=1,2

!  Take into account the remaining density
   if(option==2 .and. iloop==2)then
     if(mpi_enreg%paral_kgb==1) then
!      FFT parallelization: All the processors perform the same calculation.
!      We divided by nproc_fft in order to "remove" the xmpi_sum made in mean_fftr
       do ispden=1,dtset%nspden
         do ifft=1,nfft_tmp
           work_tot(ifft,ispden)=rhor_tot(ifft,ispden)*approp(ifft)*app_remain
         end do
       end do
       call mean_fftr(work_tot,rhosum,nfft_tmp,nfftot,dtset%nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
       rhosum(1:dtset%nspden)=rhosum(1:dtset%nspden)/mpi_enreg%nproc_fft
     else
       do ispden=1,dtset%nspden
         do ifft=1,nfft_tmp
           work(ifft,ispden)=rhor(ifft,ispden)*approp(ifft)*app_remain
         end do
       end do
       call mean_fftr(work,rhosum,nfft_tmp,nfftot,dtset%nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
     end if

!    This will be used to restore proper normalization of density
     rho_tot(1:dtset%nspden)=rhosum(1:dtset%nspden)*nfftot
   end if

   do iatom=1,dtset%natom

!    Get the covalent radius
     call atomdata_from_znucl(atom,znucl(dtset%typat(iatom))) 
     rcov = atom%rcov
!    PAW choose PAW radius instead...
     if (dtset%usepaw==1) rcov=max(rcov,pawtab(dtset%typat(iatom))%rpaw)

!    Set search range
     rcov2=rcov**2
     range=2._dp*rcov
     range2=range**2
     rcovm1=1.0_dp/rcov

!    Use range to compute an index range along R(1:3)
!    (add 1 to make sure it covers full range)
     irange(1)=1+nint((range/scale(1))*dble(n1))
     irange(2)=1+nint((range/scale(2))*dble(n2))
     irange(3)=1+nint((range/scale(3))*dble(n3))

!    Allocate ii and rrdiff
     mshift=2*maxval(irange(1:3))+1
     ABI_ALLOCATE(ii,(mshift,3))
     ABI_ALLOCATE(rrdiff,(mshift,3))

!    Consider each component in turn
     do mu=1,3

!      Convert reduced coord of given atom to [0,1)
       tau(mu)=mod(xred_old(mu,iatom)+1._dp-aint(xred_old(mu,iatom)),1._dp)

!      Use tau to find nearest grid point along R(mu)
!      (igrid=0 is the origin; shift by 1 to agree with usual index)
       igrid(mu)=nint(tau(mu)*dble(ngfft(mu)))

!      Set up a counter that explore the relevant range
!      of points around the atom
       ishift=0
       do ixp=igrid(mu)-irange(mu),igrid(mu)+irange(mu)
         ishift=ishift+1
         ii(ishift,mu)=1+mod(ngfft(mu)+mod(ixp,ngfft(mu)),ngfft(mu))
         rrdiff(ishift,mu)=dble(ixp)/dble(ngfft(mu))-tau(mu)
       end do

!      If option 2, set up quantities related with the change of atomic coordinates
       if(option==2 .and. iloop==2)then
         diff_xred(mu)=xred_new(mu,iatom)-xred_old(mu,iatom)
!        Convert to [0,1)
         diff_tau(mu)=mod(diff_xred(mu)+1._dp-aint(diff_xred(mu)),1._dp)
!        Convert to [0,ngfft)
         diff_grid(mu)=diff_tau(mu)*dble(ngfft(mu))
!        Integer part
         diff_igrid(mu)=int(diff_grid(mu))
!        Compute remainder
         diff_rem(mu)=diff_grid(mu)-diff_igrid(mu)

!        DEBUG
!        write(std_out,*)' mu,diff',mu,diff_igrid(mu),diff_rem(mu)
!        ENDDEBUG

       end if

!      End loop on mu
     end do

!    May be the atom is fixed
     atmove=1
     if(option==2 .and. iloop==2)then
       if(diff_xred(1)**2+diff_xred(2)**2+diff_xred(3)**2 < 1.0d-24)then
         atmove=0
       else
         diff_rem1=diff_rem(1)
         diff_rem2=diff_rem(2)
         diff_rem3=diff_rem(3)
       end if
     end if

!    If second loop, initialize atomic density, and the variable
!    that says whether a fft point belongs to the sphere of the atom
     if(iloop==2) then
       atmrho(:,:)=0.0_dp
       my_sphere(:)=.false.
     end if

!    Conduct triple loop over restricted range of grid points for iatom

     do ishift3=1,1+2*irange(3)
!      map back to [1,ngfft(3)] for usual fortran index in unit cell
       i3=ii(ishift3,3)
       i3m=i3-1 ; if(i3==1)i3m=n3
       i3p=i3+1 ; if(i3==n3)i3p=1

!      find vector from atom location to grid point (reduced)
       rdiff3=rrdiff(ishift3,3)

       do ishift2=1,1+2*irange(2)
         i2=ii(ishift2,2)
         i2m=i2-1 ; if(i2==1)i2m=n2
         i2p=i2+1 ; if(i2==n2)i2p=1
         index=n1*(i2-1+n2*(i3-1))
         ind3m=n1*(i2-1+n2*(i3m-1))
         ind3p=n1*(i2-1+n2*(i3p-1))
         ind2m=n1*(i2m-1+n2*(i3-1))
         ind2p=n1*(i2p-1+n2*(i3-1))
         ind2p3p=n1*(i2p-1+n2*(i3p-1))

         rdiff2=rrdiff(ishift2,2)
!        Prepare the computation of difmag2
         difmag2_part=rmet(3,3)*rdiff3**2+rmet(2,2)*rdiff2**2&
&         +2.0_dp*rmet(3,2)*rdiff3*rdiff2
         difmag2_fact=2.0_dp*(rmet(3,1)*rdiff3+rmet(2,1)*rdiff2)

         do ishift1=1,1+2*irange(1)
           rdiff1=rrdiff(ishift1,1)

!          Compute (rgrid-tau-Rprim)**2
           difmag2= difmag2_part+rdiff1*(difmag2_fact+rmet(1,1)*rdiff1)

!          Only accept contribution inside defined range
!          This condition means that x, calculated below, cannot exceed 2.0_dp
           if (difmag2<range2) then

!            Will compute contribution to appropriation function based on
!            rcov2, range2 and difmag2
             i1=ii(ishift1,1)
             ifft=i1+index

             if(iloop==1)then

!              Build appropriation function
               if (difmag2<rcov2)then
                 approp(ifft)=approp(ifft)+1.0_dp
               else
                 difmag=sqrt(difmag2)
                 xx=difmag*rcovm1
!                The following function is 1. at xx=1, 0. at xx=2, with vanishing
!                derivatives at these points.
                 approp(ifft)=approp(ifft)+((2.0_dp*xx-9.0_dp)*xx+12.0_dp)*xx-4.0_dp
               end if

             else

               if (difmag2<rcov2) then
                 fact=one
               else
                 difmag=sqrt(difmag2)
                 xx=difmag*rcovm1
                 fact=((2.0_dp*xx-9.0_dp)*xx+12.0_dp)*xx-4.0_dp
               end if

!              Build atomic density
               if(mpi_enreg%paral_kgb==1) then
                 atmrho(ifft,1:dtset%nspden)=atmrho(ifft,1:dtset%nspden) &
&                 +rhor_tot(ifft,1:dtset%nspden)*fact*approp(ifft)
               else
                 atmrho(ifft,1:dtset%nspden)=atmrho(ifft,1:dtset%nspden) &
&                 +rhor(ifft,1:dtset%nspden)*fact*approp(ifft)
               end if

!              Compute the sphere of the atom : it is different for
!              option 1 and for option 2
               i1p=i1+1 ; if(i1==n1)i1p=1
               if(option==1)then
                 i1m=i1-1 ; if(i1==1)i1m=n1
                 my_sphere(ifft)=.true.
                 my_sphere(i1p+index)=.true. ; my_sphere(i1m+index)=.true.
                 my_sphere(i1+ind2p)=.true. ; my_sphere(i1+ind2m)=.true.
                 my_sphere(i1+ind3p)=.true. ; my_sphere(i1+ind3m)=.true.
               else
                 my_sphere(ifft)=.true. ; my_sphere(i1p+index)=.true.
                 my_sphere(i1+ind2p)=.true. ; my_sphere(i1p+ind2p)=.true.
                 my_sphere(i1+ind3p)=.true. ; my_sphere(i1p+ind3p)=.true.
                 my_sphere(i1+ind2p3p)=.true. ; my_sphere(i1p+ind2p3p)=.true.
               end if

             end if

!            End of condition on the range
           end if

!          End loop on ishift1
         end do

!        End loop on ishift2
       end do

!      End loop on ishift3
     end do
!    At the end of the second loop for each atom, compute the force
!    from the atomic densities, or translate density.
!    In the first case, use a two-point finite-difference approximation,
!    since this calculation serves only to decrease the error,
!    and should not be very accurate, but fast.
!    In the second case, using a crude trilinear interpolation scheme
!    for the same reason.
!    
!    The section is skipped if option==2 and the atom is fixed
     if(iloop==2 .and. (option==1 .or. atmove==1) )then

       do i3=1,n3
         i3m=i3-1 ; if(i3==1)i3m=n3
         i3p=i3+1 ; if(i3==n3)i3p=1
!        note: diff_igrid is only set  if(option==2 .and. iloop==2)
         i3_new=i3+diff_igrid(3) ; if(i3_new > n3)i3_new=i3_new-n3
         do i2=1,n2
           i2m=i2-1 ; if(i2==1)i2m=n2
           i2p=i2+1 ; if(i2==n2)i2p=1
           i2_new=i2+diff_igrid(2) ; if(i2_new > n2)i2_new=i2_new-n2
           index=n1*(i2-1+n2*(i3-1))
           index_new=n1*(i2_new-1+n2*(i3_new-1))
           ind3m=n1*(i2-1+n2*(i3m-1))
           ind3p=n1*(i2-1+n2*(i3p-1))
           ind2m=n1*(i2m-1+n2*(i3-1))
           ind2p=n1*(i2p-1+n2*(i3-1))
           ind2m3m=n1*(i2m-1+n2*(i3m-1))
           do i1=1,n1
             ifft=i1+index
             if(my_sphere(ifft))then

               i1m=i1-1 ; if(i1==1)i1m=n1

               if(option==1)then
!                Treat option 1 : computation of residual forces
                 i1p=i1+1 ; if(i1==n1)i1p=1
!                Distinguish spin-unpolarized and spin-polarized
                 if(dtset%nspden==1)then ! Non magnetic
!                  Note that the factor needed to obtain a true finite difference
!                  estimation of the derivative will be applied afterwards, for speed
                   drho1=atmrho(i1p+index,1)-atmrho(i1m+index,1)
                   drho2=atmrho(i1+ind2p,1) -atmrho(i1+ind2m,1)
                   drho3=atmrho(i1+ind3p,1) -atmrho(i1+ind3m,1)
                   if(mpi_enreg%paral_kgb==1) then
                     vresid1=work_tot(ifft,1)
                   else
                     vresid1=work(ifft,1)
                   end if
                   gresid(1,iatom)=gresid(1,iatom)+drho1*vresid1
                   gresid(2,iatom)=gresid(2,iatom)+drho2*vresid1
                   gresid(3,iatom)=gresid(3,iatom)+drho3*vresid1
                 else if(dtset%nspden==2) then ! Collinear magnetism
                   drho1tot=atmrho(i1p+index,1)-atmrho(i1m+index,1)
                   drho2tot=atmrho(i1+ind2p,1) -atmrho(i1+ind2m,1)
                   drho3tot=atmrho(i1+ind3p,1) -atmrho(i1+ind3m,1)
                   drho1up=atmrho(i1p+index,2)-atmrho(i1m+index,2)
                   drho2up=atmrho(i1+ind2p,2) -atmrho(i1+ind2m,2)
                   drho3up=atmrho(i1+ind3p,2) -atmrho(i1+ind3m,2)
                   drho1dn=drho1tot-drho1up
                   drho2dn=drho2tot-drho2up
                   drho3dn=drho3tot-drho3up
                   if(mpi_enreg%paral_kgb==1) then
                     vresid1=work_tot(ifft,1)
                     vresid2=work_tot(ifft,2)
                   else
                     vresid1=work(ifft,1)
                     vresid2=work(ifft,2)
                   end if
                   gresid(1,iatom)=gresid(1,iatom)+drho1up*vresid1+drho1dn*vresid2
                   gresid(2,iatom)=gresid(2,iatom)+drho2up*vresid1+drho2dn*vresid2
                   gresid(3,iatom)=gresid(3,iatom)+drho3up*vresid1+drho3dn*vresid2
                 else ! Non-collinear magnetism
                   drho1tot=atmrho(i1p+index,1)-atmrho(i1m+index,1)
                   drho1mx =atmrho(i1p+index,2)-atmrho(i1m+index,2)
                   drho1my =atmrho(i1p+index,3)-atmrho(i1m+index,3)
                   drho1mz =atmrho(i1p+index,4)-atmrho(i1m+index,4)
                   drho2tot=atmrho(i1+ind2p,1) -atmrho(i1+ind2m,1)
                   drho2mx =atmrho(i1+ind2p,2) -atmrho(i1+ind2m,2)
                   drho2my =atmrho(i1+ind2p,3) -atmrho(i1+ind2m,3)
                   drho2mz =atmrho(i1+ind2p,4) -atmrho(i1+ind2m,4)
                   drho3tot=atmrho(i1+ind3p,1) -atmrho(i1+ind3m,1)
                   drho3mx =atmrho(i1+ind3p,2) -atmrho(i1+ind3m,2)
                   drho3my =atmrho(i1+ind3p,3) -atmrho(i1+ind3m,3)
                   drho3mz =atmrho(i1+ind3p,4) -atmrho(i1+ind3m,4)
                   drho11=half*(drho1tot+drho1mz)
                   drho12=half*(drho1tot-drho1mz)
                   drho13= half*drho1mx
                   drho14=-half*drho1my
                   drho21=half*(drho2tot+drho2mz)
                   drho22=half*(drho2tot-drho2mz)
                   drho23= half*drho2mx
                   drho24=-half*drho2my
                   drho31=half*(drho3tot+drho3mz)
                   drho32=half*(drho3tot-drho3mz)
                   drho33= half*drho3mx
                   drho34=-half*drho3my
                   if(mpi_enreg%paral_kgb==1) then
                     vresid1=work_tot(ifft,1)
                     vresid2=work_tot(ifft,2)
                     vresid3=work_tot(ifft,3)
                     vresid4=work_tot(ifft,4)
                   else
                     vresid1=work(ifft,1)
                     vresid2=work(ifft,2)
                     vresid3=work(ifft,3)
                     vresid4=work(ifft,4)
                   end if
                   gresid(1,iatom)=gresid(1,iatom)+drho11*vresid1+drho12*vresid2+two*(drho13*vresid3+drho14*vresid4)
                   gresid(2,iatom)=gresid(2,iatom)+drho21*vresid1+drho22*vresid2+two*(drho23*vresid3+drho24*vresid4)
                   gresid(3,iatom)=gresid(3,iatom)+drho31*vresid1+drho32*vresid2+two*(drho33*vresid3+drho34*vresid4)
                 end if
!                Treat the case option==2 now : trilinear interpolation of the density
               else
                 i1_new=i1+diff_igrid(1) ; if(i1_new > n1)i1_new=i1_new-n1
                 ifft_new=i1_new+index_new
                 do ispden=1,dtset%nspden
                   drhox00=(atmrho(i1m+index,ispden)-atmrho(i1+index,ispden))*diff_rem1 &
&                   +atmrho(i1+index,ispden)
                   drhox10=(atmrho(i1m+ind2m,ispden)-atmrho(i1+ind2m,ispden))*diff_rem1 &
&                   +atmrho(i1+ind2m,ispden)
                   drhox01=(atmrho(i1m+ind3m,ispden)-atmrho(i1+ind3m,ispden))*diff_rem1 &
&                   +atmrho(i1+ind3m,ispden)
                   drhox11=(atmrho(i1m+ind2m3m,ispden)-atmrho(i1+ind2m3m,ispden))*diff_rem1 &
&                   +atmrho(i1+ind2m3m,ispden)
                   drhoxy0=(drhox10-drhox00)*diff_rem2+drhox00
                   drhoxy1=(drhox11-drhox01)*diff_rem2+drhox01
                   drhoxyz=(drhoxy1-drhoxy0)*diff_rem3+drhoxy0
                   if(mpi_enreg%paral_kgb==1) then
                     work_tot(ifft_new,ispden)=work_tot(ifft_new,ispden)+drhoxyz
                   else
                     work(ifft_new,ispden)=work(ifft_new,ispden)+drhoxyz
                   end if
                   rho_tot(ispden)=rho_tot(ispden)+drhoxyz
                 end do
               end if

!              End condition of belonging to the sphere of influence of the atom
             end if
           end do
         end do
       end do
!      The finite-difference factor applied here also take
!      into account diverse factors
       fact=-ucvol/dble(nfftot)
       gresid(1,iatom)=gresid(1,iatom)*dble(n1)*.5_dp*fact
       gresid(2,iatom)=gresid(2,iatom)*dble(n2)*.5_dp*fact
       gresid(3,iatom)=gresid(3,iatom)*dble(n3)*.5_dp*fact
     end if

!    Update work if the atom is fixed.
     if(iloop==2 .and. option==2 .and. atmove==0)then
       if(mpi_enreg%paral_kgb==1) then
!        FFT parallelization: All the processors perform the same calculation.
!        We divided by nproc_fft in order to "remove" the xmpi_sum made in mean_fftr
         do ispden=1,dtset%nspden
           do ifft=1,nfft_tmp
             work_tot(ifft,ispden)=work_tot(ifft,ispden)+atmrho(ifft,ispden)
           end do
         end do
         call mean_fftr(atmrho,rhosum,nfft_tmp,nfftot,dtset%nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
         rhosum(1:dtset%nspden)=rhosum(1:dtset%nspden)/mpi_enreg%nproc_fft
       else
         do ispden=1,dtset%nspden
           do ifft=1,nfft_tmp
             work(ifft,ispden)=work(ifft,ispden)+atmrho(ifft,ispden)
           end do
         end do
         call mean_fftr(atmrho,rhosum,nfft_tmp,nfftot,dtset%nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
       end if

       rho_tot(1:dtset%nspden)=rho_tot(1:dtset%nspden)+rhosum(1:dtset%nspden)*nfftot
     end if

     ABI_DEALLOCATE(ii)
     ABI_DEALLOCATE(rrdiff)

!    End loop on atoms
   end do

!  DEBUG
!  if(option==2)then
!  if(iloop==1)then
!  write(std_out,*)' fresid : rhor, approp'
!  do ifft=1,n1
!  write(std_out,*)ifft,rhor(ifft,1),approp(ifft)
!  end do
!  end if
!  if(iloop==2)then
!  write(std_out,*)' fresid : rhor, approp, work(:,:)'
!  do ifft=1,n1
!  write(std_out,'(i4,3es18.8)' )ifft,rhor(ifft,1),approp(ifft),work(ifft,1)
!  end do
!  do ifft=1,nfft_tmp
!  if(work(ifft,1)<0.0_dp)then
!  write(std_out,*)' f_fft negative value',work(ifft,1),' for ifft=',ifft
!  end if
!  if(rhor(ifft,1)<0.0_dp)then
!  write(std_out,*)' rhor  negative value',rhor(ifft,1),' for ifft=',ifft
!  end if
!  end do
!  end if
!  end if
!  ENDDEBUG

   if(quit==1)exit

!  At the end of the first loop, where the appropriation function is generated,
!  invert it, to save cpu time later.
   if(iloop==1)approp(:)=1.0_dp/approp(:)

!  End first or second pass through atoms
 end do

!Restore proper normalisation of density
!(Non-collinear magnetism: n, mx,my,mz integral conservation)
 if(option==2)then
   if(mpi_enreg%paral_kgb==1) then
!    FFT parallelization: All the processors perform the same calculation.
!    We divided by nproc_fft in order to "remove" the xmpi_sum made in mean_fftr
!    Trangel: mpicomm now is optional in mean_fftr, no need to divide over nproc_fft
     call mean_fftr(rhor_tot,rhosum,nfft_tmp,nfftot,dtset%nspden)
!    rhosum(1:dtset%nspden)=rhosum(1:dtset%nspden)/mpi_enreg%nproc_fft
   else
     call mean_fftr(rhor,rhosum,nfft_tmp,nfftot,dtset%nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
   end if
!  "!OCL NOPREEX" to avoid zero division after optimization (-Of) by MM
!  (Even if nspden=1, "1.0/rho_tot" will appear on vpp fujitsu
!  OCL NOPREEX
   if(mpi_enreg%paral_kgb==1) then
     do ispden=1,dtset%nspden
       fact=rhosum(ispden)*dble(nfftot)/rho_tot(ispden)
       work_tot(:,ispden)=fact*work_tot(:,ispden)
       call pre_scatter(work(:,ispden),work_tot(:,ispden),n1,n2,n3,n4,mpi_enreg)
     end do
   else
     do ispden=1,dtset%nspden
       fact=rhosum(ispden)*dble(nfftot)/rho_tot(ispden)
       work(:,ispden)=fact*work(:,ispden)
     end do
   end if
!  DEBUG
!  Here, zero all the hard work, for checking purposes !
!  work(:,:)=rhor(:,:)
!  ENDDEBUG
 end if

 ABI_DEALLOCATE(approp)
 ABI_DEALLOCATE(atmrho)
 ABI_DEALLOCATE(my_sphere)
 if(mpi_enreg%paral_kgb==1) then
   ABI_DEALLOCATE(rhor_tot)
   ABI_DEALLOCATE(work_tot)
 end if

!DEBUG
!write(std_out,*)' fresid : exit '
!do iatom=1,dtset%natom
!write(std_out,*)iatom,gresid(1:3,iatom)
!enddo
!ENDDEBUG

 contains

   real(dp) pure function cross_fr(xx,yy,zz,aa,bb,cc)
!Define magnitude of cross product of two vectors


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cross_fr'
!End of the abilint section

   real(dp),intent(in) :: xx,yy,zz,aa,bb,cc
   cross_fr=sqrt((yy*cc-zz*bb)**2+(zz*aa-xx*cc)**2+(xx*bb-yy*aa)**2)
 end function cross_fr

end subroutine fresid
!!***
