!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_dyxc1
!! NAME
!! dfpt_dyxc1
!!
!! FUNCTION
!! Compute 2nd-order non-linear xc core-correction (part1) to the dynamical matrix.
!! In case of derivative with respect to k or electric field perturbation, 
!! the 1st-order local potential vanishes.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!  ixc= choice of exchange-correlation scheme
!!  kxc(nfft,nkxc)=first-order derivative of the xc potential
!!    if (nkxc=1) LDA kxc(:,1)= d2Exc/drho2
!!    if (nkxc=2) LDA kxc(:,1)=d2Exc/drho_up drho_up
!!                    kxc(:,2)=d2Exc/drho_up drho_dn
!!                    kxc(:,3)=d2Exc/drho_dn drho_dn
!!    if (nkxc=7) GGA kxc(:,1)= d2Exc/drho2
!!                    kxc(:,2)= 1/|grad(rho)| dExc/d|grad(rho)|
!!                    kxc(:,3)= 1/|grad(rho)| d2Exc/d|grad(rho)| drho
!!                    kxc(:,4)= 1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dExc/d|grad(rho)| )
!!                    kxc(:,5)= gradx(rho)
!!                    kxc(:,6)= grady(rho)
!!                    kxc(:,7)= gradz(rho)
!!    if (nkxc=19) spin-polarized GGA case (same as nkxc=7 with up and down components)
!!  mgfft=maximum size of 1D FFTs
!!  mpert=maximum number of ipert
!!  mpi_enreg=information about MPI parallelization
!!  mqgrid=number of grid pts in q array for f(q) spline.
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(3)=fft grid dimensions.
!!  nkxc=second dimension of the kxc array
!!   (=1 for non-spin-polarized case, =3 for spin-polarized case)
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms in cell.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  qgrid(mqgrid)=q grid for spline from 0 to qmax.
!!  qphon(3)=wavevector of the phonon
!!  rfdir(3)=array that define the directions of perturbations
!!  rfpert(mpert)=array defining the type of perturbations that have to be computed
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  timrev=1 if time-reversal preserves the q wavevector; 0 otherwise.
!!  typat(natom)=integer type for each atom in cell
!!  ucvol=unit cell volume (bohr**3).
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc,6,ntypat)=1D core charge function and five derivatives,
!!   for each type of atom, from psp
!!  xred(3,natom)=fractional coordinates for atoms in unit cell
!!
!! OUTPUT
!!  blkflgfrx1(3,natom,3,natom)=flag to indicate whether an element has been computed or not
!!  dyfrx1(2,3,natom,3,natom)=2nd-order non-linear xc
!!    core-correction (part1) part of the dynamical matrix
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      dfpt_atm2fft,dfpt_mkcore,dfpt_mkvxc,dfpt_mkvxc_noncoll,dotprod_vn,timab
!!      xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dfpt_dyxc1(atindx,blkflgfrx1,dyfrx1,gmet,gsqcut,ixc,kxc,mgfft,mpert,mpi_enreg,mqgrid,&
&          natom,nfft,ngfft,nkxc,nspden,ntypat,n1xccc,paral_kgb,psps,pawtab,&
&          ph1d,qgrid,qphon,rfdir,rfpert,rprimd,timrev,typat,ucvol,usepaw,xcccrc,xccc1d,xred,rhor,vxc)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi

 use defs_datatypes,  only : pseudopotential_type
 use m_cgtools,       only : dotprod_vn
 use m_pawtab,        only : pawtab_type
 use m_xmpi,          only : xmpi_sum

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_dyxc1'
 use interfaces_18_timing
 use interfaces_56_xc
 use interfaces_64_psp
 use interfaces_72_response, except_this_one => dfpt_dyxc1
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,mgfft,mpert,mqgrid,n1xccc,natom,nfft,nkxc,nspden,ntypat
 integer,intent(in) :: paral_kgb,timrev,usepaw
 real(dp),intent(in) :: gsqcut,ucvol
 type(pseudopotential_type),intent(in) :: psps
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),ngfft(18),rfdir(3),rfpert(mpert),typat(natom)
 real(dp),intent(in) :: gmet(3,3),kxc(nfft,nkxc)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qgrid(mqgrid),qphon(3)
 real(dp),intent(in) :: rprimd(3,3),xccc1d(n1xccc,6,ntypat),xcccrc(ntypat)
 real(dp),intent(in) :: xred(3,natom)
 integer,intent(out) :: blkflgfrx1(3,natom,3,natom)
 real(dp),intent(out) :: dyfrx1(2,3,natom,3,natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
!optional
 real(dp),optional,intent(in) :: rhor(nfft,nspden)
 real(dp),optional,intent(in) :: vxc(nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: cplex,iat1,iatom1,iatom2,idir1,idir2,ierr,ifft,my_natom,comm_atom
 integer :: n1,n2,n3,n3xccc,nfftot,option,upperdir,optnc
 logical :: paral_atom
 real(dp) :: valuei,valuer
!arrays
 integer,pointer :: my_atmtab(:)
 real(dp) :: tsec(2),gprimd_dummy(3,3)
 real(dp) :: dum_nhat(0)
 real(dp),allocatable :: rhor1(:,:),vxc10(:,:),xcccwk1(:),xcccwk2(:)
! *********************************************************************

 call timab(182,1,tsec)

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 nfftot=n1*n2*n3

!Set up parallelism over atoms
 my_natom=mpi_enreg%my_natom
 my_atmtab=>mpi_enreg%my_atmtab
 comm_atom=mpi_enreg%comm_atom
 paral_atom=(my_natom/=natom)

!Zero out the output arrays :
 blkflgfrx1(:,:,:,:)=0
 dyfrx1(:,:,:,:,:)=zero

 cplex=2-timrev ; n3xccc=nfft
 ABI_ALLOCATE(vxc10,(cplex*nfft,nspden))


!Loop on the perturbation j1
 do iat1=1,my_natom
   iatom1=iat1; if(paral_atom)iatom1=my_atmtab(iat1)
   do idir1=1,3

!    Compute the derivative of the core charge with respect to j1
     ABI_ALLOCATE(xcccwk1,(cplex*n3xccc))

!    PAW or NC with nc_xccc_gspace: 1st-order core charge in reciprocal space
     if (usepaw==1 .or. psps%nc_xccc_gspace==1) then
       call dfpt_atm2fft(atindx,cplex,gmet,gprimd_dummy,gsqcut,idir1,iatom1,&
&       mgfft,mqgrid,natom,1,nfft,ngfft,ntypat,ph1d,qgrid,&
&       qphon,typat,ucvol,usepaw,xred,psps,pawtab,atmrhor1=xcccwk1,optn2_in=1)

!      Norm-conserving psp: 1st-order core charge in real space
     else
       call dfpt_mkcore(cplex,idir1,iatom1,natom,ntypat,n1,n1xccc,&
&       n2,n3,qphon,rprimd,typat,ucvol,xcccrc,xccc1d,xcccwk1,xred)
     end if

!    Compute the corresponding potential
     option=0
     ABI_ALLOCATE(rhor1,(cplex*nfft,nspden))
     rhor1=zero
!FR SPr EB Non-collinear magnetism
     if (nspden==4.and.present(rhor).and.present(vxc)) then
       optnc=1
       call dfpt_mkvxc_noncoll(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,dum_nhat,0,dum_nhat,0,dum_nhat,0,nkxc,&
&       nspden,n3xccc,optnc,option,paral_kgb,qphon,rhor,rhor1,rprimd,0,vxc,vxc10,xcccwk1)
     else
       call dfpt_mkvxc(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,dum_nhat,0,dum_nhat,0,nkxc,&
&       nspden,n3xccc,option,paral_kgb,qphon,rhor1,rprimd,0,vxc10,xcccwk1)
     end if
     ABI_DEALLOCATE(rhor1)
     ABI_DEALLOCATE(xcccwk1)

!    vxc10 will couple with xcccwk2, that behaves like
!    a total density (ispden=1). Only the spin-up + spin-down
!    average of vxc10 is needed.
     if (nspden/=1)then
       do ifft=1,cplex*nfft
         vxc10(ifft,1)=(vxc10(ifft,1)+vxc10(ifft,2))*half
       end do
     end if

!    Loop on the perturbation j2
     do iatom2=1,iatom1
       upperdir=3
       if(iatom1==iatom2)upperdir=idir1
       do idir2=1,upperdir
         if( (rfpert(iatom1)==1 .and. rfdir(idir1) == 1) .or. &
&         (rfpert(iatom2)==1 .and. rfdir(idir2) == 1)    )then

!          Compute the derivative of the core charge with respect to j2
           ABI_ALLOCATE(xcccwk2,(cplex*n3xccc))

!          PAW or NC with nc_xccc_gspace: 1st-order core charge in reciprocal space
           if (usepaw==1 .or. psps%nc_xccc_gspace==1) then
             call dfpt_atm2fft(atindx,cplex,gmet,gprimd_dummy,gsqcut,idir2,iatom2,&
&             mgfft,mqgrid,natom,1,nfft,ngfft,ntypat,ph1d,qgrid,&
&             qphon,typat,ucvol,usepaw,xred,psps,pawtab,atmrhor1=xcccwk2,optn2_in=1)

!            Norm-conserving psp: 1st-order core charge in real space
           else
             call dfpt_mkcore(cplex,idir2,iatom2,natom,ntypat,n1,n1xccc,&
&             n2,n3,qphon,rprimd,typat,ucvol,xcccrc,xccc1d,xcccwk2,xred)
           end if

!          Get the matrix element j1,j2

           call dotprod_vn(cplex,xcccwk2,valuer,valuei,nfft,nfftot,1,2,vxc10,ucvol)

           ABI_DEALLOCATE(xcccwk2)

           dyfrx1(1,idir1,iatom1,idir2,iatom2)= valuer
           dyfrx1(2,idir1,iatom1,idir2,iatom2)= valuei
           dyfrx1(1,idir2,iatom2,idir1,iatom1)= valuer
           dyfrx1(2,idir2,iatom2,idir1,iatom1)=-valuei
           blkflgfrx1(idir1,iatom1,idir2,iatom2)=1
           blkflgfrx1(idir2,iatom2,idir1,iatom1)=1
         end if
       end do
     end do
   end do
 end do

 if (paral_atom) then
   call timab(48,1,tsec)
   call xmpi_sum(dyfrx1,comm_atom,ierr)
   call xmpi_sum(blkflgfrx1,comm_atom,ierr)
   call timab(48,2,tsec)
 end if

 ABI_DEALLOCATE(vxc10)

 call timab(182,2,tsec)

end subroutine dfpt_dyxc1
!!***
