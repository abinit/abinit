!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_dyout
!! NAME
!! dfpt_dyout
!!
!!
!! FUNCTION
!! Output of all quantities related to the 2nd-order matrix :
!! Ewald part, local and non-local frozen wf part,
!! core contributions,
!! local and non-local variational part, 2nd-order
!! matrix itself, and, for the phonon part,
!! eigenfrequencies, in Hartree, meV and cm-1.
!! Also output unformatted 2nd-order matrix for later
!! use in the Brillouin-zone interpolation
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (XG, DRH, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  becfrnl(3,natom,3*pawbec)=NL frozen contribution to Born Effective Charges (PAW only)
!!  blkflg(3,mpert,3,mpert)= ( 1 if the element of the dynamical
!!  matrix has been calculated ; 0 otherwise )
!!  carflg(3,mpert,3,mpert)= ( 1 if the element of the cartesian
!!  2DTE matrix has been calculated correctly ; 0 otherwise )
!!  ddboun=unit number for the derivative database output
!!  ddkfil(3)=components are 1 if corresponding d/dk file exists, otherwise 0
!!  (in what follows, DYMX means dynamical matrix, and D2MX means 2nd-order matrix)
!!  dyew(2,3,natom,3,natom)=Ewald part of the DYMX
!!  dyfrlo(3,3,natom)=frozen wf local part of the DYMX
!!  dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)=frozen wf nonloc part of the DYMX
!!  dyfrx1(2,3,natom,3,natom)=frozen wf nonlin. xc core corr.(1)
!!    part of the DYMX
!!  dyfrx2(3,3,natom)=frozen wf nonlin. xc core corr.(2) part of the DYMX
!!  dyfr_cplex=1 if dyfrnl is real, 2 if it is complex
!!  dyfr_nondiag=1 if dyfrnl is non diagonal with respect to atoms; 0 otherwise
!!  dyvdw(2,3,natom,3,natom*usevdw)=vdw DFT-D part of the dynamical matrix
!!  d2cart(2,3,mpert,3,mpert)=D2MX in cartesian coordinates
!!  d2cart_bbb(2,3,3,mpert,mband,mband*prtbbb)=
!!    band by band decomposition of Born effective charges
!!    (calculated from phonon-type perturbation) in cartesian coordinates
!!  d2eig0(2,3,mpert,3,mpert)=0-order eigen. station. part of the D2MX
!!  d2k0(2,3,mpert,3,mpert)=0-order kinet. station. part of the D2MX
!!  d2lo(2,3,mpert,3,mpert)=nonstation. local part of the D2MX
!!  d2loc0(2,3,mpert,3,mpert)=0-order loc station. part of the D2MX
!!  d2matr(2,3,mpert,3,mpert)=D2MX in non-cartesian coordinates
!!  d2nl(2,3,mpert,3,mpert)=nonstation. nonloc part of the D2MX
!!  d2nl0(2,3,mpert,3,mpert)=0-order nonloc station. part of the D2MX
!!  d2nl1(2,3,mpert,3,mpert)=1-order nonloc station. part of the D2MX
!!  d2ovl(2,mpert,3,mpert*usepaw)=1st-order change of WF overlap contributions to the 2DTEs (PAW)
!!  d2vn(2,3,mpert,3,mpert)=potential*dens station. part of the D2MX and without masses included)
!!  eltcore(6,6)=core contribution to the elastic tensor
!!  elteew(6+3*natom,6)=Ewald contribution to the elastic tsenor
!!  eltfrhar(6,6)=hartree contribution to the elastic tensor
!!  eltfrkin(6,6)=kinetic contribution to the elastic tensor
!!  eltfrloc(6+3*natom,6)=local psp contribution to the elastic tensor
!!  eltfrnl(6+3*natom,6)=non-local psp contribution to the elastic tensor
!!  eltfrxc(6+3*natom,6)=exchange-correlation contribution to the elastic tensor
!!  eltvdw(6+3*natom,6*usevdw)=vdw DFT-D part of the elastic tensor
!!  has_full_piezo=the full calculation of the piezoelectric tensor from electric field perturbation
!!                 is only available if nsym=1 (strain perturbation is not symmetrized)
!!  has_allddk= True if all ddk file are present on disk, so the effective charge or piezzo
!!              electric tensor are correctly computed (PAW ONLY)
!!  iout=unit number for long write-up
!!  mband=maximum number of bands
!!  mpert =maximum number of ipert
!!  natom=number of atoms
!!  ntypat=number of atom types
!!  outd2=option for the output of the 2nd-order matrix :
!!   if outd2=1, non-stationary part
!!   if outd2=2, stationary part.
!!  pawbec= flag for the computation of frozen part of Born Effective Charges (PAW only)
!!  pawpiezo= flag for the computation of frozen part of Piezoelectric tensor (PAW only)
!!  prtbbb=if 1, print the band-by-band decomposition
!!  prtvol=print volume
!!  qphon(3)=phonon wavelength, in reduced coordinates
!!  qzero=1 if zero phonon wavevector
!!  rfdir(3)=defines the directions for the perturbations
!!  rfpert(mpert)=defines the perturbations
!!  rfphon=if 1, there are phonon perturbations
!!  rfstrs=if 1,2,3 there are strain perturbations
!!  typat(natom)=integer label of each type of atom (1,2,...)
!!  usepaw=1 if PAW, 0 otherwise
!!  usevdw= flag set to 1 if vdw DFT-D semi-empirical potential is in use
!!  zion(ntypat)=charge corresponding to the atom type
!!
!! SIDE EFFECTS
!!  d2cart_bbb(2,3,3,mpert,mband,mband*prtbbb)
!!
!!
!! NOTES
!! This routine is called only by the processor me==0 .
!! In consequence, no use of message and wrtout routine.
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_dyout(becfrnl,berryopt,blkflg,carflg,ddboun,ddkfil,dyew,dyfrlo,dyfrnl,&
& dyfrx1,dyfrx2,dyfr_cplex,dyfr_nondiag,dyvdw,d2cart,d2cart_bbb,&
& d2eig0,d2k0,d2lo,d2loc0,d2matr,d2nl,d2nl0,d2nl1,d2ovl,d2vn,&
& eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,eltvdw,&
& has_full_piezo,has_allddk,iout,mband,mpert,natom,ntypat,&
& outd2,pawbec,pawpiezo,piezofrnl,prtbbb,prtvol,qphon,qzero,typat,rfdir,&
& rfpert,rfphon,rfstrs,usepaw,usevdw,zion)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_dyout'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: berryopt,ddboun,dyfr_cplex,dyfr_nondiag,iout,mband,mpert
 integer,intent(in) :: natom,ntypat,outd2,pawbec,pawpiezo,prtbbb,prtvol,qzero
 integer, intent(in) :: rfphon,rfstrs,usepaw,usevdw
!arrays
 integer,intent(in) :: blkflg(3,mpert,3,mpert),carflg(3,mpert,3,mpert)
 integer,intent(in) :: ddkfil(3),rfdir(3),rfpert(mpert),typat(natom)
 real(dp),intent(in) :: becfrnl(3,natom,3*pawbec)
 real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert),d2eig0(2,3,mpert,3,mpert)
 real(dp),intent(in) :: d2k0(2,3,mpert,3,mpert),d2lo(2,3,mpert,3,mpert)
 real(dp),intent(in) :: d2loc0(2,3,mpert,3,mpert),d2matr(2,3,mpert,3,mpert)
 real(dp),intent(in) :: d2nl(2,3,mpert,3,mpert),d2nl0(2,3,mpert,3,mpert)
 real(dp),intent(in) :: d2nl1(2,3,mpert,3,mpert),d2ovl(2,3,mpert,3,mpert*usepaw)
 real(dp),intent(in) :: d2vn(2,3,mpert,3,mpert)
 real(dp),intent(in) :: dyew(2,3,natom,3,natom),dyfrlo(3,3,natom)
 real(dp),intent(in) :: dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
 real(dp),intent(in) :: dyfrx1(2,3,natom,3,natom),dyfrx2(3,3,natom)
 real(dp),intent(in) :: dyvdw(2,3,natom,3,natom*usevdw)
 real(dp),intent(in) :: eltcore(6,6),elteew(6+3*natom,6)
 real(dp),intent(in) :: eltfrhar(6,6),eltfrkin(6,6),eltfrloc(6+3*natom,6)
 real(dp),intent(in) :: eltfrnl(6+3*natom,6),eltfrxc(6+3*natom,6)
 real(dp),intent(in) :: eltvdw(6+3*natom,6*usevdw),piezofrnl(6,3*pawpiezo)
 real(dp),intent(in) :: qphon(3),zion(ntypat)
 real(dp),intent(inout) :: d2cart_bbb(2,3,3,mpert,mband,mband*prtbbb)
 logical,intent(in) :: has_allddk,has_full_piezo

!Local variables -------------------------
!scalars
 integer :: iband,idir1,idir2,ii,ipert1,ipert2,jj,nelmts,nline
 real(dp) :: qptnrm,zi,zr
!arrays
 real(dp) :: delta(3,3)

! *********************************************************************

!DEBUG
!write(std_out,*)' dfpt_dyout : enter '
!write(std_out,*)ddkfil
!ENDDEBUG

!Long print : includes detail of every part of the 2nd-order energy
 if(prtvol>=10)then

!  In case of phonon
   if (rfphon==1)then

!    write the Ewald part of the dynamical matrix
     write(iout,*)' '
     write(iout,*)' Ewald part of the dynamical matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           write(iout,*)' '
           do ipert2=1,natom
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 dyew(1,idir1,ipert1,idir2,ipert2),&
&                 dyew(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the local frozen wf part
     write(iout,*)' '
     write(iout,*)' Frozen wf local part of the dynamical matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           write(iout,*)' '
           do ipert2=1,natom
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 if(ipert1==ipert2)then
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   dyfrlo(idir1,idir2,ipert2),zero
                 else
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   zero,zero
                 end if
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the nonlo frozen wf part
     write(iout,*)' '
     write(iout,*)' Frozen wf non-local part of the dynamical matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           write(iout,*)' '
           do ipert2=1,natom
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 if(ipert1==ipert2.or.dyfr_nondiag==1)then
                   if (dyfr_cplex==1) then
                     write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                     dyfrnl(1,idir1,idir2,ipert1,1+(ipert2-1)*dyfr_nondiag),zero
                   else
                     write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                     dyfrnl(:,idir1,idir2,ipert1,1+(ipert2-1)*dyfr_nondiag)
                   end if
                 else
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   zero,zero
                 end if
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the nonlinear xc core correction(1) part
     write(iout,*)' '
     write(iout,*)' Frozen wf xc core (1) part',&
&     ' of the dynamical matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           write(iout,*)' '
           do ipert2=1,natom
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 dyfrx1(1,idir1,ipert1,idir2,ipert2),&
&                 dyfrx1(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the nonlinear xc core correction(2) part
     write(iout,*)' '
     write(iout,*)' Frozen wf xc core (2) part',&
&     ' of the dynamical matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           write(iout,*)' '
           do ipert2=1,natom
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 if(ipert1==ipert2)then
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   dyfrx2(idir1,idir2,ipert2),zero
                 else
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   zero,zero
                 end if
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the DFT-D vdw part of the dynamical matrix
     if (usevdw==1) then
       write(iout,*)' '
       write(iout,*)' DFT-D van der Waals part of the dynamical matrix'
       write(iout,*)'    j1       j2             matrix element'
       write(iout,*)' dir pert dir pert     real part   imaginary part'
       do ipert1=1,natom
         do idir1=1,3
           if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&           .or.   outd2==1                           )then
             write(iout,*)' '
             do ipert2=1,natom
               do idir2=1,3
                 if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   dyvdw(1,idir1,ipert1,idir2,ipert2),&
&                   dyvdw(2,idir1,ipert1,idir2,ipert2)
                 end if
               end do
             end do
           end if
         end do
       end do
     end if

!    End of the phonon condition
   end if

!  In case of atom. strain/electric field perturbation (piezoelectric tensor)
   if (pawpiezo==1.and.(rfpert(natom+2)==1.or.rfstrs/=0).and.outd2==1)then
     write(iout,*)' '
     write(iout,*)' Frozen wf part of the piezoelectric tensor'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     ipert1=natom+2
     do idir1=1,3
       write(iout,*)' '
       ii=1
       do ipert2=natom+3,natom+4
         do idir2=1,3
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           piezofrnl(ii,idir1),zero
           ii=ii+1
         end do
       end do
     end do
   end if

!  In case of atom. displ/electric field perturbation (Born Effective Charges)
   if (pawbec==1.and.(rfpert(natom+2)==1.or.rfphon==1).and.outd2==1)then
     write(iout,*)' '
     write(iout,*)' Frozen wf part of the Born Effective Charges'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     ipert1 = natom+2
     do ipert2=1,natom
       do idir1=1,3
         write(iout,*)' '
         do idir2=1,3
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           becfrnl(idir2,ipert2,idir1),zero
         end do
       end do
     end do
   end if

!  In case of strain
   if (rfstrs/=0)then

!    Write the Ewald part of the elastic tensor
     write(iout,*)' '
     write(iout,*)' Ewald part of the elastic tensor in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=natom+3,natom+4
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+3*(ipert1-natom-3)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 elteew(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Write the Ewald part of the internal strain coupling parameters
     write(iout,*)' '
     write(iout,*)' Ewald part of the internal strain coupling parameters'
     write(iout,*)'  (cartesian strain, reduced atomic coordinates)'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+6+3*(ipert1-1)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 elteew(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the local frozen wf part
     write(iout,*)' '
     write(iout,*)' Frozen wf local part of the elastic tensor in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=natom+3,natom+4
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+3*(ipert1-natom-3)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrloc(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

     write(iout,*)' '
     write(iout,*)' Frozen wf local part of the internal strain coupling parameters '
     write(iout,*)'  (cartesian strain, reduced atomic coordinates)'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+6+3*(ipert1-1)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrloc(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the nonlo frozen wf part
     write(iout,*)' '
     write(iout,*)' Frozen wf nonlocal part of the elastic tensor in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=natom+3,natom+4
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+3*(ipert1-natom-3)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrnl(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

     write(iout,*)' '
     write(iout,*)' Frozen wf nonlocal part of the internal strain coupling parameters '
     write(iout,*)'  (cartesian strain, reduced atomic coordinates)'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+6+3*(ipert1-1)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrnl(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the xc part
     write(iout,*)' '
     write(iout,*)' Frozen wf xc part of the elastic tensor in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=natom+3,natom+4
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+3*(ipert1-natom-3)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrxc(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

     write(iout,*)' '
     write(iout,*)' Frozen wf xc part of the internal strain coupling parameters '
     write(iout,*)'  (cartesian strain, reduced atomic coordinates)'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+6+3*(ipert1-1)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrxc(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the kinetic frozen wf part
     write(iout,*)' '
     write(iout,*)' Frozen wf kinetic part of the elastic tensor in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=natom+3,natom+4
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+3*(ipert1-natom-3)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrkin(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the hartree frozen wf part
     write(iout,*)' '
     write(iout,*)' Frozen wf hartree part of the elastic tensor in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=natom+3,natom+4
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+3*(ipert1-natom-3)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrhar(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the psp core part
     write(iout,*)' '
     write(iout,*)' Psp core part of the elastic tensor in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=natom+3,natom+4
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+3*(ipert1-natom-3)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltcore(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the DFT-D vdw part
     if (usevdw==1) then
       write(iout,*)' '
       write(iout,*)' DFT-D van der Waals part of the elastic tensor in cartesian coordinates'
       write(iout,*)'    j1       j2             matrix element'
       write(iout,*)' dir pert dir pert     real part   imaginary part'
       do ipert1=natom+3,natom+4
         do idir1=1,3
           if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&           .or.   outd2==1                           )then
             ii=idir1+3*(ipert1-natom-3)
             write(iout,*)' '
             do ipert2=natom+3,natom+4
               do idir2=1,3
                 if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                   jj=idir2+3*(ipert2-natom-3)
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   eltvdw(ii,jj),zero
                 end if
               end do
             end do
           end if
         end do
       end do

       write(iout,*)' '
       write(iout,*)' DFT-D2 van der Waals part of the internal strain coupling parameters'
       write(iout,*)'  (cartesian strain, reduced atomic coordinates)'
       write(iout,*)'    j1       j2             matrix element'
       write(iout,*)' dir pert dir pert     real part   imaginary part'
       do ipert1=1,natom
         do idir1=1,3
           if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&           .or.   outd2==1                           )then
             ii=idir1+6+3*(ipert1-1)
             write(iout,*)' '
             do ipert2=natom+3,natom+4
               do idir2=1,3
                 if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                   jj=idir2+3*(ipert2-natom-3)
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   eltvdw(ii,jj),zero
                 end if
               end do
             end do
           end if
         end do
       end do
     end if ! usevdw

!    End of the strain condition
   end if

!  Now the local nonstationary nonfrozenwf part
   if (outd2==1)then
     write(iout,*)' '
     write(iout,*)' Non-stationary local part of the 2-order matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if ((ipert1<=natom .or.&
&         (ipert1==natom+2.and.qzero==1.and.ddkfil(idir1)/=0)).or.&
&         ((ipert1==natom+3.or.ipert1==natom+4).and.&
&         (rfpert(natom+3)==1.or.rfpert(natom+4)==1)))then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2lo(1,idir1,ipert1,idir2,ipert2),&
&                 d2lo(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the nonlocal nonstationary nonfrozenwf part
   if (outd2==1)then
     write(iout,*)' '
     write(iout,*)' Non-stationary non-local part of the 2nd-order matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if ((ipert1<=natom .or.&
&         (ipert1==natom+2.and.qzero==1.and.ddkfil(idir1)/=0)).or.&
&         ((ipert1==natom+3.or.ipert1==natom+4).and.&
&         (rfpert(natom+3)==1.or.rfpert(natom+4)==1)))then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2nl(1,idir1,ipert1,idir2,ipert2),&
&                 d2nl(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the overlap change nonstationnary nonfrozenwf part (PAW only)
   if (outd2==1.and.usepaw==1)then
     write(iout,*)' '
     write(iout,*)' PAW: Non-stationary WF-overlap part of the 2nd-order matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if ((ipert1<=natom .or.&
&         (ipert1==natom+2.and.qzero==1.and.ddkfil(idir1)/=0)).or.&
&         ((ipert1==natom+3.or.ipert1==natom+4).and.&
&         (rfpert(natom+3)==1.or.rfpert(natom+4)==1)))then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2ovl(1,idir1,ipert1,idir2,ipert2),&
&                 d2ovl(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the 0-order local stationary nonfrozenwf part
   if (outd2==2)then
     write(iout,*)' '
     write(iout,*)' Stationary 0-order local part of the 2nd-order matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if(rfpert(ipert1)==1.and.rfdir(idir1)==1)then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2loc0(1,idir1,ipert1,idir2,ipert2),&
&                 d2loc0(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the stationary 0-order kinetic nonfrozenwf part
   if (outd2==2)then
     write(iout,*)' '
     write(iout,*)' Stationary 0-order kinetic part of the 2nd-order matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if(rfpert(ipert1)==1.and.rfdir(idir1)==1)then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2k0(1,idir1,ipert1,idir2,ipert2),&
&                 d2k0(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the stationary 0-order eigenvalue nonfrozenwf part
   if (outd2==2)then
     write(iout,*)' '
     write(iout,*)' Stationary 0-order eigenvalue part of the'&
&     ,' 2nd-order matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if(rfpert(ipert1)==1.and.rfdir(idir1)==1)then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2eig0(1,idir1,ipert1,idir2,ipert2),&
&                 d2eig0(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the stationary potential-density nonfrozenwf part
   if (outd2==2)then
     write(iout,*)' '
     write(iout,*)' Station. potential-density part of the ',&
&     ' 2nd-order matrix'
     write(iout,*)'  (Note : include some xc core-correction) '
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if(rfpert(ipert1)==1.and.rfdir(idir1)==1)then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2vn(1,idir1,ipert1,idir2,ipert2),&
&                 d2vn(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the stationary 0-order nonloc nonfrozenwf part
   if (outd2==2)then
     write(iout,*)' '
     write(iout,*)' Stationary 0-order nonlocal part of the 2-order'&
&     ,' matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if(rfpert(ipert1)==1.and.rfdir(idir1)==1)then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2nl0(1,idir1,ipert1,idir2,ipert2),&
&                 d2nl0(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the stationary 1-order nonloc nonfrozenwf part
   if (outd2==2)then
     write(iout,*)' '
     write(iout,*)' Stationary 1-order nonlocal part of the'&
&     ,' 2nd-order matrix'
     write(iout,*)' (or the ddk wf part of it, in case of',&
&     ' an electric field perturbation )'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if(rfpert(ipert1)==1.and.rfdir(idir1)==1)then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2nl1(1,idir1,ipert1,idir2,ipert2),&
&                 d2nl1(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  End of the long print out condition
 end if


!Derivative database initialisation

!Calculation of the number of elements
 nelmts=0
 do ipert1=1,mpert
   do idir1=1,3
     do ipert2=1,mpert
       do idir2=1,3
         nelmts=nelmts+blkflg(idir1,ipert1,idir2,ipert2)
       end do
     end do
   end do
 end do

 if(outd2==2)then
   write(ddboun, '(/,a,i8)' ) ' 2nd derivatives (stationary) - # elements :',nelmts
 else if(outd2==1)then
   write(ddboun, '(/,a,i8)' ) ' 2nd derivatives (non-stat.)  - # elements :',nelmts
 end if

!Phonon wavevector
 qptnrm=1.0_dp

!Note : if qptnrm should assume another value, it should
!be checked if the f6.1 format is OK.
 write(ddboun, '(a,3es16.8,f6.1)' ) ' qpt',(qphon(ii),ii=1,3),qptnrm

!Now the whole 2nd-order matrix, but not in cartesian coordinates,
!and masses not included
 write(iout,*)' '
 write(iout,*)' 2nd-order matrix (non-cartesian coordinates,',' masses not included,'
 write(iout,*)'  asr not included )'
 if(rfstrs/=0) then
   write(iout,*)' cartesian coordinates for strain terms (1/ucvol factor '
   write(iout,*)'  for elastic tensor components not included) '
 end if
 write(iout,*)'    j1       j2             matrix element'
 write(iout,*)' dir pert dir pert     real part     imaginary part'
 nline=1
 do ipert1=1,mpert
   do idir1=1,3
     if(nline/=0)write(iout,*)' '
     nline=0
     do ipert2=1,mpert
       do idir2=1,3
         if(blkflg(idir1,ipert1,idir2,ipert2)==1)then
           nline=nline+1
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           d2matr(1,idir1,ipert1,idir2,ipert2),&
&           d2matr(2,idir1,ipert1,idir2,ipert2)
           write(ddboun, '(4i4,2d22.14)' )idir1,ipert1,idir2,ipert2,&
&           d2matr(1,idir1,ipert1,idir2,ipert2),&
&           d2matr(2,idir1,ipert1,idir2,ipert2)
         end if
       end do
     end do
   end do
 end do

!Now the dynamical matrix
 if(rfphon==1)then
   write(iout,*)' '
   write(iout,*)' Dynamical matrix, in cartesian coordinates,'
   write(iout,*)'  if specified in the inputs, asr has been imposed'
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   nline=1
   do ipert1=1,natom
     do idir1=1,3
       if(nline/=0)write(iout,*)' '
       nline=0
       do ipert2=1,natom
         do idir2=1,3
           if(carflg(idir1,ipert1,idir2,ipert2)==1)then
             nline=nline+1
             write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&             d2cart(1,idir1,ipert1,idir2,ipert2),&
&             d2cart(2,idir1,ipert1,idir2,ipert2)
           end if
         end do
       end do
     end do
   end do
 end if

!Now the dielectric tensor ! normal case
 if(rfpert(natom+2)==1)then

   write(iout,*)' '
   write(iout,*)' Dielectric tensor, in cartesian coordinates,'
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   ipert1=natom+2
   ipert2=natom+2
   nline=1
   do idir1=1,3
     if(nline/=0)write(iout,*)' '
     nline=0
     do idir2=1,3
       if(carflg(idir1,ipert1,idir2,ipert2)==1)then
         nline=nline+1
         write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&         d2cart(1,idir1,ipert1,idir2,ipert2),&
&         d2cart(2,idir1,ipert1,idir2,ipert2)
       end if
     end do
   end do

   if (prtbbb == 1) then

     delta(:,:) = zero
     delta(1,1) = one ; delta(2,2) = one ; delta(3,3) = one

     write(iout,*)
     write(iout,*)'Band by band decomposition of the dielectric tensor'
     write(iout,*)' '

     write(iout,*)' Vacuum polarization'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part    imaginary part'
     nline=1
     do idir1=1,3
       if(nline/=0)write(iout,*)' '
       nline=0
       do idir2=1,3
         nline=nline+1
         write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&         delta(idir2,idir1),zero
       end do
     end do

     do iband = 1,mband
       write(iout,*)' '
       write(iout,*)' Dielectric tensor, in cartesian coordinates, for band',iband
       write(iout,*)'    j1       j2             matrix element'
       write(iout,*)' dir pert dir pert     real part    imaginary part'
       ipert1 = natom + 2
       ipert2 = natom + 2
       nline=1
       do idir1=1,3
         if(nline/=0)write(iout,*)' '
         nline=0
         do idir2=1,3
           if(carflg(idir1,ipert1,idir2,ipert2)==1)then
!            substract vacuum polarization
             if (idir1 == idir2) then
               d2cart_bbb(1,idir1,idir2,ipert2,iband,iband) = &
&               d2cart_bbb(1,idir1,idir2,ipert2,iband,iband) - 1
             end if
             nline=nline+1
             write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&             d2cart_bbb(1,idir1,idir2,ipert2,iband,iband),&
&             d2cart_bbb(2,idir1,idir2,ipert2,iband,iband)
           end if
         end do
       end do
     end do !iband

   end if !prtbbb

 end if ! end natom+2 dielectric output

!Now the effective charges
!In case of the stationary calculation
 if(outd2==2 .and. rfpert(natom+2)==1 .and.rfphon==1)then
   write(iout,*)' '
   write(iout,*)' Effective charges, in cartesian coordinates,'
   write(iout,*)'  if specified in the inputs, asr has been imposed'
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   ipert1=natom+2
   nline=1
   do idir1=1,3
     if(nline/=0)write(iout,*)' '
     nline=0
     do ipert2=1,natom
       do idir2=1,3
         if(carflg(idir1,ipert1,idir2,ipert2)==1)then
           nline=nline+1
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           d2cart(1,idir1,ipert1,idir2,ipert2),&
&           d2cart(2,idir1,ipert1,idir2,ipert2)
         end if
       end do
     end do
   end do
 end if

!Now in case of the non-stationary calculation
 if(outd2==1 .and. rfpert(natom+2)==1)then
   write(iout,*)' '
   if(usepaw==1.and..not.(has_allddk))then
     write(iout,*)' Warning: Born effectives charges are not correctly computed'
     write(iout,*)' you need all ddk perturbations!'
   end if
   write(iout,*)' Effective charges, in cartesian coordinates,'
   write(iout,*)' (from electric field response) '
   write(iout,*)'  if specified in the inputs, asr has been imposed'
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   ipert2=natom+2
   nline=1
   do idir2=1,3
     if(nline/=0)write(iout,*)' '
     nline=0
     do ipert1=1,natom
       do idir1=1,3
         if(carflg(idir1,ipert1,idir2,ipert2)==1)then
           nline=nline+1
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           d2cart(1,idir1,ipert1,idir2,ipert2),&
&           d2cart(2,idir1,ipert1,idir2,ipert2)
         end if
       end do
     end do
   end do
 end if

 if(outd2==1 .and. rfphon==1 .and. qzero==1&
& .and. ( (ddkfil(1)/=0.or.ddkfil(2)/=0.or.ddkfil(3)/=0) .or.   &
& berryopt==4 .or. berryopt==6 .or. berryopt==7 .or. berryopt==14 .or. berryopt==16 .or. berryopt==17 ) )then  !!HONG  need to test for fixed E and D
   write(iout,*)' '
   if(usepaw==1.and..not.(has_allddk))then
     write(iout,*)' Warning: Born effectives charges are not correctly computed'
     write(iout,*)' you need all ddk perturbations!'
   end if
   write(iout,*)' Effective charges, in cartesian coordinates,'
   write(iout,*)' (from phonon response) '
   write(iout,*)'  if specified in the inputs, asr has been imposed'
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   nline=1
   do ipert2=1,natom
     do idir2=1,3
       if(nline/=0)write(iout,*)' '
       nline=0
       ipert1=natom+2
       do idir1=1,3
         if(carflg(idir1,ipert1,idir2,ipert2)==1)then
           nline=nline+1
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           d2cart(1,idir1,ipert1,idir2,ipert2),&
&           d2cart(2,idir1,ipert1,idir2,ipert2)
         end if
       end do
     end do
   end do
   write(iout,*)' '
   write(iout,*)' '
   write(iout,*)' '

   if (prtbbb == 1) then

     write(iout,*)'Band by band decomposition of the Born effective charges'
     write(iout,*)' '
     write(iout,*)'Ionic charges in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part    imaginary part'
     zr = zero
     zi = zero
     do ipert2=1,natom
       do idir2=1,3
         if(nline/=0)write(iout,*)' '
         nline=0
         ipert1=natom+2
         do idir1=1,3
           zr = zero
           if (idir1 == idir2) then
             zr = zion(typat(ipert2))
           end if
           nline=nline+1
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           zr,zi
         end do
       end do
     end do

     do iband = 1,mband
       write(iout,*)' '
       write(iout,*)' Effective charges, in cartesian coordinates, for band',iband
       write(iout,*)' (from phonon response) '
       write(iout,*)'  if specified in the inputs, asr has been imposed'
       write(iout,*)'    j1       j2             matrix element'
       write(iout,*)' dir pert dir pert     real part    imaginary part'
       nline=1
       do ipert2=1,natom
         do idir2=1,3
           if(nline/=0)write(iout,*)' '
           nline=0
           ipert1=natom+2
           do idir1=1,3
             if(carflg(idir1,ipert1,idir2,ipert2)==1)then
               nline=nline+1
               write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&               d2cart_bbb(1,idir1,idir2,ipert2,iband,iband),&
&               d2cart_bbb(2,idir1,idir2,ipert2,iband,iband)
             end if
           end do
         end do
       end do
     end do !iband
   end if !prtbbb
 end if ! end of print effective charges

!Now the elastic tensor
 if(rfstrs/=0) then
   write(iout,*)' '
   write(iout,*)' Rigid-atom elastic tensor , in cartesian coordinates,'
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   nline=1
   do ipert1=natom+3,natom+4
     do idir1=1,3
       if(nline/=0)write(iout,*)' '
       nline=0
       do ipert2=natom+3,natom+4
         do idir2=1,3
           if(carflg(idir1,ipert1,idir2,ipert2)==1)then
             nline=nline+1
             write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&             d2cart(1,idir1,ipert1,idir2,ipert2),&
&             d2cart(2,idir1,ipert1,idir2,ipert2)
           end if
         end do
       end do
     end do
   end do
 end if

!Now the internal strain coupling parameters
 if(rfstrs/=0) then
   write(iout,*)' '
   write(iout,*)' Internal strain coupling parameters, in cartesian coordinates,'
   write(iout,*)'  zero average net force deriv. has been imposed  '
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   nline=1
   do ipert1=1,natom
     do idir1=1,3
       if(nline/=0)write(iout,*)' '
       nline=0
       do ipert2=natom+3,natom+4
         do idir2=1,3
           if(carflg(idir1,ipert1,idir2,ipert2)==1)then
             nline=nline+1
             write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&             d2cart(1,idir1,ipert1,idir2,ipert2),&
&             d2cart(2,idir1,ipert1,idir2,ipert2)
           end if
         end do
       end do
     end do
   end do
 end if

!Now the piezoelectric tensor
 if(rfstrs/=0 .and. (ddkfil(1)/=0.or.ddkfil(2)/=0.or.ddkfil(3)/=0))then
   write(iout,*)' '
   if(usepaw==1.and..not.(has_allddk))then
     write(iout,*)' Warning: Rigid-atom proper piezoelectric tensor is not correctly computed'
     write(iout,*)' you need all ddk perturbations!'
   end if
   write(iout,*)' Rigid-atom proper piezoelectric tensor, in cartesian coordinates,'
   write(iout,*)' (from strain response)'
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   nline=1
   ipert1=natom+2
   do idir1=1,3
     if(nline/=0)write(iout,*)' '
     nline=0
     do ipert2=natom+3,natom+4
       do idir2=1,3
         if(carflg(idir1,ipert1,idir2,ipert2)==1)then
           nline=nline+1
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           d2cart(1,idir1,ipert1,idir2,ipert2),&
&           d2cart(2,idir1,ipert1,idir2,ipert2)
         end if
       end do
     end do
   end do
 end if

!Now the piezoelectric tensor
 if(outd2==1 .and. (pawpiezo==1.and.rfpert(natom+2)==1)&
& .and. (ddkfil(1)/=0.or.ddkfil(2)/=0.or.ddkfil(3)/=0)) then
   write(iout,*)' '
   if(usepaw==1.and..not.(has_allddk))then
     write(iout,*)' Warning: Rigid-atom proper piezoelectric tensor is not correctly computed'
     write(iout,*)' you need all ddk perturbations!'
   end if
   if(usepaw==1.and..not.has_full_piezo)then
     write(iout,*)' Warning: The rigid-atom proper piezoelectric tensor'
     write(iout,*)' from  electric field response requires nsym=1'
   end if
   if (has_full_piezo) then
     write(iout,*)' Rigid-atom proper piezoelectric tensor, in cartesian coordinates,'
     write(iout,*)' (from electric field response)'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part    imaginary part'
     nline=1
     ipert1=natom+2
     do idir1=1,3
       if(nline/=0)write(iout,*)' '
       nline=0
       do ipert2=natom+3,natom+4
         do idir2=1,3
           if(carflg(idir2,ipert2,idir1,ipert1)==1)then
             nline=nline+1
             write(iout,'(2(i4,i5),2(1x,f20.10))')idir2,ipert2,idir1,ipert1,&
&             d2cart(1,idir2,ipert2,idir1,ipert1),&
&             d2cart(2,idir2,ipert2,idir1,ipert1)
           end if
         end do
       end do
     end do
   end if
 end if

 !write(std_out,*)' dfpt_dyout : exit '

end subroutine dfpt_dyout
!!***

! LocalWords:  piezo
