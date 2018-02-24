!{\src2tex{textfont=tt}}
!!****f* ABINIT/mag_constr
!! NAME
!! mag_constr
!!
!! FUNCTION
!! This routine is called to compute the potential corresponding to constrained magnetic moments.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (ILuk, MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  natom=number of atoms
!!  spinat=fixed magnetic moments vectors
!!  nspden = number of spin densities (1 2 or 4)
!!  magconon=constraining option (on/off); 1=fix only the direction, 2=fix the direction and size
!!  magcon_lambda=the size of the penalty terms
!!  rprimd=lattice vectors (dimensionful)
!!  mpi_enreg=mpi structure with communicator info
!!  nfft=number of points in standard fft grid
!!  ngfft=FFT grid dimensions
!!  ntypat=number of types of atoms
!!  ratsph=radii for muffin tin spheres of each atom
!!  rhor=density in real space
!!  typat=types of atoms
!!  xred=reduced atomic positions
!!
!! OUTPUT
!!  Vmagconstr=the constraining potential
!!
!! PARENTS
!!      energy,rhotov,setvtr
!!
!! CHILDREN
!!      calcdensph,metric,ptabs_fourdp,timab,xmpi_sum
!!
!! NOTES
!!  based on html notes for the VASP implementation at 
!!  http://cms.mpi.univie.ac.at/vasp/vasp/Constraining_direction_magnetic_moments.html
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mag_constr(natom,spinat,nspden,magconon,magcon_lambda,rprimd, &
                      mpi_enreg,nfft,ngfft,ntypat,ratsph,rhor, &
                      typat,Vmagconstr,xred)

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 
 use m_mpinfo,  only : ptabs_fourdp

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mag_constr'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_54_abiutil
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,magconon,nfft,nspden
 integer,intent(in) :: ntypat
 real(dp),intent(in) :: magcon_lambda
 real(dp),intent(out) :: Vmagconstr(nfft,nspden)
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in)  :: typat(natom)
 integer,intent(in)  :: ngfft(18)
 real(dp),intent(in) :: ratsph(ntypat)
 real(dp),intent(in) :: rhor(nfft,nspden)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: spinat(3,natom)
 real(dp),intent(in) :: xred(3,natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: ishift=5
 integer :: iatom, ierr
 integer :: cplex1=1
 integer :: n1a, n1b, n3a, n3b, n2a, n2b
 integer :: n1, n2, n3
 integer :: ifft_local
 integer ::  i1,i2,i3,ix,iy,iz,izloc,nd3
 real(dp) :: arg,intgden_proj,r2atsph,rr1,rr2,rr3,fsm,ratsm,ratsm2,difx,dify,difz,r2,rx,ry,rz
 real(dp) :: norm
 real(dp),parameter :: delta=0.99_dp
!arrays
 real(dp) :: intgden(nspden,natom)
 real(dp) :: spinat_norm(3,natom)
 !real(dp),allocatable :: cmm_x(:,:),cmm_y(:,:),cmm_z(:,:)
 real(dp):: cmm_x,cmm_y,cmm_z
 real(dp) :: gprimd(3,3),rmet(3,3),gmet(3,3),ucvol
 real(dp) :: tsec(2)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)

! ***********************************************************************************************

!We need the metric because it is needed in calcdensph.F90
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!We need the integrated magnetic moments and the smoothing function
 call calcdensph(gmet,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,std_out,ratsph,rhor,rprimd,typat,ucvol,xred,1,cplex1,intgden=intgden)

 n1 = ngfft(1)
 n2 = ngfft(2)
 n3 = ngfft(3)
 nd3=n3/mpi_enreg%nproc_fft

 ratsm = 0.05_dp ! default value for the smearing region radius - may become input variable later
 Vmagconstr = zero

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Loop over atoms
!-------------------------------------------
 do iatom=1,natom
   
   norm = sqrt(sum(spinat(:,iatom)**2))
   spinat_norm(:,iatom) = zero
   if (norm > tol10) then
     spinat_norm(:,iatom) = spinat(:,iatom) / norm
   else if (magconon == 1) then
!    if spinat = 0 and we are imposing the direction only, skip this atom
     cycle 
   end if
   
!  Define a "box" around the atom
   r2atsph=1.0000001_dp*ratsph(typat(iatom))**2
   rr1=sqrt(r2atsph*gmet(1,1))
   rr2=sqrt(r2atsph*gmet(2,2))
   rr3=sqrt(r2atsph*gmet(3,3))

   n1a=int((xred(1,iatom)-rr1+ishift)*n1+delta)-ishift*n1
   n1b=int((xred(1,iatom)+rr1+ishift)*n1      )-ishift*n1
   n2a=int((xred(2,iatom)-rr2+ishift)*n2+delta)-ishift*n2
   n2b=int((xred(2,iatom)+rr2+ishift)*n2      )-ishift*n2
   n3a=int((xred(3,iatom)-rr3+ishift)*n3+delta)-ishift*n3
   n3b=int((xred(3,iatom)+rr3+ishift)*n3      )-ishift*n3

   ratsm2 = -(ratsm**2 - 2*ratsph(typat(iatom))*ratsm)

   do i3=n3a,n3b
     iz=mod(i3+ishift*n3,n3)
     if(fftn3_distrib(iz+1)==mpi_enreg%me_fft) then
       izloc = ffti3_local(iz+1) - 1
       difz=dble(i3)/dble(n3)-xred(3,iatom)
       do i2=n2a,n2b
         iy=mod(i2+ishift*n2,n2)
         dify=dble(i2)/dble(n2)-xred(2,iatom)
         do i1=n1a,n1b
           ix=mod(i1+ishift*n1,n1)
           difx=dble(i1)/dble(n1)-xred(1,iatom)
           rx=difx*rprimd(1,1)+dify*rprimd(1,2)+difz*rprimd(1,3)
           ry=difx*rprimd(2,1)+dify*rprimd(2,2)+difz*rprimd(2,3)
           rz=difx*rprimd(3,1)+dify*rprimd(3,2)+difz*rprimd(3,3)
           r2=rx**2+ry**2+rz**2


!          Identify the fft indexes of the rectangular grid around the atom
           if (r2 > r2atsph) cycle

           fsm = radsmear(r2, r2atsph, ratsm2)

           ifft_local=1+ix+n1*(iy+n2*izloc)

!          Calculate the x- and y-components of the square bracket term
           cmm_x = zero
           cmm_y = zero
           cmm_z = zero
           intgden_proj = zero
           if (nspden == 4) then
             if (magconon==1) then
!              Calculate the scalar product of the fixed mag. mom. vector and calculated mag. mom. vector
!              This is actually the size of the projection of the calc. mag. mom. vector on the fixed mag. mom. vector
               intgden_proj=spinat_norm(1,iatom)*intgden(2,iatom)+ &
&               spinat_norm(2,iatom)*intgden(3,iatom)+ &
&               spinat_norm(3,iatom)*intgden(4,iatom)

               cmm_x=intgden(2,iatom)
               cmm_x=cmm_x-spinat_norm(1,iatom)*intgden_proj

               cmm_y=intgden(3,iatom)
               cmm_y=cmm_y-spinat_norm(2,iatom)*intgden_proj

             else if (magconon==2 .and. nspden == 4) then
               cmm_x=intgden(2,iatom)-spinat(1,iatom)
               cmm_y=intgden(3,iatom)-spinat(2,iatom)
             end if

!            do loop on spirals - should be added later on. arg = scalar product k.R (with 2pi factor?)
             arg=0

!            Calculate the constraining potential for x- and y- components of the mag. mom. vector
!            Eric Bousquet has derived the relationship between spin components and potential spin matrix elements:
!            1 = up up     = +z
!            2 = down down = -z
!            3 = up down   = +x
!            4 = down up   = -y
             Vmagconstr(ifft_local,3)=Vmagconstr(ifft_local,3) + &
&             2*magcon_lambda*fsm*cmm_x
!            & 2*magcon_lambda*fsm*(cmm_x*cos(arg)+cmm_y*sin(arg))
             Vmagconstr(ifft_local,4)=Vmagconstr(ifft_local,4) - &
&             2*magcon_lambda*fsm*cmm_y
!            & 2*magcon_lambda*fsm*(cmm_y*cos(arg)+cmm_x*sin(arg))
!            end do loop on spirals
           end if ! nspden 4

!          Calculate the z-component of the square bracket term
           if (magconon==1) then
             if (nspden == 4) then
               ! m_z - spinat_z * <m | spinat>
               cmm_z = intgden(4,iatom) - spinat_norm(3,iatom)*intgden_proj
             else if (nspden == 2) then
               ! this will be just a sign +/- : are we in the same direction as spinat_z?
               !    need something more continuous??? To make sure the gradient pushes the state towards FM/AFM?
               cmm_z = -sign(one, (intgden(1,iatom)-intgden(2,iatom))*spinat_norm(3,iatom))
             end if
           else if (magconon==2) then
             if (nspden == 4) then
               cmm_z=intgden(4,iatom)-spinat(3,iatom)
             else if (nspden == 2) then
               ! this is up spins - down spins - requested moment ~ 0
               ! EB: note that intgden comes from calcdensph, which, in nspden=2 case, returns
               ! intgden(1)=rho_up=n+m
               ! intgden(2)=rho_dn=n-m
               ! Then, is the following line be
               ! cmm_z=half*(intgden(1,iatom)-intgden(2,iatom)) - spinat(3,iatom) 
               ! ??
               cmm_z=intgden(1,iatom)-intgden(2,iatom) - spinat(3,iatom)
             end if
           end if

!          Calculate the constraining potential for z-component of the mag. mom. vector
           Vmagconstr(ifft_local,1)=Vmagconstr(ifft_local,1) + 2*magcon_lambda*fsm*cmm_z
           Vmagconstr(ifft_local,2)=Vmagconstr(ifft_local,2) - 2*magcon_lambda*fsm*cmm_z
!          end do loop on spirals

         end do  ! i1
       end do  ! i2 
     end if  ! if this is my fft slice
   end do ! i3

!  end loop over atoms
 end do

!MPI parallelization
!TODO: test if xmpi_sum does the correct stuff for a slice of Vmagconstr
 if(mpi_enreg%nproc_fft>1)then
   call timab(48,1,tsec)
   call xmpi_sum(Vmagconstr,mpi_enreg%comm_fft,ierr)
   call timab(48,2,tsec)
 end if

! write (201,*) '# potential 1'
! write (201,*) Vmagconstr(:,1)

! write (202,*) '# potential 2'
! write (202,*) Vmagconstr(:,2)

! if (nspden > 2) then
!   write (203,*) '# potential 3'
!   write (203,*) Vmagconstr(:,3)

!   write (204,*) '# potential 4'
!   write (204,*) Vmagconstr(:,4)
! end if

end subroutine mag_constr

!!***
