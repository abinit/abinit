!!****m* ABINIT/m_read_plowannier
!! NAME
!!  m_read_plowannier
!!
!! FUNCTION
!!  Read Wannier coefficient in the file forlb.ovlp for ucrpa calculation
!!  this file was typically created in a DFT run with usedmft=1 and nbandkss -1
!!
!! COPYRIGHT
!! Copyright (C) 2006-2020 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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

MODULE m_read_plowannier

 use defs_basis
 use m_abicore
 use m_errors

 use m_io_tools,      only : open_file
 use m_crystal,       only : crystal_t
 use m_bz_mesh,       only : kmesh_t, get_BZ_item
 use m_pawang,        only : pawang_type

 implicit none

 private

 public :: read_plowannier
!!***

contains

!!****m* m_read_plowannier/read_plowannier
!! NAME
!! read_plowannier
!!
!! FUNCTION
!!  Read Wannier coefficient in the file forlb.ovlp for ucrpa calculation
!!  this file was typically created in a DFT run with usedmft=1 and nbandkss -1
!!
!! COPYRIGHT
!! Copyright (C) 2006-2020 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!!
!! INPUTS
!! Cryst<cryst_t>= data type gathering info on symmetries and unit cell
!!    %natom=number of atoms
!!    %nsym=number of symmetries
!!    %xred(3,natom)=reduced coordinated of atoms
!!    %typat(natom)=type of each atom
!!    %rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!    %timrev= 2 if time reversal can be used, 1 otherwise
!! Kmesh <kmesh_t>= datatype gathering parameters related to the k-point sampling
!!    %nibz=number of k-points in the IBZ
!!    %nbz=number of k-points in the BZ
!!    %bz(3,nbz)=reduced coordinates for k-points in the full Brillouin zone
!!    %ibz(3,nibz)=reduced coordinates for k-points in the irreducible wedge
!!    %tab(nbz)=mapping between a kpt in the BZ (array bz) and the irred point in the array ibz
!!    %tabi(nbz)= -1 if inversion is needed to obtain this particular kpt in the BZ, 1 means identity
!!    %tabo(nbz)= for each point in the BZ, the index of the symmetry operation S in reciprocal
!!      space which rotates k_IBZ onto \pm k_BZ (depending on tabi)
!!    %tabp(nbz)= For each k_BZ, it gives the phase factors associated to non-symmorphic operations, i.e
!!      e^{-i 2 \pi k_IBZ \cdot R{^-1}t} == e{-i 2\pi k_BZ cdot t} where :
!!      \transpose R{-1}=S and (S k_IBZ) = \pm k_BZ (depending on ktabi)
!!    %tabr(nfftot,nbz) For each point r on the real mesh and for each k-point in the BZ, tabr
!!      gives the index of (R^-1 (r-t)) in the FFT array where R=\transpose S^{-1} and k_BZ=S k_IBZ.
!!      t is the fractional translation associated to R
!!  luwindow: T if ucrpa_window is activated, F if not.
!!  prtvol: integer to give the amount of printing.
!! nsppol : number of spin polarization.
!! Pawang<pawang_type> angular mesh discretization and related data:
!!
!! OUTPUT
!! bandinf, bandsup : lower and upper bands for define Wannier functions
!! coeffW_BZ(nsppol,bandinf:bandsup,nvz,2*lcor+1)
!! lcor : angular momentum for correlated orbitals
!! itypatcor : correlated species
!!
!! PARENTS
!!      cchi0,cchi0q0,prep_calc_ucrpa
!!
!! CHILDREN
!!      get_bz_item,wrtout
!!
!! SOURCE

subroutine read_plowannier(cryst,bandinf,bandsup,coeffW_BZ,itypatcor,Kmesh,lcor,luwindow,nspinor,nsppol,pawang,prtvol,ucrpa_bands)

!Arguments ------------------------------------
!types and arrays
 type(kmesh_t),intent(in) :: Kmesh
 type(crystal_t),intent(in) :: Cryst
 complex(dpc), allocatable, intent(inout) :: coeffW_BZ(:,:,:,:,:,:)
 type(Pawang_type),intent(in) :: Pawang
!scalars
 logical, intent(inout) :: luwindow
 integer, intent(out) :: bandinf,bandsup,itypatcor,lcor
 integer, intent(in) :: nspinor,nsppol,prtvol
 integer, intent(in) :: ucrpa_bands(2)

!Local variables-------------------------------
 character(len=500) :: message,msg
 integer :: at_indx,ik_ibz,band1,m1,m2,spin,ik_bz,dummy,isym,itim,iat,indx,ispinor,unt
 real(dp) :: xx,yy
 real(dp) :: kbz(3)
 complex(dpc), allocatable :: coeffW_IBZ(:,:,:,:,:,:)
! *********************************************************************
 write(message,*) "Read wannier in iBZ"
 call wrtout(std_out,message,'COLL')

 if (open_file('forlb.ovlp',msg,newunit=unt,form='formatted',status='unknown') /= 0) then
   MSG_ERROR(msg)
 end if
 rewind(unt)
 read(unt,*) message
 read(unt,*) message, lcor,itypatcor
 read(unt,*) message, bandinf,bandsup
 write(std_out,*) 'read from forlb.ovlp',lcor, bandinf,bandsup
 write(std_out,*) "for wannier", bandinf,bandsup
 if(prtvol>0) then
 endif

 if(.not.luwindow.and.(ucrpa_bands(1)/=bandinf.or.ucrpa_bands(2)/=bandsup)) then
   write(msg,'(a,a)')' Bands used for Wannier construction and cRPA construction differ',&
& 'It might be physically correct and it is possible with the current implementation'
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
!   MSG_ERROR(msg)
 end if

!Do not dead the bandinf, bandinf redondance information

 ABI_ALLOCATE(coeffW_IBZ,(Cryst%nattyp(itypatcor),nsppol,bandinf:bandsup,Kmesh%nibz,nspinor,2*lcor+1))
 coeffW_IBZ=czero
 do spin=1,nsppol
   do ik_ibz=1,Kmesh%nibz
     !read k
      read(unt,*)
      do band1=bandinf,bandsup
     !read band
        read(unt,*)
     !read projection
        do ispinor=1,nspinor
          do iat=1,Cryst%nattyp(itypatcor)
            do m1=1,2*lcor+1
              read(unt,*) dummy,dummy,dummy,dummy,xx,yy
              coeffW_IBZ(iat,spin,band1,ik_ibz,ispinor,m1)=cmplx(xx,yy)
            end do
          end do
        end do
      end do
   end do
 end do
 close(unt)

 ABI_ALLOCATE(coeffW_BZ,(Cryst%nattyp(itypatcor),nsppol,bandinf:bandsup,Kmesh%nbz,nspinor,2*lcor+1))
 coeffW_BZ=czero

 if (Kmesh%nbz==Kmesh%nibz) then
   coeffW_BZ=coeffW_IBZ
 else if (Cryst%nsym==1) then
   write(message,*) "Reconstruct in full BZ"
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   do ik_bz=1,Kmesh%nbz
!     if(prtvol>=10) write(6,*) "ik",ik_bz,Kmesh%tab(ik_bz),Kmesh%tabi(ik_bz),Kmesh%tabo(ik_bz)
     if(Kmesh%tabi(ik_bz)==1) then
       coeffW_BZ(:,:,:,ik_bz,:,:)=coeffW_IBZ(:,:,:,Kmesh%tab(ik_bz),:,:)
     else if(Kmesh%tabi(ik_bz)==-1) then
       coeffW_BZ(:,:,:,ik_bz,:,:)=conjg(coeffW_IBZ(:,:,:,Kmesh%tab(ik_bz),:,:))
     endif
!     if(prtvol>=10)write(6,*) "coeffW 21",coeffW_BZ(1,1,bandinf,ik_bz,:)
!     if(prtvol>=10)write(6,*) "coeffW 21",coeffW_BZ(1,1,bandinf+1,ik_bz,:)
!     if(prtvol>=10)write(6,*) "coeffW 25",coeffW_BZ(1,1,bandsup,ik_bz,:)
   enddo
 else if (Cryst%nsym>1) then
   write(message,*) "Reconstruct in full BZ (nsym=0)"
   call wrtout(std_out,message,'COLL')
   do ik_bz=1,Kmesh%nbz
   !write(6,*) ik_bz,Kmesh%nbz
     call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym,itim)
     do indx=1,cryst%nattyp(itypatcor)
       iat=cryst%atindx1(indx) ! correct index for the full atom list
       at_indx=cryst%indsym(4,isym,iat) !! see eg sym_matlu and m_crystal
       do spin=1,nsppol
         do ispinor=1,nspinor
           do m1=1,2*lcor+1
             do m2=1,2*lcor+1
             coeffW_BZ(indx,spin,:,ik_bz,ispinor,m1)= coeffW_BZ(indx,spin,:,ik_bz,ispinor,m1)&
&                     +coeffW_IBZ(at_indx,spin,:,ik_ibz,ispinor,m2)*pawang%zarot(m2,m1,lcor+1,isym)
             enddo
!             write(6,'(20f7.3)') (pawang%zarot(m1,m2,lcor+1,isym),m2=1,2*lcor+1)
           enddo
         enddo
       enddo
     enddo
!     do m1=1,2*lcor+1
!      write(6,*) "coeffW_IBZ",coeffW_IBZ(1,8,ik_ibz,m1)
!     enddo
!     do m1=1,2*lcor+1
!      write(6,*) "coeffW_ BZ",coeffW_BZ(1,8,ik_bz,m1)
!     enddo
   enddo
 end if
 ABI_DEALLOCATE(coeffW_IBZ)


end subroutine read_plowannier
!!***

END MODULE m_read_plowannier
!!***
