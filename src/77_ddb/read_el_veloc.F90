!{\src2tex{textfont=tt}}
!!****f* ABINIT/read_el_veloc
!!
!! NAME
!! read_el_veloc
!!
!! FUNCTION
!! This routine reads the velocities of the electronic GS
!! for all kpts and bands
!! then maps them into the FS kpt states
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (JPCroc) based on conducti
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nkpt_in = number of kpoints according to parent routine 
!! nband_in = number of bands according to parent routine
!! nsppol_in = number of spin polarizations
!!
!! OUTPUT
!! el_veloc(nkpt_in,nband_in,3)
!!  
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      destroy_kptrank,get_rank_1kpt,hdr_free,inpgkk,mkkptrank
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine read_el_veloc(nband_in,nkpt_in,kpt_in,nsppol_in,elph_tr_ds)

 use defs_basis
 use defs_abitypes
 use defs_elphon
 use m_errors
 use m_profiling_abi
 use m_kptrank
 use m_hdr

 use m_io_tools,        only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'read_el_veloc'
 use interfaces_72_response
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer, intent(in) :: nband_in,nkpt_in,nsppol_in
 type(elph_tr_type), intent(inout) :: elph_tr_ds
 real(dp), intent(in) :: kpt_in(3,nkpt_in)

!Local variables-------------------------------
!scalars
 integer :: bd2tot_index
 integer :: iband,ii,ikpt, ikpt_ddk
 integer :: isppol,l1,mband
 integer :: bantot1
 integer :: unit_ddk
 integer :: symrankkpt
 character(len=fnlen) :: filnam1,filnam2,filnam3
 character(len=500) :: msg
 type(hdr_type) :: hdr1
 type(kptrank_type) :: kptrank_t

!arrays
 real(dp) :: im_el_veloc(3)
 real(dp),allocatable :: eig1_k(:,:)
 real(dp),allocatable :: eigen11(:),eigen12(:),eigen13(:)

! *********************************************************************************

!Read data file name
!TODO: this should be standardized and read in anaddb always, not
!conditionally. Otherwise when new files are added to the anaddb files
!file...  Catastrophe!

 write(std_out,*)'enter read_el_veloc '

!Read data file
 if (open_file(elph_tr_ds%ddkfilename,msg,newunit=unit_ddk,form='formatted') /= 0) then
   MSG_ERROR(msg)
 end if

 rewind(unit_ddk)
 read(unit_ddk,'(a)')filnam1       ! first ddk file
 read(unit_ddk,'(a)')filnam2       ! second ddk file
 read(unit_ddk,'(a)')filnam3       ! third ddk file
 close (unit_ddk)

 bantot1 = 2*nband_in**2*nkpt_in*nsppol_in

 call inpgkk(eigen11,filnam1,hdr1)
 call hdr_free(hdr1)

 call inpgkk(eigen12,filnam2,hdr1)
 call hdr_free(hdr1)

!we use the hdr1 from the last call - should add some consistency
!testing here, we are trusting users not to mix different ddk files...
 call inpgkk(eigen13,filnam3,hdr1)

!Extract info from the header
 if(hdr1%nsppol /= nsppol_in) then
   MSG_ERROR('nsspol /= input nsppol')
 end if

!Get mband, as the maximum value of nband(nkpt)
 mband=maxval(hdr1%nband(1:hdr1%nkpt))
 if (mband /= nband_in) then
   MSG_ERROR('nband_in input to read_el_veloc is inconsistent with mband')
 end if

 write(std_out,*)
 write(std_out,*)                     'readings from read_el_veloc header'
 write(std_out,'(a,i8)')              ' natom                =',hdr1%natom
 write(std_out,'(a,3i8)')             ' nkpt,nband_in,mband  =',hdr1%nkpt,nband_in,mband
 write(std_out,'(a, f10.5,a)' )      ' ecut                 =',hdr1%ecut,' Ha'
 write(std_out,'(a,e15.5,a,e15.5,a)' )' fermie               =',hdr1%fermie,' Ha ',hdr1%fermie*Ha_eV,' eV'

 ABI_ALLOCATE(eig1_k,(2*nband_in**2,3))
 bd2tot_index = 0
 elph_tr_ds%el_veloc=zero

!need correspondence between the DDK kpoints and the kpt_phon
 call mkkptrank (hdr1%kptns,hdr1%nkpt,kptrank_t)

 do isppol=1,nsppol_in
   im_el_veloc(:)=zero
   do ikpt=1,nkpt_in
     call get_rank_1kpt (kpt_in(:,ikpt),symrankkpt, kptrank_t)
     ikpt_ddk = kptrank_t%invrank(symrankkpt)
     if (ikpt_ddk == -1) then
       write(std_out,*)'read_el_veloc ******** error in correspondence between ddk and gkk kpoint sets'
       write(std_out,*)' kpt sets in gkk and ddk files must agree.'
       MSG_ERROR("Aborting now")
     end if
     bd2tot_index=2*nband_in**2*(ikpt_ddk-1)

!    first derivative eigenvalues for k-point
     eig1_k(:,1)=eigen11(1+bd2tot_index:2*nband_in**2+bd2tot_index)
     eig1_k(:,2)=eigen12(1+bd2tot_index:2*nband_in**2+bd2tot_index)
     eig1_k(:,3)=eigen13(1+bd2tot_index:2*nband_in**2+bd2tot_index)

!    turn el_veloc to cartesian coordinates
     do iband=1,nband_in
       do l1=1,3
         do ii=1,3
           elph_tr_ds%el_veloc(ikpt,iband,l1,isppol)=elph_tr_ds%el_veloc(ikpt,iband,l1,isppol)+&
&           hdr1%rprimd(l1,ii)*eig1_k(2*iband-1+(iband-1)*2*nband_in,ii)/two_pi
           im_el_veloc(l1)=im_el_veloc(l1)+&
&           hdr1%rprimd(l1,ii)*eig1_k(2*iband+(iband-1)*2*nband_in,ii)/two_pi
         end do
       end do ! l1
     end do
   end do
 end do ! end isppol

 call destroy_kptrank (kptrank_t)
 ABI_DEALLOCATE(eig1_k)
 ABI_DEALLOCATE(eigen11)
 ABI_DEALLOCATE(eigen12)
 ABI_DEALLOCATE(eigen13)

 call hdr_free(hdr1)

 write(std_out,*)'out of read_el_veloc '

end subroutine read_el_veloc
!!***
