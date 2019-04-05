!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_hightemp
!! NAME
!!  m_hightemp
!!
!! FUNCTION
!!  This module provides routines to run computations at very high temperatures
!!  with orbitals.
!!
!! COPYRIGHT
!!  Copyright (C) 2018-2019 ABINIT group (??)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
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

module m_hightemp
  use defs_basis
  use defs_abitypes
  use m_io_tools
  use m_errors
  use m_geometry
  use m_specialmsg
  use m_mpinfo,       only : ptabs_fourdp, proc_distrb_cycle

  implicit none

  public :: mkrho_hightemp
  public :: prt_eigocc
  private :: free_transfactor
contains

!!****f* ABINIT/prt_eigocc
!! NAME
!! prt_eigocc
!!
!! FUNCTION
!! Printing in a _EIGOCC file many informations like Fermi energy, Average Vxc, Unit cell volume,
!! electronic temperature (tsmear when occopt=3), eigenvalues and occupation for each k-point,
!! wheight of the k-point, number of bands, number of k-points...
!! This file is intended to be used for custom DOS computation with external tools for example.
!!
!! INPUTS
!! eigen(mband*nkpt*nsppol)=eigenvalues (hartree)
!! fermie=fermi energy (Hartree)
!! fnameabo_eig=filename for printing of the eigenenergies
!! iout=unit number for formatted output file
!! kptns(3,nkpt)=k points in reduced coordinates
!! mband=maximum number of bands
!! nband(nkpt)=number of bands at each k point
!! nkpt=number of k points
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! occ(maxval(nband(:))*nkpt*nsppol)=occupancies for each band and k point
!! rprimd(3,3)=Lattive vectors in Bohr
!! tsmear=smearing width (or temperature)
!! vxcavg=average of vxc potential
!! wtk(nkpt)=k-point weights
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine prt_eigocc(eigen,fermie,fnameabo_eig,iout,kptns,&
  & mband,nband,nkpt,nsppol,occ,rprimd,tsmear,vxcavg,wtk)

  ! Arguments -------------------------------
  ! Scalars
  integer,intent(in) :: iout,mband,nkpt,nsppol
  real(dp),intent(in) :: fermie,tsmear,vxcavg
  character(len=*),intent(in) :: fnameabo_eig
  ! Arrays
  integer,intent(in) :: nband(nkpt*nsppol)
  real(dp),intent(in) :: eigen(mband*nkpt*nsppol),kptns(3,nkpt)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: occ(mband*nkpt*nsppol)
  real(dp),intent(in) :: wtk(nkpt)

  ! Local variables -------------------------
  ! Scalars
  integer :: band_index,iband,ii,ikpt,isppol,nband_k,temp_unit
  real(dp) :: ucvol
  character(len=200) :: fnameabo_eigocc
  character(len=39) :: kind_of_output
  character(len=500) :: msg
  ! Arrays
  real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)

  ! *********************************************************************

  band_index=0
  call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

  fnameabo_eigocc = trim(fnameabo_eig) // trim('OCC')
  write(msg,'(a,a)') ' prt_eigocc : about to open file ',trim(fnameabo_eigocc)
  call wrtout(iout,msg,'COLL')
  if (open_file(fnameabo_eigocc,msg,newunit=temp_unit,status='unknown',form='formatted') /= 0) then
    MSG_ERROR(msg)
  end if
  rewind(temp_unit)

  write(msg, '(a,f12.6,a,f12.6)') &
    & ' Fermi (or HOMO) energy (hartree) =',fermie,'   Average Vxc (hartree)=',vxcavg
  call wrtout(temp_unit,msg,'COLL')
  write(msg, '(a,e16.8,a,f16.8,a)') &
    & ' Unit cell volume ucvol= ',ucvol,' bohr^3, electronic temperature= ',tsmear,' Ha'
  call wrtout(temp_unit,msg,'COLL')

  ! Loop over spins
  do isppol=1,nsppol
    write(msg, '(a,i6,a)') ' Eigenvalues (hartree) for nkpt=',nkpt,'k points:'
    call wrtout(temp_unit,msg,'COLL')

    ! Loop over k-points
    do ikpt=1,nkpt
      nband_k=nband(ikpt+(isppol-1)*nkpt)
      write(msg, '(a,i6,a,i6,a,f9.5,a,3f8.4,a)') &
        & ' kpt#',ikpt,', nband=',nband_k,', wtk=',wtk(ikpt)+tol10,', kpt=', &
        & kptns(1:3,ikpt)+tol10,' (reduced coord)'
      call wrtout(temp_unit,msg,'COLL')

      ! Loop over bands
      do ii=0,(nband_k-1)/4
        write(msg, '(4(i6,f13.8,f12.8,a,1x))') &
          & (iband,eigen(iband+band_index),occ(iband+band_index),',',iband=1+ii*4,min(nband_k,4+ii*4))
        call wrtout(temp_unit,msg,'COLL')
      end do

      band_index=band_index+nband_k
    end do ! do ikpt=1,nkpt
  end do ! do isppol=1,nsppol

  close(temp_unit)
end subroutine prt_eigocc

!!****f* ABINIT/free_transfactor
!! NAME
!! free_transfactor
!!
!! FUNCTION
!! Compute the translation factor $U_0$ that appears in the density of states of free electrons.
!!
!! TODO
!! INPUTS
!!
!! TODO
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine free_transfactor()
end subroutine free_transfactor

!!****f* ABINIT/mkrho_hightemp
!! NAME
!! mkrho_hightemp
!!
!! FUNCTION
!! Make the electronic density replacing high energy bands by free electron states.
!!
!! TODO
!! INPUTS
!!
!! TODO
!! OUTPUT
!!
!! NOTES
!! TODO Parallelization is still to be done.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
subroutine mkrho_hightemp(dtset,mpi_enreg,npwarr,occ,rhor)

  ! Arguments -------------------------------
  ! Scalars
  type(MPI_type),intent(in) :: mpi_enreg
  type(dataset_type),intent(in) :: dtset
  ! Arrays
  integer,intent(in) :: npwarr(dtset%nkpt)
  real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(out) :: rhor(dtset%nfft,dtset%nspden)

  ! Local variables -------------------------
  ! Scalars
  integer :: band_cut,bdtot_index,iband,ikpt,isppol,me,n1,n2,n3,n4,n5,n6,nband_k,npw_k
  ! Arrays
  real(dp),allocatable :: rhoaug(:,:,:)

  ! *********************************************************************

  band_cut = 100
  me=mpi_enreg%me_kpt
  bdtot_index=0
  ! n4,n5,n6 are FFT dimensions, modified to avoid cache trashing
  n1=dtset%ngfft(1);n2=dtset%ngfft(2);n3=dtset%ngfft(3)
  n4=dtset%ngfft(4);n5=dtset%ngfft(5);n6=dtset%ngfft(6)
  ABI_ALLOCATE(rhoaug,(n4,n5,n6))

  ! Loop over spins
  do isppol=1,dtset%nsppol
    rhoaug(:,:,:)=zero

    ! Loop over k-points
    do ikpt=1,dtset%nkpt
      nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
      npw_k=npwarr(ikpt)

      if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
        bdtot_index=bdtot_index+nband_k
        cycle
      end if
      ! Loop over bands
      do iband=1,nband_k

        ! Treating only occupied bands
        if(abs(occ(iband+bdtot_index))>tol8) then

          if(iband<=band_cut) then
            ! call fourwf(1,rhoaug(:,:,:),cwavef(:,:,1),dummy,wfraug,gbound,gbound,&
            !   & istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,&
            !   & npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight,weight_i,&
            !   & use_ndo=use_nondiag_occup_dmft,fofginb=cwavefb(:,:,1),&
            !   & use_gpu_cuda=dtset%use_gpu_cuda)
          end if
        end if
      end do ! End loop over bands
    end do ! End loop over k-points
  end do ! End loop over spins
end subroutine mkrho_hightemp

end module m_hightemp
