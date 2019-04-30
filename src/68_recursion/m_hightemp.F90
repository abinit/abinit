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
  use m_special_funcs
  use m_specialmsg
  use m_mpinfo,         only : ptabs_fourdp, proc_distrb_cycle
  use m_numeric_tools,  only : simpson, simpson_int

  implicit none

  type,public :: hightemp_type
    logical :: enabled
    integer :: bcut,nbcut
    real(dp) :: ebcut,int_energycontrib,int_rhocontrib,u0,ucvol
  contains
    procedure :: compute_int_energycontrib,compute_int_rhocontrib
    procedure :: compute_obj,compute_u0,init
    final :: finalize
  end type hightemp_type

  public :: freedos,hightemp_addtorho,prt_eigocc
contains

  !!****f* ABINIT/m_hightemp/init
  !! NAME
  !! init
  !!
  !! FUNCTION
  !! Initialize hightemp_type object, memory allocation of arrays...
  !!
  !! INPUTS
  !! this=hightemp_type object concerned
  !! activate=condition to activate or not the hightemp methods
  !!
  !! OUTPUT
  !! this=hightemp_type object concerned
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine init(this,activate,mband,nbcut,rprimd)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    logical,intent(in) :: activate
    integer,intent(in) :: mband,nbcut
    ! Arrays
    real(dp),intent(in) :: rprimd(3,3)

    ! Local variables -------------------------
    ! Arrays
    real(dp) :: gprimd(3,3),rmet(3,3), gmet(3,3)

    ! *********************************************************************

    this%enabled = activate
    if(this%enabled) then
      call metric(gmet,gprimd,-1,rmet,rprimd,this%ucvol)
      this%bcut=mband
      this%nbcut=nbcut
      this%ebcut=zero
      this%int_energycontrib=zero
      this%int_rhocontrib=zero
      this%u0=zero
    end if
  end subroutine init

  !!****f* ABINIT/m_hightemp/finalize
  !! NAME
  !! finalize
  !!
  !! FUNCTION
  !! Finalize hightemp_type object, deallocation of arrays...
  !!
  !! INPUTS
  !! this=hightemp_type object concerned
  !!
  !! OUTPUT
  !! this=hightemp_type object concerned
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine finalize(this)

    ! Arguments -------------------------------
    ! Scalars
    type(hightemp_type),intent(inout) :: this

    ! DEALLOCATE THINGS
  end subroutine finalize

  !!****f* ABINIT/m_hightemp/compute_obj
  !! NAME
  !! compute_obj
  !!
  !! FUNCTION
  !! Compute differents needed quantities for the object hightemp_type
  !!
  !! INPUTS
  !! this=hightemp_type object concerned
  !! eigen(mband*nkpt*nsppol)=eigenvalues (hartree)
  !! eknk(mband*nkpt*nsppol)=kinetic energies (hartree)
  !! fermie=fermi energy (Hartree)
  !! mband=maximum number of bands
  !! nband(nkpt)=number of bands at each k point
  !! nkpt=number of k points
  !! nsppol=1 for unpolarized, 2 for spin-polarized
  !! tsmear=smearing width (or temperature)
  !! wtk(nkpt)=k-point weights
  !!
  !! OUTPUT
  !! this=hightemp_type object concerned
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine compute_obj(this,eigen,eknk,fermie,mband,&
    & nband,nkpt,nsppol,tsmear,wtk)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    integer,intent(in) :: mband,nkpt,nsppol
    real(dp),intent(in) :: fermie,tsmear
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
    real(dp),intent(in) :: eknk(mband*nkpt*nsppol)
    real(dp),intent(in) :: wtk(nkpt)

    ! *********************************************************************

    call this%compute_u0(eigen,eknk,mband,nband,nkpt,nsppol,wtk)
    call this%compute_int_rhocontrib(fermie,1024,tsmear)
    call this%compute_int_energycontrib(fermie,1024,tsmear)

    write(0,*) this%u0, this%int_rhocontrib
  end subroutine compute_obj

  !!****f* ABINIT/m_hightemp/compute_u0
  !! NAME
  !! compute_u0
  !!
  !! FUNCTION
  !! Compute the translation factor $U_0$ that appears in the density of states of free electrons.
  !!
  !! INPUTS
  !! this=hightemp_type object concerned
  !! eigen(mband*nkpt*nsppol)=eigenvalues (hartree)
  !! eknk(mband*nkpt*nsppol)=kinetic energies (hartree)
  !! mband=maximum number of bands
  !! nband(nkpt)=number of bands at each k point
  !! nkpt=number of k points
  !! nsppol=1 for unpolarized, 2 for spin-polarized
  !! wtk(nkpt)=k-point weights
  !!
  !! OUTPUT
  !! this=hightemp_type object concerned
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine compute_u0(this,eigen,eknk,mband,nband,nkpt,nsppol,wtk)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    integer,intent(in) :: mband,nkpt,nsppol
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
    real(dp),intent(in) :: eknk(mband*nkpt*nsppol)
    real(dp),intent(in) :: wtk(nkpt)

    ! Local variables -------------------------
    ! Scalars
    integer :: bdtot_index,iband,ikpt,isppol,nband_k,niter
    ! Arrays
    real(dp) :: eig_n(mband),ek_n(mband)

    ! *********************************************************************

    eig_n(:)=zero
    ek_n(:)=zero
    bdtot_index=1
    do isppol=1,nsppol
      do ikpt=1,nkpt
        nband_k=nband(ikpt+(isppol-1)*nkpt)
        do iband=1,nband_k
          eig_n(iband)=eig_n(iband)+wtk(ikpt)*eigen(bdtot_index)
          ek_n(iband)=ek_n(iband)+wtk(ikpt)*eknk(bdtot_index)
          bdtot_index=bdtot_index+1
        end do
      end do
    end do

    this%u0=zero
    niter=0
    do iband=this%bcut-this%bcut/10,this%bcut
      this%u0=this%u0+(eig_n(iband)-ek_n(iband))
      niter=niter+1
    end do
    this%u0=this%u0/niter
    this%ebcut=eig_n(this%bcut)
  end subroutine compute_u0

  !!****f* ABINIT/m_hightemp/compute_int_rhocontrib
  !! NAME
  !! compute_int_rhocontrib
  !!
  !! FUNCTION
  !! Compute the value of the integral corresponding to the residual of density after the band cut
  !! I = \int_{Ec}^{\Infty}f(\epsilon)\frac{\sqrt{2}}{\pi^2}\sqrt{\epsilon - U_0}d \epsilon
  !!
  !! INPUTS
  !! this=hightemp_type object concerned
  !! fermie=fermi energy (Hartree)
  !! mrgrid=number of grid points to compute the integral
  !! tsmear=smearing width (or temperature)
  !!
  !! OUTPUT
  !! this=hightemp_type object concerned
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine compute_int_rhocontrib(this,fermie,mrgrid,tsmear)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    integer,intent(in) :: mrgrid
    real(dp),intent(in) :: fermie,tsmear

    ! Local variables -------------------------
    ! Scalars
    integer :: ii
    real(dp) :: ix,step,ucvol
    ! Arrays
    real(dp) :: values(mrgrid)

    ! *********************************************************************

    step=(1/this%ebcut)/mrgrid
    do ii=1,mrgrid
      ix=(ii)*step
      values(ii)=fermi_dirac(1./ix,fermie,tsmear)*freedos(1/ix,this%u0,this%ucvol)/(ix*ix)/this%ucvol
    end do
    this%int_rhocontrib=simpson(step,values)
  end subroutine compute_int_rhocontrib

  !!****f* ABINIT/m_hightemp/compute_int_energycontrib
  !! NAME
  !! compute_int_energycontrib
  !!
  !! FUNCTION
  !! Compute the value of the integral corresponding to the residual of energy after the band cut
  !! I = \Omega\int_{Ec}^{\Infty}f(\epsilon)\frac{\sqrt{2}}{\pi^2}(\epsilon - U_0)^{3/2}d \epsilon
  !!
  !! INPUTS
  !! this=hightemp_type object concerned
  !! fermie=fermi energy (Hartree)
  !! mrgrid=number of grid points to compute the integral
  !! tsmear=smearing width (or temperature)
  !!
  !! OUTPUT
  !! this=hightemp_type object concerned
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine compute_int_energycontrib(this,fermie,mrgrid,tsmear)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    integer,intent(in) :: mrgrid
    real(dp),intent(in) :: fermie,tsmear

    ! Local variables -------------------------
    ! Scalars
    integer :: ii
    real(dp) :: ix,step,ucvol
    ! Arrays
    real(dp) :: values(mrgrid)

    ! *********************************************************************

    step=(1/this%ebcut)/mrgrid
    do ii=1,mrgrid
      ix=(ii)*step
      values(ii)=fermi_dirac(1./ix,fermie,tsmear)*freedos(1/ix,this%u0,this%ucvol)*(1/ix-this%u0)/(ix*ix)
    end do
    this%int_energycontrib=simpson(step,values)

    write(0,*) this%int_energycontrib
  end subroutine compute_int_energycontrib

  ! *********************************************************************

  !!****f* ABINIT/m_hightemp/free_dos
  !! NAME
  !! free_dos
  !!
  !! FUNCTION
  !! Returns the free density of states for a given energy (in Hartree)
  !!
  !! INPUTS
  !! this=hightemp_type object concerned
  !!
  !! OUTPUT
  !! this=hightemp_type object concerned
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function freedos(energy,u0,ucvol)

    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: energy,u0,ucvol
    real(dp) :: freedos

    ! *********************************************************************

    freedos=sqrt(2.)*ucvol*sqrt(energy-u0)/(PI*PI)
  end function freedos

  !!****f* ABINIT/m_hightemp/prt_eigocc
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
  !! rprimd(3,3)=Lattice vectors in Bohr
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

  !!****f* ABINIT/m_hightemp/hightemp_addtorho
  !! NAME
  !! hightemp_addtorho
  !!
  !! FUNCTION
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
  subroutine hightemp_addtorho(int_rhocontrib,nfft,nspden,rhor)

    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: nfft,nspden
    real(dp),intent(in) :: int_rhocontrib
    ! Arrays
    real(dp),intent(inout) :: rhor(nfft,nspden)

    ! *********************************************************************

    if(nspden==1) then
      rhor(:,:)=rhor(:,:)+int_rhocontrib
    else if(nspden==2) then
      rhor(:,:)=rhor(:,:)+.5*int_rhocontrib
    end if

    ! SHOULD WE ADD A CONTRIBUTION TO RHOG ?
  end subroutine hightemp_addtorho

  !!****f* ABINIT/m_hightemp/hightemp_addtorho
  !! NAME
  !! hightemp_addtorho
  !!
  !! FUNCTION
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
  subroutine hightemp_addtoenergy(int_energycontrib,etotal)

    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: int_energycontrib
    real(dp),intent(inout) :: etotal

    ! *********************************************************************

    etotal=etotal+int_energycontrib
  end subroutine hightemp_addtoenergy

end module m_hightemp
