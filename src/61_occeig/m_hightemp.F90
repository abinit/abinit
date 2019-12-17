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
!!  Copyright (C) 2018-2019 ABINIT group (A. Blanchet)
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
  use m_energies,       only : energies_type
  use m_mpinfo,         only : ptabs_fourdp, proc_distrb_cycle
  use m_numeric_tools,  only : simpson, simpson_int

  implicit none

  type,public :: hightemp_type
    integer :: bcut,ioptden,nbcut
    real(dp) :: ebcut,edc_kin_freeel,e_kin_freeel,e_ent_freeel
    real(dp) :: nfreeel,e_shiftfactor,ucvol
  contains
    procedure :: compute_e_kin_freeel,compute_nfreeel
    procedure :: compute_e_shiftfactor,init,destroy
  end type hightemp_type

  ! type(hightemp_type),save,pointer :: hightemp=>null()
  public :: hightemp_dosfreeel
  public :: hightemp_prt_eigocc,hightemp_getnfreeel
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
  !! mband=maximum number of bands
  !! nbcut=input value of on how many states would you like to average the band energy
  !! rprimd(3,3)=dimensional primitive translations in real space (bohr)
  !!
  !! OUTPUT
  !! this=hightemp_type object concerned
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine init(this,ioptden,mband,nbcut,rprimd)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    integer,intent(in) :: ioptden,mband,nbcut
    ! Arrays
    real(dp),intent(in) :: rprimd(3,3)

    ! Local variables -------------------------
    ! Arrays
    real(dp) :: gprimd(3,3),rmet(3,3), gmet(3,3)

    ! *********************************************************************

    this%ioptden=ioptden
    this%bcut=mband
    this%nbcut=nbcut
    this%ebcut=zero
    this%edc_kin_freeel=zero
    this%e_kin_freeel=zero
    this%e_ent_freeel=zero
    this%nfreeel=zero
    this%e_shiftfactor=zero
    call metric(gmet,gprimd,-1,rmet,rprimd,this%ucvol)
  end subroutine init

  subroutine destroy(this)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this

    ! *********************************************************************

    this%ioptden=0
    this%bcut=0
    this%nbcut=0
    this%ebcut=zero
    this%edc_kin_freeel=zero
    this%e_kin_freeel=zero
    this%e_ent_freeel=zero
    this%nfreeel=zero
    this%e_shiftfactor=zero
    this%ucvol=zero
  end subroutine destroy

  !!****f* ABINIT/m_hightemp/compute_e_shiftfactor
  !! NAME
  !! compute_e_shiftfactor
  !!
  !! FUNCTION
  !! Compute the shift factor $U_0$ that appears in the density of states of free electrons.
  !!
  !! INPUTS
  !! this=hightemp_type object concerned
  !! eigen(mband*nkpt*nsppol)=eigenvalues (hartree)
  !! eknk(mband*nkpt*nsppol)=kinetic energies (hartree)
  !! mband=maximum number of bands
  !! nkpt=number of k points
  !! nsppol=1 for unpolarized, 2 for spin-polarized
  !!
  !! OUTPUT
  !! this=hightemp_type object concerned
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine compute_e_shiftfactor(this,eigen,eknk,mband,mpi_enreg,nkpt,nsppol)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    integer,intent(in) :: mband,nkpt,nsppol
    type(MPI_type), intent(inout) :: mpi_enreg
    ! Arrays
    real(dp),intent(in) :: eigen(mband*nkpt*nsppol)
    real(dp),intent(in) :: eknk(mband*nkpt*nsppol)

    ! Local variables -------------------------
    ! Scalars
    integer :: ii,krow,niter
    ! Arrays
    real(dp) :: eigentemp(mband*nkpt*nsppol)
    real(dp) :: eknktemp(mband*nkpt*nsppol)
    logical :: mk(mband*nkpt*nsppol)

    ! *********************************************************************

    eigentemp(:)=zero
    eknktemp(:)=zero
    mk(:)=.true.

    ! Sorting eigen and eknk in ascending energy order.
    do ii=1,this%bcut*nkpt*nsppol
      krow=minloc(eigen,dim=1,mask=mk)
      mk(minloc(eigen,dim=1,mask=mk))=.false.
      eigentemp(ii)=eigen(krow)
      eknktemp(ii)=eknk(krow)
    end do

    ! Doing the average over the 25 lasts states...
    niter=0
    this%e_shiftfactor=zero
    do ii=this%bcut*nkpt*nsppol-25,this%bcut*nkpt*nsppol
      this%e_shiftfactor=this%e_shiftfactor+(eigentemp(ii)-eknktemp(ii))
      niter=niter+1
    end do
    this%e_shiftfactor=this%e_shiftfactor/niter
    this%ebcut=eigentemp(this%bcut*nkpt*nsppol)
  end subroutine compute_e_shiftfactor

  !!****f* ABINIT/m_hightemp/compute_nfreeel
  !! NAME
  !! compute_nfreeel
  !!
  !! FUNCTION
  !! Compute the value of the integral corresponding to the missing free electrons after energy cut.
  !! $I = \int_{Ec}^{\Infty}f(\epsilon)\frac{\sqrt{2}\Omega}{\pi^2}\sqrt{\epsilon - U_0}d \epsilon$
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
  subroutine compute_nfreeel(this,fermie,mrgrid,tsmear)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    integer,intent(in) :: mrgrid
    real(dp),intent(in) :: fermie,tsmear

    ! *********************************************************************

    call hightemp_getnfreeel(this%ebcut,this%e_ent_freeel,fermie,mrgrid,&
    & this%nfreeel,tsmear,this%e_shiftfactor,this%ucvol)
  end subroutine compute_nfreeel

  !!****f* ABINIT/m_hightemp/compute_e_kin_freeel
  !! NAME
  !! compute_e_kin_freeel
  !!
  !! FUNCTION
  !! Compute the value of the integral corresponding to the residual of energy after the band cut
  !! $I = \Omega\int_{Ec}^{\Infty}f(\epsilon)\frac{\sqrt{2}}{\pi^2}(\epsilon - U_0)^{3/2}d \epsilon$
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
  subroutine compute_e_kin_freeel(this,fermie,mrgrid,nfftf,nspden,tsmear,vtrial)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    integer,intent(in) :: mrgrid,nfftf,nspden
    real(dp),intent(in) :: fermie,tsmear
    ! Arrays
    real(dp),intent(in) :: vtrial(nfftf,nspden)

    ! Local variables -------------------------
    ! Scalars
    integer :: ii,ifft,ispden
    real(dp) :: ix,step
    ! Arrays
    real(dp),dimension(:),allocatable :: values

    ! *********************************************************************
    step=1e-5

    ! Dynamic array find size
    ix=this%ebcut
    ii=0
    do while(fermi_dirac(ix,fermie,tsmear)>tol16)
      ii=ii+1
      ix=ix+step
    end do
    ABI_ALLOCATE(values,(ii))

    ! Computation of e_kin_freeel
    ix=this%ebcut
    ii=0
    do while(fermi_dirac(ix,fermie,tsmear)>tol16)
      ii=ii+1
      values(ii)=fermi_dirac(ix,fermie,tsmear)*&
      & hightemp_dosfreeel(ix,this%e_shiftfactor,this%ucvol)*(ix-this%e_shiftfactor)
      ix=ix+step
    end do
    if (ii>1) then
      this%e_kin_freeel=simpson(step,values)
    end if

    ! Computation of edc_kin_freeel
    this%edc_kin_freeel=zero
    do ispden=1,nspden
      do ifft=1,nfftf
        this%edc_kin_freeel=this%edc_kin_freeel+vtrial(ifft,ispden)
      end do
    end do
    ! Verifier la constante (/nspden**2)
    this%edc_kin_freeel=this%edc_kin_freeel*this%nfreeel/nspden/nfftf/nspden
  end subroutine compute_e_kin_freeel

  ! *********************************************************************

  !!****f* ABINIT/m_hightemp/hightemp_free_dos
  !! NAME
  !! hightemp_free_dos
  !!
  !! FUNCTION
  !! Returns the free particle density of states for a given energy (in Hartree)
  !!
  !! INPUTS
  !! energy=get the value of the free particle density of states at this energy
  !! e_shiftfactor=energy shift factor
  !! ucvol=unit cell volume (bohr^3)
  !!
  !! OUTPUT
  !! freedos=value of free particle density of states at given energy
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function hightemp_dosfreeel(energy,e_shiftfactor,ucvol)

    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: energy,e_shiftfactor,ucvol
    real(dp) :: hightemp_dosfreeel

    ! *********************************************************************

    hightemp_dosfreeel=sqrt(2.)*ucvol*sqrt(energy-e_shiftfactor)/(PI*PI)
  end function hightemp_dosfreeel

  !!****f* ABINIT/m_hightemp/hightemp_prt_eigocc
  !! NAME
  !! hightemp_prt_eigocc
  !!
  !! FUNCTION
  !! Printing in a _EIGOCC file many informations like Fermi energy, Unit cell volume,
  !! electronic temperature (tsmear when occopt=3), eigenvalues and occupation for each k-point,
  !! wheight of the k-point, number of bands, number of k-points and more...
  !! This file is intended to be used for custom DOS computation with external tools for example.
  !!
  !! INPUTS
  !! eigen(mband*nkpt*nsppol)=eigenvalues (hartree)
  !! etotal=total energy
  !! energies=object which contains all energies of the computation
  !! fnameabo_eig=filename for printing of the eigenenergies
  !! iout=unit number for formatted output file
  !! kptns(3,nkpt)=k points in reduced coordinates
  !! mband=maximum number of bands
  !! nband(nkpt)=number of bands at each k point
  !! nkpt=number of k points
  !! nsppol=1 for unpolarized, 2 for spin-polarized
  !! occ(maxval(nband(:))*nkpt*nsppol)=occupancies for each band and k point
  !! rprimd(3,3)=lattice vectors in Bohr
  !! strten(6)=components of the stress tensor (hartree/bohr^3)
  !! tsmear=smearing width (or temperature)
  !! usepaw= 0 for non paw calculation; =1 for paw calculation
  !! wtk(nkpt)=k-point weights
  !!
  !! OUTPUT
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine hightemp_prt_eigocc(eigen,etotal,energies,fnameabo_eig,iout,kptns,&
    & mband,nband,nkpt,nsppol,occ,rprimd,strten,tsmear,usepaw,wtk)

    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: iout,mband,nkpt,nsppol,usepaw
    real(dp),intent(in) :: etotal,tsmear
    character(len=*),intent(in) :: fnameabo_eig
    type(energies_type),intent(in) :: energies
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: eigen(mband*nkpt*nsppol),kptns(3,nkpt)
    real(dp),intent(in) :: rprimd(3,3),strten(6)
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

    ! write(msg, '(a,f12.6,a,f12.6,a,ES16.8,a)') &
    !   & ' Fermi (or HOMO) energy (hartree)= ',fermie,'   Average Vxc (hartree)=',vxcavg,&
    !   & ', Total energy= ',etotal, ' Ha'
    ! call wrtout(temp_unit,msg,'COLL')
    ! write(msg, '(a,e16.8,a,f16.8,a,a,ES12.4,a)') &
    !   & ' Unit cell volume ucvol= ',ucvol,' bohr^3, Electronic temperature= ',tsmear,' Ha',&
    !   & ', Pressure= ',-(strten(1)+strten(2)+strten(3))*HaBohr3_GPa/3.0_dp,' GPa'
    ! call wrtout(temp_unit,msg,'COLL')

    write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
      & ' Chemical potential = ',energies%e_fermie,' Ha         XC energy          = ',energies%e_xc,&
      & ' Ha         Total free energy  = ',etotal, ' Ha'
    call wrtout(temp_unit,msg,'COLL')

    write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
    & ' Kinetic energy     = ',energies%e_kinetic,' Ha         Kin. free el. E    = ',energies%e_kin_freeel,&
    & ' Ha         Total kin. energy  = ',energies%e_kinetic+energies%e_kin_freeel, ' Ha'
    call wrtout(temp_unit,msg,'COLL')

    write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
    & ' Hartree energy     = ',energies%e_hartree,' Ha         Ewald energy       = ',energies%e_ewald,&
    & ' Ha         E. Shift factor U0 = ',energies%e_shiftfactor, ' Ha'
    call wrtout(temp_unit,msg,'COLL')

    if (usepaw==0) then
      write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
        & ' PsPCore energy     = ',energies%e_corepsp,' Ha         Loc. Psp. energy   = ',energies%e_localpsp,&
        & ' Ha         NL PsP energy      = ',energies%e_nlpsp_vfock, ' Ha'
    else if (usepaw==1) then
      write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
        & ' PsPCore energy     = ',energies%e_corepsp,' Ha         Loc. Psp. energy   = ',energies%e_localpsp,&
        & ' Ha         PAW PsP energy     = ',energies%e_paw, ' Ha'
    end if
    call wrtout(temp_unit,msg,'COLL')

    write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
      & ' Entropy energy     = ',energies%e_entropy,' Ha         Band energy        = ',energies%e_eigenvalues,&
      & ' Ha         Internal energy    = ',etotal-energies%e_entropy, ' Ha'
    call wrtout(temp_unit,msg,'COLL')
    write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
      & ' Unit cell vol      = ',ucvol,' Bohr^3     Elec. temperature  = ',tsmear,&
      & ' Ha         Pressure           = ',-(strten(1)+strten(2)+strten(3))*HaBohr3_GPa/3.0_dp,' GPa'
    call wrtout(temp_unit,msg,'COLL')

    ! Loop over spins
    do isppol=1,nsppol
      ! write(msg, '(a,i6,a)') ' Eigenvalues (hartree) for nkpt=',nkpt,'k points:'
      write(msg, '(a,i12,a,i12)') ' Number of kpts     = ',nkpt,'            Number of bands    = ',&
      & mband
      call wrtout(temp_unit,msg,'COLL')

      ! Loop over k-points
      do ikpt=1,nkpt
        nband_k=nband(ikpt+(isppol-1)*nkpt)
        write(msg, '(a,i6,a,i6,a,ES12.5,a,3f8.4,a)') &
          & ' ------- ikpt = ',ikpt,'      nband = ',nband_k,'     ikpt weight = ',wtk(ikpt)+tol10,&
          & '   Reduced coordinates : ',kptns(1:3,ikpt)+tol10,' -------'
        call wrtout(temp_unit,msg,'COLL')

        ! Loop over bands
        do ii=0,(nband_k-1)/4
          write(msg, '(4(i6,ES12.4,ES13.5,a,1x))') &
            & (iband,eigen(iband+band_index),occ(iband+band_index),',',iband=1+ii*4,min(nband_k,4+ii*4))
          call wrtout(temp_unit,msg,'COLL')
        end do

        band_index=band_index+nband_k
      end do ! do ikpt=1,nkpt
    end do ! do isppol=1,nsppol

    close(temp_unit)
  end subroutine hightemp_prt_eigocc

  !!****f* ABINIT/m_hightemp/hightemp_getnfreeel
  !! NAME
  !! hightemp_getnfreeel
  !!
  !! FUNCTION
  !! Compute the value of the integral corresponding to the missing free electrons after energy cut.
  !! $I = \int_{Ec}^{\Infty}f(\epsilon)\frac{\sqrt{2}\Omega}{\pi^2}\sqrt{\epsilon - U_0}d \epsilon$
  !!
  !! INPUTS
  !! this=hightemp_type object concerned
  !! fermie=fermi energy (Hartree)
  !! mrgrid=number of grid points to compute the integral
  !! tsmear=smearing width (or temperature)
  !! e_shiftfactor=energy shift factor
  !! ucvol=unit cell volume (bohr^3)
  !!
  !! OUTPUT
  !! nfreeel=number of free electrons after the energy cut with given fermi level
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine hightemp_getnfreeel(ebcut,entropy,fermie,mrgrid,nfreeel,tsmear,e_shiftfactor,ucvol)

    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: mrgrid
    real(dp),intent(in) :: ebcut,fermie,tsmear,e_shiftfactor
    real(dp),intent(out) :: entropy,nfreeel

    ! Local variables -------------------------
    ! Scalars
    integer :: ii
    real(dp) :: ix,step,ucvol
    ! Arrays
    real(dp),dimension(:),allocatable :: valuesnel,valuesent

    ! *********************************************************************
    step=1e-4
    nfreeel=zero
    entropy=zero

    ! Dynamic array find size
    ix=ebcut
    ii=0
    do while(fermi_dirac(ix,fermie,tsmear)>1e-8)
      ii=ii+1
      ix=ix+step
    end do

    ABI_ALLOCATE(valuesnel,(ii))
    ABI_ALLOCATE(valuesent,(ii))

    ix=ebcut
    ii=0
    do while(fermi_dirac(ix,fermie,tsmear)>1e-8)
      ii=ii+1
      valuesnel(ii)=fermi_dirac(ix,fermie,tsmear)*hightemp_dosfreeel(ix,e_shiftfactor,ucvol)
      valuesent(ii)=(fermi_dirac(ix,fermie,tsmear)*log(fermi_dirac(ix,fermie,tsmear))+&
      & (1.-fermi_dirac(ix,fermie,tsmear))*log(1.-fermi_dirac(ix,fermie,tsmear)))*&
      & hightemp_dosfreeel(ix,e_shiftfactor,ucvol)
      ix=ix+step
    end do
    if (ii>1) then
      nfreeel=simpson(step,valuesnel)
      entropy=simpson(step,valuesent)
    end if
  end subroutine hightemp_getnfreeel

end module m_hightemp
