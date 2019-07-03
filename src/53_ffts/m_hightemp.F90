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
    logical :: enabled
    integer :: bcut,nbcut
    real(dp) :: ebcut,energycontrib,entropycontrib,nfreeel,u0,ucvol
  contains
    procedure :: compute_energycontrib,compute_nfreeel
    procedure :: compute_u0,init
    final :: finalize
  end type hightemp_type

  public ::freedos,hightemp_addtoenergy,hightemp_addtorho,prt_eigocc,hightemp_getnfreeel
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
  subroutine init(this,mband,nbcut,rprimd)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    integer,intent(in) :: mband,nbcut
    ! Arrays
    real(dp),intent(in) :: rprimd(3,3)

    ! Local variables -------------------------
    ! Arrays
    real(dp) :: gprimd(3,3),rmet(3,3), gmet(3,3)

    ! *********************************************************************

    this%bcut=mband
    this%nbcut=nbcut
    this%ebcut=zero
    this%energycontrib=zero
    this%entropycontrib=zero
    this%nfreeel=zero
    this%u0=zero
    call metric(gmet,gprimd,-1,rmet,rprimd,this%ucvol)
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

  !!****f* ABINIT/m_hightemp/compute_u0
  !! NAME
  !! compute_u0
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
  subroutine compute_u0(this,eigen,eknk,mband,nkpt,nsppol)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    integer,intent(in) :: mband,nkpt,nsppol
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
    this%u0=zero
    do ii=this%bcut*nkpt*nsppol-25,this%bcut*nkpt*nsppol
      this%u0=this%u0+(eigentemp(ii)-eknktemp(ii))
      niter=niter+1
    end do
    this%u0=this%u0/niter
    this%ebcut=eigentemp(this%bcut*nkpt*nsppol)
  end subroutine compute_u0

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

    call hightemp_getnfreeel(this%ebcut,this%entropycontrib,fermie,mrgrid,&
    & this%nfreeel,tsmear,this%u0,this%ucvol)
  end subroutine compute_nfreeel

  !!****f* ABINIT/m_hightemp/compute_energycontrib
  !! NAME
  !! compute_energycontrib
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
  subroutine compute_energycontrib(this,fermie,mrgrid,tsmear)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    integer,intent(in) :: mrgrid
    real(dp),intent(in) :: fermie,tsmear

    ! Local variables -------------------------
    ! Scalars
    integer :: ii
    real(dp) :: ix,step
    ! Arrays
    real(dp) :: values(mrgrid)

    ! *********************************************************************

    step=(1/this%ebcut)/mrgrid
    do ii=1,mrgrid
      ix=(ii)*step
      values(ii)=fermi_dirac(1./ix,fermie,tsmear)*freedos(1/ix,this%u0,this%ucvol)*(1/ix-this%u0)/(ix*ix)
    end do
    this%energycontrib=simpson(step,values)
  end subroutine compute_energycontrib

  ! *********************************************************************

  !!****f* ABINIT/m_hightemp/free_dos
  !! NAME
  !! free_dos
  !!
  !! FUNCTION
  !! Returns the free particle density of states for a given energy (in Hartree)
  !!
  !! INPUTS
  !! energy=get the value of the free particle density of states at this energy
  !! u0=energy shift factor
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
  subroutine prt_eigocc(eigen,etotal,energies,fnameabo_eig,iout,kptns,&
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
      & ' Kinetic energy     = ',energies%e_kinetic,' Ha         Hartree energy     = ',energies%e_hartree,&
      & ' Ha         Ewald energy       = ',energies%e_ewald, ' Ha'
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
  end subroutine prt_eigocc

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
  !! u0=energy shift factor
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
  subroutine hightemp_getnfreeel(ebcut,entropy,fermie,mrgrid,nfreeel,tsmear,u0,ucvol)

    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: mrgrid
    real(dp),intent(in) :: ebcut,fermie,tsmear,u0
    real(dp),intent(out) :: entropy,nfreeel

    ! Local variables -------------------------
    ! Scalars
    integer :: ii
    real(dp) :: ix,step,ucvol
    ! Arrays
    real(dp) :: valuesnel(mrgrid),valuesent(mrgrid)

    ! *********************************************************************

    step=(1/ebcut)/mrgrid
    do ii=1,mrgrid
      ix=(ii)*step
      valuesnel(ii)=fermi_dirac(1./ix,fermie,tsmear)*freedos(1/ix,u0,ucvol)/(ix*ix)
      valuesent(ii)=(fermi_dirac(1./ix,fermie,tsmear)*max(log(fermi_dirac(1./ix,fermie,tsmear)),-0.001_dp*huge(0.0_dp))+&
      & (1.-fermi_dirac(1./ix,fermie,tsmear))*log(1.-fermi_dirac(1./ix,fermie,tsmear)))*&
      & freedos(1/ix,u0,ucvol)/(ix*ix)
    end do

    nfreeel=simpson(step,valuesnel)
    entropy=simpson(step,valuesent)
  end subroutine hightemp_getnfreeel

  !!****f* ABINIT/m_hightemp/hightemp_addtorho
  !! NAME
  !! hightemp_addtorho
  !!
  !! FUNCTION
  !! Add the free continous part corresponding to free particle model to the existing electronic density.
  !! $n = 2 \Sum_{i=1}^{N_c} f(\epsilon_i)|\psi_i(\mathbb{r})|^2 +
  !! \Omega\int_{Ec}^{\Infty}f(\epsilon)\frac{\sqrt{2}}{\pi^2}(\epsilon - U_0)^{3/2}d \epsilon$
  !!
  !! INPUTS
  !! int_rhocontrib=computed free particle part integral which is the number of free electrons.
  !! nfft=number of FFT points.
  !! nspden=number of spin-density components
  !! rhor(cplex*nfft,nspden)=array for electron density in electrons/bohr**3.
  !! ucvol=unit cell volume (bohr^3)
  !!
  !! OUTPUT
  !! rhor(cplex*nfft,nspden)=array for electron density in electrons/bohr**3.
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine hightemp_addtorho(nfreeel,nfft,nspden,rhor,ucvol)

    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: nfft,nspden
    real(dp),intent(in) :: nfreeel,ucvol
    ! Arrays
    real(dp),intent(inout) :: rhor(nfft,nspden)

    ! *********************************************************************

    rhor(:,:)=rhor(:,:)+nfreeel/ucvol/nspden
  end subroutine hightemp_addtorho

  !!****f* ABINIT/m_hightemp/hightemp_addtorho
  !! NAME
  !! hightemp_addtorho
  !!
  !! FUNCTION
  !! Add the free continous part corresponding to free particle model to the energy.
  !!
  !! INPUTS
  !! energycontrib=computed free particle part integral corresponding to energy of free electrons
  !! etotal=energy in which we want to add the free part
  !!
  !! OUTPUT
  !! etotal=energy in which we want to add the free part
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine hightemp_addtoenergy(energycontrib,etotal)

    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: energycontrib
    real(dp),intent(inout) :: etotal

    ! *********************************************************************

    etotal=etotal+energycontrib
  end subroutine hightemp_addtoenergy

  !!****f* ABINIT/m_hightemp/hightemp_addtostress
  !! NAME
  !! hightemp_addtostress
  !!
  !! FUNCTION
  !! Add the free continous part corresponding to free particle model to stress tensor.
  !!
  !! INPUTS
  !! energycontrib=computed free particle part integral corresponding to energy of free electrons
  !!
  !! OUTPUT
  !! strten(6)=components of the stress tensor (hartree/bohr^3) for the
  !!    6 unique components of this symmetric 3x3 tensor:
  !!    Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).
  !!    The diagonal components of the returned stress tensor are
  !!    CORRECTED for the Pulay stress.
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine hightemp_addtostress(energycontrib,strten)

    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: energycontrib
    ! Arrays
    real(dp),intent(inout) :: strten(6)

    ! *********************************************************************

    strten(1:3) = strten(1:3)+-(2./3.)*energycontrib
  end subroutine hightemp_addtostress

end module m_hightemp
