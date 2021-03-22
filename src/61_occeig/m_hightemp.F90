!!****m* ABINIT/m_hightemp
!! NAME
!! m_hightemp
!!
!! FUNCTION
!! This module provides routines to run computations at very high temperatures
!! with limited number of bands. High energy orbitals are represented as
!! pure planewaves.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2020 ABINIT group (A. Blanchet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! TODO
!! 1) Add contribution to conductivity.
!! 2) Smooth the contributions.
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
  use m_xmpi
  use m_energies,       only : energies_type
  use m_gsphere,        only : getkpgnorm
  use m_kg,             only : mkkin
  use m_mpinfo,         only : ptabs_fourdp, proc_distrb_cycle
  use m_numeric_tools,  only : simpson, simpson_int
  use m_spacepar,       only : meanvalue_g

  implicit none
  public :: dip12,djp12,dip32,djp32,hightemp_e_heg
  public :: hightemp_dosfreeel
  public :: hightemp_prt_eigocc
  public :: prt_cg_from_wf
  !!***

  !----------------------------------------------------------------------

  !!****t* m_hightemp/hightemp_type
  !! NAME
  !! hightemp_type
  !!
  !! FUNCTION
  !! Store hightemp functions and parameters.
  !!
  !! SOURCE
  type,public :: hightemp_type
    integer :: bcut,nbcut,nfftf,nspden,version
    real(dp) :: ebcut,edc_kin_freeel,e_kin_freeel,ent_freeel
    real(dp) :: nfreeel,e_shiftfactor,ucvol
    real(dp),allocatable :: vtrial(:,:)
    logical :: prt_cg
  contains
    procedure :: compute_efreeel
    procedure :: compute_ent_freeel
    procedure :: compute_nfreeel
    procedure :: compute_pw_avg_std
    procedure :: compute_e_shiftfactor
    procedure :: init
    procedure :: destroy
  end type hightemp_type
  !!***

contains

  !!****f* ABINIT/m_hightemp/init
  !! NAME
  !!  init
  !!
  !! FUNCTION
  !!  Initialize hightemp_type object, memory allocation of arrays...
  !!
  !! INPUTS
  !!  this=hightemp_type object concerned
  !!  mband=maximum number of bands
  !!  nbcut=number of states used to average the constant potential value
  !!  prt_cg=debug input to print wavefunctions planewaves coefficients.
  !!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
  !!  version=hightemp implementation version
  !!
  !! OUTPUT
  !!  this=hightemp_type object concerned
  !!
  !! PARENTS
  !!  m_gstate
  !!
  !! CHILDREN
  !!  metric
  !!
  !! SOURCE
  subroutine init(this,mband,nbcut,nfftf,nspden,prt_cg,rprimd,version)
    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    integer,intent(in) :: mband,nbcut,nfftf,nspden,version,prt_cg
    ! Arrays
    real(dp),intent(in) :: rprimd(3,3)

    ! Local variables -------------------------
    ! Arrays
    real(dp) :: gprimd(3,3),rmet(3,3), gmet(3,3)

    ! *********************************************************************

    this%bcut=mband
    this%nbcut=nbcut
    this%version=version
    ABI_MALLOC(this%vtrial,(nfftf,nspden))
    this%vtrial(:,:)=zero
    this%nfftf=nfftf
    this%nspden=nspden
    this%ebcut=zero
    this%edc_kin_freeel=zero
    this%e_kin_freeel=zero
    this%ent_freeel=zero
    this%nfreeel=zero
    this%e_shiftfactor=zero
    if(prt_cg==1) then
      this%prt_cg=.true.
    else
      this%prt_cg=.false.
    end if
    call metric(gmet,gprimd,-1,rmet,rprimd,this%ucvol)
  end subroutine init
  !!***

  !!****f* ABINIT/m_hightemp/destroy
  !! NAME
  !!  destroy
  !!
  !! FUNCTION
  !!  Destroy hightemp_type object, memory deallocation of arrays...
  !!
  !! INPUTS
  !!  this=hightemp_type object concerned
  !!
  !! OUTPUT
  !!  this=hightemp_type object concerned
  !!
  !! PARENTS
  !!  m_gstate
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine destroy(this)

    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this

    ! *********************************************************************

    this%vtrial(:,:)=zero
    ABI_FREE(this%vtrial)
    this%nfftf=0
    this%nspden=0
    this%bcut=0
    this%nbcut=0
    this%version=1
    this%ebcut=zero
    this%edc_kin_freeel=zero
    this%e_kin_freeel=zero
    this%ent_freeel=zero
    this%nfreeel=zero
    this%e_shiftfactor=zero
    this%prt_cg=.false.
    this%ucvol=zero
  end subroutine destroy
  !!***

  !!****f* ABINIT/m_hightemp/compute_e_shiftfactor
  !! NAME
  !!  compute_e_shiftfactor
  !!
  !! FUNCTION
  !!  Compute the energy shift factor $U_0$ corresponding to constant
  !!  potential contribution.
  !!
  !! INPUTS
  !!  this=hightemp_type object concerned
  !!  eigen(mband*nkpt*nsppol)=eigenvalues (hartree)
  !!  eknk(mband*nkpt*nsppol)=kinetic energies (hartree)
  !!  mband=maximum number of bands
  !!  nband(nkpt*nsppol)=desired number of bands at each k point
  !!  nkpt=number of k points
  !!  nsppol=1 for unpolarized, 2 for spin-polarized
  !!  wtk(nkpt)=k point weights
  !!
  !! OUTPUT
  !!  this=hightemp_type object concerned
  !!
  !! PARENTS
  !!  m_vtorho
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine compute_e_shiftfactor(this,eigen,eknk,mband,nband,nkpt,nsppol,wtk)
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
    integer :: band_index,ii,ikpt,isppol,nband_k

    ! *********************************************************************

    ! Compute U_0 from the sum of local
    ! potentials (vtrial), averaging over all space.
    ! Simplest and most precise way to evaluate U_0.
    if(this%version==1) then
      this%e_shiftfactor=sum(this%vtrial)/(this%nfftf*this%nspden)
    end if

    ! Compute U_0^{HEG} from the difference between
    ! eigenvalues and Fermi gas energies, averaged
    ! over lasts nbcut bands.
    if(this%version==2) then
      this%e_shiftfactor=zero
      band_index=0
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          do ii=nband_k-this%nbcut+1,nband_k
            this%e_shiftfactor=this%e_shiftfactor+&
            & wtk(ikpt)*(eigen(band_index+ii)-hightemp_e_heg(dble(ii),this%ucvol))
          end do
          band_index=band_index+nband_k
        end do
      end do
      this%e_shiftfactor=this%e_shiftfactor/this%nbcut
    end if

    ! Compute U_0^K from the difference between
    ! eigenvalues and kinetic energies, averaged
    ! over lasts nbcut bands.
    if(this%version==3) then
      this%e_shiftfactor=zero
      this%ebcut=0
      band_index=0
      do isppol=1,nsppol
        do ikpt=1,nkpt
          nband_k=nband(ikpt+(isppol-1)*nkpt)
          do ii=nband_k-this%nbcut+1,nband_k
            this%e_shiftfactor=this%e_shiftfactor+&
            & wtk(ikpt)*(eigen(band_index+ii)-eknk(band_index+ii))
          end do
          this%ebcut=this%ebcut+wtk(ikpt)*eigen(band_index+nband_k)
          band_index=band_index+nband_k
        end do
      end do
      this%e_shiftfactor=this%e_shiftfactor/this%nbcut
    end if
  end subroutine compute_e_shiftfactor
  !!***

  !!****f* ABINIT/m_hightemp/compute_nfreeel
  !! NAME
  !!  compute_nfreeel
  !!
  !! FUNCTION
  !!  Compute the value of the integral corresponding to the missing
  !!  free electrons contribution after band cut, with an order 1/2
  !!  incomplete Fermi-Dirac integral.
  !!
  !! INPUTS
  !!  this=hightemp_type object concerned
  !!  fermie=chemical potential (Hartree)
  !!  nelect=number of electrons per unit cell
  !!  tsmear=smearing width (or temperature)
  !!
  !! OUTPUT
  !!  this=hightemp_type object concerned
  !!
  !! PARENTS
  !!  m_gstate,m_vtorho,m_occ
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine compute_nfreeel(this,fermie,nelect,tsmear)
    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: fermie,tsmear
    real(dp),intent(inout) :: nelect
    class(hightemp_type),intent(inout) :: this

    ! Local variables -------------------------
    ! Scalars
    integer :: ifftf,ispden
    real(dp) :: factor,gamma,xcut

    ! *********************************************************************

    factor=sqrt(2.)/(PI*PI)*this%ucvol*tsmear**(1.5)

    ! Compute hightemp contribution to nelect integrating
    ! over accessible states from bcut to infinity with
    ! order 1/2 incomplete Fermi-Dirac integral.
    if(this%version==1.or.this%version==2) then
      gamma=(fermie-this%e_shiftfactor)/tsmear
      xcut=hightemp_e_heg(dble(this%bcut),this%ucvol)/tsmear
      nelect=nelect+factor*djp12(xcut,gamma)
    end if

    ! Compute hightemp contribution to nelect integrating
    ! over energy from ebcut to infinity with order 1/2
    ! incomplete Fermi-Dirac integral.
    if(this%version==3) then
      gamma=(fermie-this%e_shiftfactor)/tsmear
      xcut=(this%ebcut-this%e_shiftfactor)/tsmear
      nelect=nelect+factor*djp12(xcut,gamma)
    end if

    ! Compute hightemp contribution to nelect using a sum
    ! of Fermi gas contributions for each point of the fftf grid.
    ! Warning: This is not yet operational. Work in progress.
    if(this%version==4) then
      do ifftf=1,this%nfftf
        do ispden=1,this%nspden
          gamma=(fermie-this%vtrial(ifftf,ispden))/tsmear
          xcut=hightemp_e_heg(dble(this%bcut),this%ucvol)/tsmear
          nelect=nelect+factor*djp12(xcut,gamma)/(this%nfftf*this%nspden)
        end do
      end do
    end if
  end subroutine compute_nfreeel
  !!***

  !!****f* ABINIT/m_hightemp/compute_efreeel
  !! NAME
  !!  compute_efreeel
  !!
  !! FUNCTION
  !!  Compute the value of the integral corresponding to the missing
  !!  kinetic energy contribution of free electrons after band cut,
  !!  with an order 3/2 incomplete Fermi-Dirac integral.
  !!
  !! INPUTS
  !!  this=hightemp_type object concerned
  !!  fermie=chemical potential (Hartree)
  !!  nfftf=number of points in the fine FFT mesh (for this processor)
  !!  nspden=number of spin-density components
  !!  tsmear=smearing width (or temperature)
  !!  vtrial(nfft,nspden)= trial potential (Hartree)
  !!
  !! OUTPUT
  !!  this=hightemp_type object concerned
  !!
  !! PARENTS
  !!  m_gstate,m_vtorho
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine compute_efreeel(this,fermie,nfftf,nspden,tsmear,vtrial)
    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    real(dp),intent(in) :: fermie,tsmear
    integer,intent(in) :: nfftf,nspden
    ! Arrays
    real(dp),intent(in) :: vtrial(nfftf,nspden)

    ! Local variables -------------------------
    ! Scalars
    integer :: ifftf,ispden
    real(dp) :: factor,gamma,xcut

    ! *********************************************************************

    this%e_kin_freeel=zero
    factor=sqrt(2.)/(PI*PI)*this%ucvol*tsmear**(2.5)

    ! Compute hightemp contribution to kinetic energy integrating
    ! over accessible states from bcut to infinity with
    ! order 3/2 incomplete Fermi-Dirac integral.
    if(this%version==1.or.this%version==2) then
      gamma=(fermie-this%e_shiftfactor)/tsmear
      xcut=hightemp_e_heg(dble(this%bcut),this%ucvol)/tsmear
      this%e_kin_freeel=this%e_kin_freeel+factor*djp32(xcut,gamma)
    end if

    ! Compute hightemp contribution to kinetic energy integrating
    ! over energy from ebcut to infinity with order 3/2
    ! incomplete Fermi-Dirac integral.
    if(this%version==3) then
      gamma=(fermie-this%e_shiftfactor)/tsmear
      xcut=(this%ebcut-this%e_shiftfactor)/tsmear
      this%e_kin_freeel=this%e_kin_freeel+factor*djp32(xcut,gamma)
    end if

    ! Compute hightemp contribution to kinetic energy using a sum
    ! of Fermi gas contributions for each point of the fftf grid.
    ! Warning: This is not yet operational. Work in progress.
    if(this%version==4) then
      do ifftf=1,this%nfftf
        do ispden=1,this%nspden
          gamma=(fermie-this%vtrial(ifftf,ispden))/tsmear
          xcut=hightemp_e_heg(dble(this%bcut),this%ucvol)/tsmear
          this%e_kin_freeel=this%e_kin_freeel+factor*djp32(xcut,gamma)/&
          & (this%nfftf*this%nspden)
        end do
      end do
    end if

    ! Compute the double counting term of the contribution to
    ! kinetic energy from the vtrial potential.
    this%edc_kin_freeel=zero
    do ispden=1,nspden
      do ifftf=1,nfftf
        this%edc_kin_freeel=this%edc_kin_freeel+vtrial(ifftf,ispden)
      end do
    end do
    this%edc_kin_freeel=this%edc_kin_freeel*this%nfreeel/nspden/nfftf/nspden
  end subroutine compute_efreeel
  !!***

  !!****f* ABINIT/m_hightemp/compute_ent_freeel
  !! NAME
  !!  compute_ent_freeel
  !!
  !! FUNCTION
  !!  Compute the value of the integral corresponding to the missing
  !!  entropy contribution of free electrons after band cut using
  !!  incomplete Fermi-Dirac integrals.
  !!
  !! INPUTS
  !!  this=hightemp_type object concerned
  !!  fermie=chemical potential (Hartree)
  !!  tsmear=smearing width (or temperature)
  !!
  !! OUTPUT
  !!  this=hightemp_type object concerned
  !!
  !! PARENTS
  !!  m_vtorho
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine compute_ent_freeel(this,fermie,tsmear)
    ! Arguments -------------------------------
    ! Scalars
    class(hightemp_type),intent(inout) :: this
    real(dp),intent(in) :: fermie,tsmear

    ! Local variables -------------------------
    ! Scalars
    integer :: ii,ifftf,ispden
    real(dp) :: ix,step,factor,fn,gamma,minocc
    ! Arrays
    real(dp),dimension(:),allocatable :: valuesent

    ! *********************************************************************

    this%ent_freeel=zero

    ! Compute hightemp contribution to the entropy integrating
    ! over accessible states with Fermi-Dirac complete integrals and
    ! substracting 0 to bcut contribution with numeric integration.
    if(this%version==1.or.this%version==2) then
      factor=sqrt(2.)/(PI*PI)*this%ucvol*tsmear**(2.5)
      gamma=(fermie-this%e_shiftfactor)/tsmear
      ABI_MALLOC(valuesent,(this%bcut+1))

      !$OMP PARALLEL DO PRIVATE(fn,ix) SHARED(valuesent)
      do ii=1,this%bcut+1
        ix=dble(ii)-one
        fn=fermi_dirac(hightemp_e_heg(ix,this%ucvol)+this%e_shiftfactor,fermie,tsmear)
        if(fn>tol16.and.(one-fn)>tol16) then
          valuesent(ii)=-two*(fn*log(fn)+(one-fn)*log(one-fn))
        else
          valuesent(ii)=zero
        end if
      end do
      !$OMP END PARALLEL DO

      ! We need at least 6 elements in valuesent to call simpson function.
      if(size(valuesent)>=6) then
        this%ent_freeel=5./3.*factor*dip32(gamma)/tsmear-&
        & gamma*factor*dip12(gamma)/tsmear-&
        simpson(one,valuesent)
      end if
      ABI_FREE(valuesent)
    end if

    ! Compute hightemp contribution to the entropy integrating
    ! over energy with Fermi-Dirac complete integrals and
    ! substracting 0 to bcut contribution with numeric integration.
    if(this%version==3) then
      factor=sqrt(2.)/(PI*PI)*this%ucvol*tsmear**(2.5)
      gamma=(fermie-this%e_shiftfactor)/tsmear
      ABI_MALLOC(valuesent,(this%bcut+1))

      step=(this%ebcut-this%e_shiftfactor)/(this%bcut)
      !$OMP PARALLEL DO PRIVATE(fn,ix) SHARED(valuesent)
      do ii=1,this%bcut+1
        ix=this%e_shiftfactor+(dble(ii)-one)*step
        fn=fermi_dirac(ix,fermie,tsmear)
        if(fn>tol16.and.(one-fn)>tol16) then
          valuesent(ii)=-(fn*log(fn)+(one-fn)*log(one-fn))*&
          & hightemp_dosfreeel(ix,this%e_shiftfactor,this%ucvol)
        else
          valuesent(ii)=zero
        end if
      end do
      !$OMP END PARALLEL DO

      ! We need at least 6 elements in valuesent to call simpson function.
      if(size(valuesent)>=6) then
        this%ent_freeel=5./3.*factor*dip32(gamma)/tsmear-&
        & gamma*factor*dip12(gamma)/tsmear-&
        simpson(step,valuesent)
      end if
      ABI_FREE(valuesent)
    end if

    ! Compute hightemp contribution to the entropy using a sum
    ! of Fermi gas contributions for each point of the fftf grid,
    ! as we do for version=1 and version=2.
    ! Warning: This is not yet operational. Work in progress.
    if(this%version==4) then
      step=one
      do ifftf=1,this%nfftf
        do ispden=1,this%nspden
          ! Dynamic array find size
          ix=dble(this%bcut)
          ii=0
          fn=fermi_dirac(hightemp_e_heg(ix,this%ucvol)+&
          & this%vtrial(ifftf,ispden),fermie,tsmear)
          minocc=tol16
          do while(fn>minocc)
            fn=fermi_dirac(hightemp_e_heg(ix,this%ucvol)+&
            & this%vtrial(ifftf,ispden),fermie,tsmear)
            ii=ii+1
            ix=ix+step
          end do
          ABI_MALLOC(valuesent,(ii))
          ix=dble(this%bcut)
          ii=0
          fn=fermi_dirac(hightemp_e_heg(ix,this%ucvol)+&
          & this%vtrial(ifftf,ispden),fermie,tsmear)
          do while(fn>minocc)
            fn=fermi_dirac(hightemp_e_heg(ix,this%ucvol)+&
            & this%vtrial(ifftf,ispden),fermie,tsmear)
            ii=ii+1
            valuesent(ii)=-2*(fn*log(fn)+(1.-fn)*log(1.-fn))
            ix=ix+step
          end do
          if (ii>1) then
            this%ent_freeel=this%ent_freeel+simpson(step,valuesent)/&
            & (this%nfftf*this%nspden)
          end if
          ABI_FREE(valuesent)
        end do
      end do
    end if
  end subroutine compute_ent_freeel
  !!***

  !!****f* ABINIT/m_hightemp/compute_pw_avg_std
  !! NAME
  !!  compute_pw_avg_std
  !!
  !! FUNCTION
  !!  Debug subroutine to print wave functions coefficients in a precise way.
  !!
  !! INPUTS
  !!  this=hightemp_type object concerned
  !!  cg(2,mpw*dtset%nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
  !!  eig_k(nband_k)=array for holding eigenvalues (hartree)
  !!  ek_k(nband_k)=contribution from each band to kinetic energy, at this k-point
  !!  fnameabo=filename for printing of the planewaves coefficients
  !!  gprimd(3,3)=dimensional reciprocal space primitive translations
  !!  icg=shift to be applied on the location of data in the array cg
  !!  ikpt=number of the k-point
  !!  istwf_k=parameter that describes the storage of wfs
  !!  kg_k(3,npw)=integer coordinates of G vectors in basis sphere
  !!  kinpw(npw_k)=(modified) kinetic energy for each plane wave (Hartree)
  !!  kpt(3,nkpt)=reduced coordinates of k points.
  !!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
  !!  mpi_enreg=information about MPI parallelization
  !!  nband_k=number of bands at each k point
  !!  nkpt=number of k points.
  !!  npw_k=number of plane waves at this k point
  !!  nspinor=Number of spinor components
  !!  wtk(nkpt)=weight assigned to each k point
  !!
  !! OUTPUT
  !!  this=hightemp_type object concerned
  !!
  !! PARENTS
  !!  m_vtorho
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine compute_pw_avg_std(this,cg,eig_k,ek_k,fnameabo,&
  & gprimd,icg,ikpt,istwf_k,kg_k,kinpw,kpt,mcg,mpi_enreg,nband_k,&
  & nkpt,npw_k,nspinor,wtk)
    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: icg,ikpt,istwf_k,mcg,nband_k,nkpt,npw_k,nspinor
    real(dp),intent(in) :: wtk
    class(hightemp_type),intent(inout) :: this
    type(MPI_type),intent(in) :: mpi_enreg
    character(len=*),intent(in) :: fnameabo
    ! Arrays
    integer,intent(in) :: kg_k(3,npw_k)
    real(dp),intent(in) :: eig_k(nband_k),ek_k(nband_k),gprimd(3,3)
    real(dp),intent(in) :: kpt(3,nkpt),kinpw(npw_k),cg(2,mcg)

    ! Local variables -------------------------
    ! Scalars
    integer :: blocksize,i1,iband,iblock,iblocksize,ipw,nblockbd,unit
    real(dp) :: kpg1,kpg2,kpg3
    character(len=50) :: filenameoutpw
    ! Arrays
    real(dp) :: cgnk(2,npw_k*nspinor),cgnk2(npw_k*nspinor)

    ! *********************************************************************

    unit=96+mpi_enreg%me
    ! Debug: Printing Eigenvalues, Kinetic and PW grid if requested
    if(this%prt_cg.and.(mpi_enreg%me==mpi_enreg%me_kpt)) then
      write(filenameoutpw,'(A,I5.5)') '_PW_MESH_k',ikpt
      open(file=trim(fnameabo)//trim(filenameoutpw),unit=unit)
      do ipw=1,npw_k
        kpg1=kpt(1,ikpt)+dble(kg_k(1,ipw))
        kpg2=kpt(2,ikpt)+dble(kg_k(2,ipw))
        kpg3=kpt(3,ikpt)+dble(kg_k(3,ipw))

        write(unit,'(i14,ES13.5,ES13.5,ES13.5,ES13.5)')&
        & ipw,&
        & gprimd(1,1)*kpg1+gprimd(1,2)*kpg2+gprimd(1,3)*kpg3,&
        & gprimd(2,1)*kpg1+gprimd(2,2)*kpg2+gprimd(2,3)*kpg3,&
        & gprimd(3,1)*kpg1+gprimd(3,2)*kpg2+gprimd(3,3)*kpg3,&
        & kinpw(ipw)
      end do
      close(unit)

      write(filenameoutpw,'(A,I5.5)') '_PW_EIG_k',ikpt
      open(file=trim(fnameabo)//trim(filenameoutpw),unit=unit)
      write(unit,'(ES12.5)') eig_k
      close(unit)

      write(filenameoutpw,'(A,I5.5)') '_PW_KIN_k',ikpt
      open(file=trim(fnameabo)//trim(filenameoutpw),unit=unit)
      write(unit,'(ES12.5)') ek_k
      close(unit)
    end if

    ! Parallelism over FFT and/or bands: define sizes and tabs
    nblockbd=nband_k/(mpi_enreg%nproc_band*mpi_enreg%bandpp)
    blocksize=nband_k/nblockbd

    ! Loop over bands or blocks of bands. Note that in sequential mode
    ! iblock=iband, nblockbd=nband_k and blocksize=1
    do iblock=1,nblockbd
      do iblocksize=1,blocksize
        iband=(iblock-1)*blocksize+iblocksize

        ! Not sure of this debug block if we parallelise over bands.
        ! Forcing disable if npband > 1.
        if(this%prt_cg.and.mpi_enreg%nproc_band==1) then
          cgnk(:,:)=cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
          cgnk2(:)=(cgnk(1,:)*cgnk(1,:)+cgnk(2,:)*cgnk(2,:))

          write(filenameoutpw, '(A,I5.5,A,I5.5)') '_PW_k',ikpt,'_b',iband
          open(file=trim(fnameabo)//trim(filenameoutpw),unit=unit)

          ! if(istwf_k==1)then
          !   do ipw=1,npw_k
          !     write(unit,'(i14,ES14.6,ES14.6,i14)') ipw,kinpw(ipw),&
          !     & sqrt(cgnk(1,ipw)*cgnk(1,ipw)+cgnk(2,ipw)*cgnk(2,ipw))
          !   end do
          ! else if(istwf_k>=2)then
          !   i1=1
          !   if(istwf_k==2.and.mpi_enreg%me_g0==1)then ! MPIWF need to know which proc has G=0
          !     write(unit,'(i14,ES14.6,ES14.6,i14)') ipw,kinpw(ipw),&
          !     & sqrt(cgnk(1,1)*cgnk(1,1))
          !     i1=2
          !   end if
          !   do ipw=i1,npw_k
          !     write(unit,'(i14,ES14.6,ES14.6,i14)') ipw,kinpw(ipw),&
          !     & sqrt(cgnk(1,ipw)*cgnk(1,ipw)+cgnk(2,ipw)*cgnk(2,ipw))
          !     ! write(unit,'(i14,ES14.6,ES14.6,i14)') ipw,kinpw(ipw),&
          !     ! & sqrt(cgnk(1,ipw)*cgnk(1,ipw)+cgnk(2,ipw)*cgnk(2,ipw))
          !   end do
          ! end if

          do ipw=1,npw_k
            write(unit,'(i14,ES14.6,ES14.6,i14)') ipw,&
            & kinpw(ipw),sqrt(cgnk2(ipw))
          end do

          close(unit)
        end if
      end do
    end do
  end subroutine compute_pw_avg_std
  !!***

  !----------------------------------------------------------------------

  !!****f* ABINIT/m_hightemp/dip12
  !! NAME
  !!  dip12
  !!
  !! FUNCTION
  !!  Returns the complete Fermi integral of order 1/2.
  !!
  !! INPUTS
  !!  gamma=complete Fermi integral argument
  !!
  !! OUTPUT
  !!  dip12=resulting function
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function dip12(gamma)
    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: gamma

    ! Local variables -------------------------
    ! Scalars
    real(dp) :: d,dip12,dy

    ! *********************************************************************

    if (gamma.lt.3.) then
      dy=exp(gamma)
      if (gamma+1.9375.LE.0) then
        dip12=dy*&
        & (1.-dy*(0.35355283-dy*(0.19242767-dy*(0.12456909-dy*&
        & (8.5114507E-02-dy*4.551794E-02)))))
      else
        d=gamma-0.5
        dip12=dy*(0.677695804-d*(0.187773135+d*(2.16197521E-02-d*&
        & (9.23703807E-03+d*&
        & (1.71735167E-03-d*(6.07913775E-04+d*&
        & (1.1448629E-04-d*&
        & (4.544432E-05+d*(6.4719368E-06-d*(3.794983E-06+d*&
        & (1.7338029E-07-d*&
        & (3.5546516E-07-d*(3.7329191E-08+d*&
        & (3.3097822E-08-d*&
        & (8.3190193E-09+d*(2.2752769E-09-d*(7.836005E-10+d*&
        & (7.519551E-11-d*2.960006E-11))))))))))))))))))
      end if
    else if (gamma.lt.20.) then
      if (gamma.lt.10.) then
        d=gamma-6.5
        dip12=12.839811+d*&
        & (2.844774+d*(0.114920926-d*(3.43733039E-03-d*&
        & (2.3980356E-04-d*&
        & (2.0201888E-05-d*(1.5219883E-06-d*&
        & (6.2770524E-08+d*&
        & (4.8830336E-09-d*(2.1031164E-09-d*(5.785753E-10-d*&
        & (7.233066E-11-d*1.230727E-12)))))))))))
      else
        d=gamma-14.5
        dip12=41.7799227+d*&
        & (4.2881461+d*(7.45407825E-02-d*(8.79243296E-04-d*&
        & (2.38288861E-05-d*&
        & (8.82474867E-07-d*(3.82865217E-08-d*&
        & (1.9274292E-09-d*(1.42248669E-10-d*8.17019813E-12))))))))
      end if
    else
      d=1./gamma
      dy=gamma*sqrt(gamma)/1.329340388
      dip12=dy*(1.-d*(9.354E-07-d*(1.2338391-d*(6.77931E-03-d*&
      & 1.17871643))))
    end if
    dip12=dip12*0.88622692
  end function dip12
  !!***

  !!****f* ABINIT/m_hightemp/djp12
  !! NAME
  !!  djp12
  !!
  !! FUNCTION
  !!  Returns the incomplete Fermi integral of order 1/2.
  !!
  !! INPUTS
  !!  xcut=lower bound of the incomplete Fermi integral
  !!  gamma=incomplete Fermi integral argument
  !!
  !! OUTPUT
  !!  djp12=resulting function
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function djp12(xcut,gamma)
    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: xcut,gamma
    real(dp) :: djp12

    ! Local variables -------------------------
    ! Scalars
    real(dp) :: d2h,db,dc,dd,de,df1,df2,df3,dh
    real(dp) :: ds,dt,dv,dw,dxm,dxp
    integer :: i,ind,iq,k,nm,np,nq
    ! Arrays
    real(dp) :: dq(5),df(101),dy(101)

    ! *********************************************************************

    dh=0.2D+0
    d2h=0.4D+0
    nm=101
    ind=0
    dq=(/1.D+0,2.828427124D+0,5.196152423D+0,8.D+0,1.118033989D+1/)

    djp12=0.D+0
    dxm=gamma-1.5D+1
    if (xcut.gt.dxm) then
      if (ind.eq.0) then
        do i=1,nm
          dy(i)=-1.5D+1+(i-1)*dh
          df(i)=1.D+0+dexp(dy(i))
        end do
        ind=1
      end if
      dxp=gamma+5.D+0
      if (xcut.lt.dxp) then
        dc=dxp
      else
        dc=xcut
      end if
      db=dexp(gamma-dc)
      dt=db
      do iq=1,5
        dd=iq*dc
        ds=dsqrt(dd)
        dw=1.+.3275911*ds
        dw=1.D+0/dw
        dv=dw*(.2258368458D+0+&
        & dw*(-.2521286676D+0+dw*(1.2596951294D+0+&
        & dw*(-1.2878224530D+0+dw*(.9406460699D+0)))))
        dv=dv+ds
        de=dt*dv/dq(iq)
        djp12=djp12+de
        if (dabs(de).lt.(1.D-07*djp12)) exit
        dt=-dt*db
      end do
      if (xcut.ge.dxp) return
      np=(dxp-xcut)/dh
      np=2*(np/2)
      np=nm-np
      nq=(15.-gamma)/dh
      nq=1+2*(nq/2)
      if (np.lt.nq) np=nq
      if (np.le.nm) then
        df3=0.D+0
        dt=dy(np)+gamma
        dv=(dt-xcut)/2.D+0
        df3=0.D0
        if (dt.ge.1.D-13) df3=dsqrt(dt)
        if (dabs(dv).ge.1.D-13) then
          df1=dsqrt(xcut)
          dt=df1+df3
          dw=(dv+dv)/(dt*dt)
          df2=dw*dw
          df2=df2+df2
          db=df2*(df2+7.D+0)+7.D+1
          dc=7.D+0*(1.D+1-df2)
          dc=dc*dw
          dd=-df2*(df2-2.8D+1)+1.4D+2
          dd=dd+dd
          ds=dt*((db-dc)/(1.D+0+dexp(xcut-gamma))+&
          & dd/(1.D+0+dexp(xcut+dv-gamma))+(db+dc)/df(np))
          ds=ds*dv/4.2D+2
          djp12=djp12+ds
        end if
        if (np.ne.nm) then
          ds=0.D+0
          np=np+2
          do k=np,nm,2
            df1=df3
            df3=dsqrt(dy(k)+gamma)
            dt=df1+df3
            dw=d2h/(dt*dt)
            df2=dw*dw
            df2=df2+df2
            db=df2*(df2+7.D+0)+7.D+1
            dc=7.D+0*(1.D+1-df2)
            dc=dc*dw
            dd=-df2*(df2-2.8D+1)+1.4D+2
            dd=dd+dd
            ds=ds+dt*((db-dc)/df(k-2)+dd/df(k-1)+(db+dc)/df(k))
          end do
          ds=ds*dh/4.2D+2
          djp12=djp12+ds
        end if
        if (xcut.ge.dxm) return
      end if
    end if
    djp12=dip12(gamma)-xcut*dsqrt(xcut)/1.5D+0
  end function djp12
  !!***

  !!****f* ABINIT/m_hightemp/dip32
  !! NAME
  !!  dip32
  !!
  !! FUNCTION
  !!  Returns the complete Fermi integral of order 3/2.
  !!
  !! INPUTS
  !!  gamma=complete Fermi integral argument
  !!
  !! OUTPUT
  !!  dip32=resulting function
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function dip32(gamma)
    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: gamma

    ! Local variables -------------------------
    ! Scalars
    real(dp) :: d,dip32,dval

    ! *********************************************************************

    if (gamma.GT.1.75) then
      dval=gamma*gamma*SQRT(gamma)
      if (gamma.LT.4.5) then
        d=gamma-3.125
        dip32=(1.27623+0.596065*gamma+0.3*dval)*&
        & (1.0055558385-d*(5.23889494E-03+d*&
        & (3.13523144E-03-d*(3.06124286E-03-d*&
        & (1.3644667E-03-d*&
        & (4.1528384E-04-d*(8.901188E-05-d*(1.079979E-05+d*&
        & (2.29058E-06-d*(2.58985E-06-d*7.30909E-07))))))))))
      else if (gamma.LT.12.) then
        if (gamma.LT.8.) then
          d=gamma-6.25
          dip32=(2.01508+0.425775*gamma+0.3*dval)*&
          & (1.000387131-d*(3.93626295E-04+d*&
          & (2.55710115E-04-d*&
          & (1.57383494E-04-d*(5.0286036E-05-d*&
          & (1.2073559865E-05-d*&
          & (2.4909523213E-06-d*(5.244328548E-07-d*&
          & 8.0884033896E-08))))))))
        else
          d=gamma-10.
          dip32=0.3*dval*&
          & (1.064687247-d*(1.22972303E-02-d*(1.8362121E-03-d*&
          & (2.433558E-04-d*(3.018186E-05-d*(3.5694E-06-d*&
          & (4.11212E-07-d*(5.2151E-08-d*5.8424E-09))))))))
        end if
      else
        d=1./gamma
        dip32=0.30090111127*dval*&
        & (1.-d*(2.863E-06-d*(6.168876549-d*&
        & (1.740553E-02+d*(1.425257+d*2.95887)))))
      end if
    else if (gamma+0.75.LE.0) then
      d=EXP(gamma)
      dip32=d*&
      & (1.-d*(1.76775246E-01-d*(6.4124584E-02-d*&
      & (3.1027055E-02-d*(1.6797637E-02-d*&
      & (8.212636E-03-d*(2.384106E-03)))))))
    else
      d=gamma-0.5
      dip32=EXP(gamma)*(0.846691-0.128948*gamma)*&
      & (1.034064158+d*(2.778947E-03-d*&
      & (3.572502805E-02+d*(3.0411645E-03-d*&
      & (1.7380548E-03+d*(2.7756776E-04-d*&
      & (8.08302E-05+d*(1.59606E-05-d*&
      & (3.8144E-06+d*7.4446E-07)))))))))
    end if
    dip32=dip32*1.32934038
  end function dip32
  !!***

  !!****f* ABINIT/m_hightemp/djp32
  !! NAME
  !!  djp32
  !!
  !! FUNCTION
  !!  Returns the incomplete Fermi integral of order 3/2.
  !!
  !! INPUTS
  !!  xcut=lower bound of the incomplete Fermi integral
  !!  gamma=incomplete Fermi integral argument
  !!
  !! OUTPUT
  !!  djp32=resulting function
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function djp32(xcut,gamma)
    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: xcut,gamma
    real(dp) :: djp32

    ! Local variables -------------------------
    ! Scalars
    real(dp) :: d2h,db,dc,dd,de,df1,df2,df3,dh
    real(dp) :: ds,dt,dv,dw,dx1,dx2
    real(dp) :: dx3,dxm,dxp
    integer :: i,ind,iq,k,nm,np,nq
    ! Arrays
    real(dp) :: dq(5),df(101),dy(101)

    ! *********************************************************************

    dh=0.2D+0
    d2h=0.4D+0
    nm=101
    ind=0
    dq=(/1.D+0,5.656854228D+0,1.558845727D+1,3.2D+1,5.590169945D+0/)

    djp32=0.D+0
    dxm=gamma-1.5D+1
    if (xcut.GT.dxm) then
      if (ind.EQ.0) then
        do i=1,nm
          dy(i)=-1.5D+1+(i-1)*dh
          df(i)=1.D+0+DEXP(dy(i))
        end do
        ind=1
      end if
      dxp=gamma+5.D+0
      if (xcut.LT.dxp) then
        dc=dxp
      else
        dc=xcut
      end if
      db=DEXP(gamma-dc)
      dt=db
      do iq=1,5
        dd=iq*dc
        ds=DSQRT(dd)
        dw=1.+.3275911*ds
        dw=1.D+0/dw
        dv=dw*(.2258368458D+0+&
        & dw*(-.2521286676D+0+dw*(1.2596951294D+0+&
        & dw*(-1.2878224530D+0+dw*(.9406460699D+0)))))
        dv=dv+ds
        dv=1.5D+0*dv+ds*dd
        de=dt*dv/dq(iq)
        djp32=djp32+de
        if (DABS(de).LT.(1.D-07*djp32)) exit
        dt=-dt*db
      end do
      if (xcut.GE.dxp) return
      np=(dxp-xcut)/dh
      np=2*(np/2)
      np=nm-np
      nq=(15.-gamma)/dh
      nq=1+2*(nq/2)
      if (np.LT.nq) np=nq
      if (np.LE.nm) then
        df3=0.D+0
        dt=dy(np)+gamma
        dv=(dt-xcut)/2.D+0
        df3=DSQRT(dt)
        dx3=dt
        if (DABS(dv).GE.1.D-13) then
          df1=DSQRT(xcut)
          dt=df1+df3
          dw=(dv+dv)/(dt*dt)
          df2=dw*dw
          df2=df2+df2
          db=df2*(df2+7.D+0)+7.D+1
          dc=7.D+0*(1.D+1-df2)
          dc=dc*dw
          dd=-df2*(df2-2.8D+1)+1.4D+2
          dd=dd+dd
          ds=dt*((db-dc)*xcut/(1.D+0+DEXP(xcut-gamma))+dd*(xcut+dv)&
          & /(1.D+0+DEXP(xcut+dv-gamma))+(db+dc)*(dy(np)+gamma)/df(np))
          ds=ds*dv/4.2D+2
          djp32=djp32+ds
        end if
        if (np.NE.nm) then
          ds=0.D+0
          np=np+2
          do k=np,nm,2
            dx1=dx3
            df1=df3
            dx2=dy(k-1)+gamma
            dx3=dy(k)+gamma
            df3=DSQRT(dx3)
            dt=df1+df3
            dw=d2h/(dt*dt)
            df2=dw*dw
            df2=df2+df2
            db=df2*(df2+7.D+0)+7.D+1
            dc=7.D+0*(1.D+1-df2)
            dc=dc*dw
            dd=-df2*(df2-2.8D+1)+1.4D+2
            dd=dd+dd
            ds=ds+dt*((db-dc)*dx1/df(k-2)+dd*dx2/df(k-1)&
            & +(db+dc)*dx3/df(k))
          end do
          ds=ds*dh/4.2D+2
          djp32=djp32+ds
        end if
        if (xcut.GE.dxm) return
      end if
    end if
    djp32=dip32(gamma)-xcut*xcut*DSQRT(xcut)/2.5D+0
  end function djp32
  !!***

  !!****f* ABINIT/m_hightemp/hightemp_dosfreeel
  !! NAME
  !!  hightemp_dosfreeel
  !!
  !! FUNCTION
  !!  Returns the free particle density of states for a given energy.
  !!
  !! INPUTS
  !!  energy=get the value of the free particle density of states at this energy
  !!  e_shiftfactor=energy shift factor
  !!  ucvol=unit cell volume (bohr^3)
  !!
  !! OUTPUT
  !!  hightemp_dosfreeel=value of free particle density of states at given energy
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
  !!***

  !!****f* ABINIT/m_hightemp/hightemp_e_heg
  !! NAME
  !!  hightemp_e_heg
  !!
  !! FUNCTION
  !!  Returns the energy of the Fermi gas for a given number of
  !!  accessible states.
  !!
  !! INPUTS
  !!  iband=number of accessible states
  !!  ucvol=unit cell volume (bohr^3)
  !!
  !! OUTPUT
  !!  hightemp_e_heg=energy of homogeneous electron gas for a given number of accessible states
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  function hightemp_e_heg(iband,ucvol)
    ! Arguments -------------------------------
    ! Scalars
    real(dp),intent(in) :: iband,ucvol
    real(dp) :: hightemp_e_heg

    ! *********************************************************************

    hightemp_e_heg=.5*(iband*6*PI*PI/ucvol)**(2./3.)
  end function hightemp_e_heg
  !!***

  !!****f* ABINIT/m_hightemp/hightemp_prt_eigocc
  !! NAME
  !!  hightemp_prt_eigocc
  !!
  !! FUNCTION
  !!  Printing in a _EIGOCC file many informations like Fermi energy, Unit cell volume,
  !!  electronic temperature (tsmear when occopt=3), eigenvalues and occupation for each k-point,
  !!  weight of the k-point, number of bands, number of k-points and more...
  !!  This file is intended to be used for custom DOS computation with external tools for example.
  !!
  !! INPUTS
  !!  e_kin_freeel=kinetic energy contribution of free electrons after band cut
  !!  e_shiftfactor=energy shift factor
  !!  eigen(mband*nkpt*nsppol)=eigenvalues (hartree)
  !!  etotal=total energy
  !!  energies=object containing all energies
  !!  fnameabo=filename for printing of the planewaves coefficients
  !!  iout=out file number
  !!  iter=MD iteration step
  !!  kptns(3,nkpt),kptns0(3,nkpt0)=k point sets (reduced coordinates)
  !!  mband=maximum number of bands
  !!  nband(nkpt*nsppol)=desired number of bands at each k point
  !!  nfreeel=number free electrons
  !!  nkpt=number of k points
  !!  nsppol=1 for unpolarized, 2 for spin-polarized
  !!  occ(mband*nkpt*nsppol) = occupation number for each band and k
  !!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
  !!  tsmear=temperature (Hartree)
  !!  usepaw=1 if PAW, 0 otherwise
  !!  wtk(nkpt)=k point weights
  !!  strten(6)=components of the stress tensor (hartree/bohr^3)
  !!  istep=SCF electronic iteration number
  !!
  !! OUTPUT
  !!
  !! PARENTS
  !!
  !! CHILDREN
  !!
  !! SOURCE
  subroutine hightemp_prt_eigocc(e_kin_freeel,e_shiftfactor,eigen,etotal,energies,fnameabo,&
  & iout,iter,kptns,mband,nband,nfreeel,nkpt,nsppol,occ,rprimd,tsmear,usepaw,wtk,&
  & strten,istep) ! Optional arguments
    ! Arguments -------------------------------
    ! Scalars
    integer,intent(in) :: iout,iter,mband,nkpt,nsppol,usepaw
    real(dp),intent(in) :: e_kin_freeel,e_shiftfactor,etotal,nfreeel,tsmear
    character(len=*),intent(in) :: fnameabo
    type(energies_type),intent(in) :: energies
    integer,intent(in),optional :: istep
    ! Arrays
    integer,intent(in) :: nband(nkpt*nsppol)
    real(dp),intent(in) :: eigen(mband*nkpt*nsppol),kptns(3,nkpt)
    real(dp),intent(in) :: rprimd(3,3)
    real(dp),intent(in) :: occ(mband*nkpt*nsppol)
    real(dp),intent(in) :: wtk(nkpt)
    real(dp),intent(in),optional :: strten(6)

    ! Local variables -------------------------
    ! Scalars
    integer :: band_index,iband,ii,ikpt,isppol,nband_k,temp_unit
    real(dp) :: ucvol,pressure
    character(len=200) :: fnameabo_eigocc
    character(len=500) :: msg
    ! Arrays
    real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)

    ! *********************************************************************

    band_index=0
    call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

    fnameabo_eigocc=trim(fnameabo) // trim('_EIGOCC')
    write(msg,'(a,a)') ' prt_eigocc : about to open file ',trim(fnameabo_eigocc)
    call wrtout(iout,msg,'COLL')
    if (open_file(fnameabo_eigocc,msg,newunit=temp_unit,status='unknown',form='formatted') /= 0) then
      ABI_ERROR(msg)
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
    & ' Kinetic energy     = ',energies%e_kinetic,' Ha         Kin. free el. E    = ',e_kin_freeel,&
    & ' Ha         Total kin. energy  = ',energies%e_kinetic+e_kin_freeel, ' Ha'
    call wrtout(temp_unit,msg,'COLL')

    write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
    & ' Hartree energy     = ',energies%e_hartree,' Ha         Ewald energy       = ',energies%e_ewald,&
    & ' Ha         E. Shift factor U0 = ',e_shiftfactor, ' Ha'
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
    if(present(strten)) then
      pressure=-(strten(1)+strten(2)+strten(3))*HaBohr3_GPa/3.0_dp
    else
      pressure=zero
    end if
    write(msg, '(a,ES12.5,a,ES12.5,a,ES15.8,a)') &
    & ' Unit cell vol      = ',ucvol,' Bohr^3     Elec. temperature  = ',tsmear,&
    & ' Ha         Pressure           = ',pressure,' GPa'
    call wrtout(temp_unit,msg,'COLL')

    write(msg, '(a,i12,a,i12,a,ES15.8)') &
    & ' N. elec. iter      = ',istep,'            MD iteration       = ',iter,&
    & '            N. free electrons  = ',nfreeel
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
  !!***

  subroutine prt_cg_from_wf()
    write(0,*) 'TEST'
    ! !BLANCHET TEMPORARY
    ! cgshift=(iband-1)*npw_k*nspinor + (cspinor-1)*npw_k
    ! write(0,*) 'TEST'
    ! write(filenameoutpw,'(A,I5.5)') '_PW_MESH_k',ckpt
    ! open(file=trim(wfk_fname)//trim(filenameoutpw),unit=96)
    ! do ipw=1,npw_k
    !   kpg1=kpt(1,ikpt)+dble(kg_k(1,ipw))
    !   kpg2=kpt(2,ikpt)+dble(kg_k(2,ipw))
    !   kpg3=kpt(3,ikpt)+dble(kg_k(3,ipw))
    !   write(96,'(i14,ES13.5,ES13.5,ES13.5,ES13.5)')&
    !   & ipw,&
    !   & gprimd(1,1)*kpg1+gprimd(1,2)*kpg2+gprimd(1,3)*kpg3,&
    !   & gprimd(2,1)*kpg1+gprimd(2,2)*kpg2+gprimd(2,3)*kpg3,&
    !   & gprimd(3,1)*kpg1+gprimd(3,2)*kpg2+gprimd(3,3)*kpg3,&
    !   & half*(two_pi*kpgnorm(ipw))**2
    ! end do
    ! close(96)
    ! ABI_MALLOC(cgnk,(2,npw_k*nspinor))
    ! ABI_MALLOC(cgnk2,(npw_k*nspinor))
    ! do iband=1,nband(ckpt)
    !   cgshift=(iband-1)*npw_k*nspinor + (cspinor-1)*npw_k
    !   cgnk(:,:)=cg_k(:,cgshift+1:cgshift+npw_k)
    !   cgnk2(:)=(cgnk(1,:)*cgnk(1,:)+cgnk(2,:)*cgnk(2,:))
    !   write(filenameoutpw, '(A,I5.5,A,I5.5)') '_PW_k',ckpt,'_b',cband
    !   open(file=trim(wfk_fname)//trim(filenameoutpw),unit=96)
    !   do ipw=1,npw_k
    !     write(96,'(i14,ES14.6,ES14.6,i14)') ipw,&
    !     & half*(two_pi*kpgnorm(ipw))**2,sqrt(cgnk2(ipw))
    !   end do
    !   close(96)
    !   ABI_FREE(cgnk)
    !   ABI_FREE(cgnk2)
    ! end do
  end subroutine prt_cg_from_wf
end module m_hightemp
!!***
