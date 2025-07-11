
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_phdos

  use defs_basis
!FB  use m_nctk
  use m_errors
  use m_abicore
  use m_phonons
  use m_xmpi
  use m_io_tools
  use m_ifc,              only : ifc_type
  use m_crystal,          only : crystal_t
  use m_ddb,              only : ddb_type
  use m_tdep_phi2,        only : Eigen_type, tdep_write_yaml, tdep_write_dij, tdep_calc_dij,&
&                                tdep_init_phi2, tdep_destroy_phi2, Phi2_type
  use m_tdep_qpt,         only : Qpoints_type
  use m_tdep_readwrite,   only : Input_type, MPI_enreg_type
  use m_tdep_latt,        only : Lattice_type
  use m_tdep_sym,         only : Symetries_type
  use m_tdep_shell,       only : Shell_type
  use m_tdep_abitypes,    only : Qbz_type, tdep_ifc2phi2, tdep_read_ifc, tdep_write_ifc, &
&                                tdep_write_ddb,tdep_init_ifc

  implicit none

  public :: tdep_calc_phdos
  public :: tdep_calc_thermo
  public :: tdep_calc_elastic

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine tdep_calc_phdos(Crystal,DDB,Eigen2nd_MP,Eigen2nd_path,Ifc,Invar,Lattice,MPIdata,natom,&
&                          natom_unitcell,Phi2,PHdos,Qbz,Qpt,Rlatt4abi,Shell2at,Sym)

  implicit none

  integer, intent(in) :: natom,natom_unitcell
  double precision, intent(in) :: Rlatt4abi(3,natom_unitcell,natom)
  type(Input_type),intent(in) :: Invar
  type(phdos_t),intent(out) :: PHdos
  type(Phi2_type),intent(in) :: Phi2
  type(ifc_type),intent(inout) :: Ifc
  type(Lattice_type),intent(in) :: Lattice
  type(Symetries_type),intent(in) :: Sym
  type(crystal_t),intent(inout) :: Crystal
  type(Qbz_type),intent(in) :: Qbz
  type(Qpoints_type),intent(in) :: Qpt
  type(ddb_type),intent(inout) :: DDB
  type(MPI_enreg_type), intent(in) :: MPIdata
  type(Shell_type),intent(in) :: Shell2at
  type(Eigen_type),intent(inout) :: Eigen2nd_path
  type(Eigen_type),intent(inout) :: Eigen2nd_MP

  integer :: prtdos,iqpt,iq_ibz,iomega,iatom
  integer :: dos_ngqpt(3)
  integer :: count_wminmax(2)
  character (len=fnlen):: phdos_fname
  double precision :: dossmear,integ,domega
  double precision :: dos_qshift(3)
  double precision, allocatable :: displ(:,:)
  character(len=500) :: message
  real(dp) :: wminmax(2)
  type(ifc_type) :: Ifc_tmp
  type(Phi2_type) :: Phi2_tmp

!FB  integer :: jatom,ii,jj
!FB  double precision, allocatable :: omega (:)
!FB  double complex  , allocatable :: dij   (:,:)
!FB  double complex  , allocatable :: eigenV(:,:)

  write(Invar%stdout,'(a)')' '
  write(Invar%stdout,'(a)') ' #############################################################################'
  write(Invar%stdout,'(a)') ' ################### vibrational Density OF States (vDOS) ####################'
  write(Invar%stdout,'(a)') ' #############################################################################'
  write(Invar%stdout,'(a)') ' See the vdos.dat and TDEP_PHDOS* files'

! Copy Phi2_tmp to Ifc%atmfrc
! ===========================
  call tdep_init_phi2(Phi2_tmp,Invar%loto,natom)
  Phi2_tmp%SR =Phi2%SR
  if (Invar%loto) then
    Phi2_tmp%Tot=Phi2%Tot
    Phi2_tmp%LR =Phi2%LR
  end if
  call tdep_ifc2phi2(Ifc%dipdip,Ifc,Invar,Lattice,natom_unitcell,0,Phi2_tmp,Rlatt4abi,Shell2at,Sym)

! Write Ifc%atmfrc in the ifc_out.dat file
! ========================================
  if (MPIdata%iam_master) call tdep_write_ifc(Crystal,Ifc,Invar,natom_unitcell,0)

! For test purpose : read the previous IFC from ifc_out.dat and write it in ifc_check.dat
! =======================================================================================
  if (Invar%readifc.eq.2) then
    call tdep_init_ifc(Crystal,DDB,Ifc_tmp,Invar,Lattice,MPIdata,Phi2_tmp,Rlatt4Abi,Shell2at,Sym)
    if (MPIdata%iam_master) then
!     Read IFC from ifc_out.dat (readifc=2)
      call tdep_read_ifc(Ifc_tmp,Invar,natom_unitcell)
!     Copy Ifc_tmp%atmfrc to Phi2_tmp
      call tdep_ifc2phi2(Ifc_tmp%dipdip,Ifc_tmp,Invar,Lattice,natom_unitcell,1,Phi2_tmp,Rlatt4abi,Shell2at,Sym)
!     Copy Phi2_tmp to Ifc_tmp%atmfrc
      call tdep_ifc2phi2(Ifc_tmp%dipdip,Ifc_tmp,Invar,Lattice,natom_unitcell,0,Phi2_tmp,Rlatt4abi,Shell2at,Sym)
!     Write IFC in ifc_check.dat (for check)
      call tdep_write_ifc(Crystal,Ifc_tmp,Invar,natom_unitcell,1)

!     Write the Phi2-tmp.dat file
      if (Invar%debug) then
        write(Invar%stdout,'(a)') ' See the Phi2-tmp.dat file corresponding to the ifc_out.dat/Phi2 file'
        open(unit=55,file=trim(Invar%output_prefix)//'_Phi2-tmp.dat')
        do iatom=1,3*natom
          write(55,'(10000(f10.6,1x))') Phi2_tmp%SR(iatom,:)
        end do
        close(55)
      end if
    end if
    call Ifc_tmp%free()
  end if
  call tdep_destroy_phi2(Phi2_tmp,Invar%loto)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ON THE FINE GRID !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the DOS
! ===============
  prtdos=1 !Gaussian
!  prtdos=2 !Tetra
  dossmear=4.5d-6
!  dossmear=4.5d-5
  dos_qshift(:)=0.5d0
  dos_ngqpt(:)=Invar%ngqpt2(:)
  write(Invar%stdout,'(a)') ' Compute the vDOS'
  ! Only 1 shift in q-mesh
  wminmax = zero
  do
    call PHdos%init(Crystal,Ifc,prtdos,Invar%dosdeltae,dossmear,dos_ngqpt,1,dos_qshift, &
      "freq_displ", wminmax, count_wminmax, XMPI_WORLD,prtout=.false.)
     if (all(count_wminmax == 0)) exit
     wminmax(1) = wminmax(1) - abs(wminmax(1)) * 0.05
     wminmax(2) = wminmax(2) + abs(wminmax(2)) * 0.05
     call phdos%free()
     write(message, "(a, 2f8.5)")"Initial frequency mesh not large enough. Recomputing PHDOS with wmin, wmax: ",wminmax
     call wrtout(std_out, message)
  end do
  write(Invar%stdout,'(a)') ' ------- achieved'
  write(Invar%stdout,'(a)') ' (Please, pay attention to convergency wrt the BZ mesh : the ngqpt2 input variable)'

! Print the DOS
! =============
  phdos_fname = trim(Invar%output_prefix)//"_PHDOS"
  if (MPIdata%iam_master) call phdos%print(phdos_fname)
  domega=(Invar%dosdeltae*Ha_meV)
  integ=0.d0
  do iomega=1,PHdos%nomega
    integ=integ + domega*PHdos%phdos(iomega)
  end do
  PHdos%phdos(:)=PHdos%phdos(:)/integ
  if (MPIdata%iam_master) then
    open(unit=56,file=trim(Invar%output_prefix)//'_vdos.dat')
    do iomega=1,PHdos%nomega
      if (Invar%enunit.eq.0) write(56,'(2(f18.6,1x))') PHdos%omega(iomega)*Ha_eV*1000,PHdos%phdos(iomega)
      if (Invar%enunit.eq.1) write(56,'(2(f18.6,1x))') PHdos%omega(iomega)*Ha_cmm1   ,PHdos%phdos(iomega)
      if (Invar%enunit.eq.2) write(56,'(2(f18.6,1x))') PHdos%omega(iomega)*1000      ,PHdos%phdos(iomega)
      if (Invar%enunit.eq.3) write(56,'(2(f18.6,1x))') PHdos%omega(iomega)*Ha_THz    ,PHdos%phdos(iomega)
    end do
    close(56)
  end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ON THE PATH GRID !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute the spectrum and the dynamical matrix
! =============================================
! Compute the frequencies (1)
! =======================
!FB  Eigen2nd_path%eigenval=zero ; Eigen2nd_path%eigenvec=zero ; Eigen2nd_path%dynmat=zero
!FB  ABI_MALLOC(dij   ,(3*Invar%natom_unitcell,3*Invar%natom_unitcell))
!FB  ABI_MALLOC(eigenV,(3*Invar%natom_unitcell,3*Invar%natom_unitcell))
!FB  ABI_MALLOC(omega,(3*Invar%natom_unitcell))
!FB  do iqpt=1,Qpt%nqpt
!FB    omega=zero ; eigenV=zero ; dij=zero
!FB    call tdep_calc_dij(dij,eigenV,iqpt,Invar,omega,Phi2%Tot,Qpt%qpt_cart(:,iqpt),Rlatt_cart)
!FB    Eigen2nd_path%eigenval(:,iqpt)= omega(:)
!FB    do iatom=1,Invar%natom_unitcell
!FB      do ii=1,3
!FB        do jatom=1,Invar%natom_unitcell
!FB          do jj=1,3
!FB            Eigen2nd_path%eigenvec(1,ii,iatom,jj,jatom,iqpt)= real(eigenV(ii+3*(iatom-1),jj+3*(jatom-1)))
!FB            Eigen2nd_path%eigenvec(2,ii,iatom,jj,jatom,iqpt)=aimag(eigenV(ii+3*(iatom-1),jj+3*(jatom-1)))
!FB            Eigen2nd_path%dynmat(1,ii,iatom,jj,jatom,iqpt)  = real(dij(ii+3*(iatom-1),jj+3*(jatom-1)))
!FB            Eigen2nd_path%dynmat(2,ii,iatom,jj,jatom,iqpt)  =aimag(dij(ii+3*(iatom-1),jj+3*(jatom-1)))
!FB          end do
!FB        end do
!FB      end do
!FB    end do
!FB  end do
!FB  ABI_FREE(dij)
!FB  ABI_FREE(eigenV)
!FB  ABI_FREE(omega)
!FB! Write the Dij, eigenvalues and eigenvectors (in ASCII)
!FB! ======================================================
!FB  if (MPIdata%iam_master) then
!FB    open(unit=51,file=trim(Invar%output_prefix)//'eigenvectors-path-1.dat')
!FB    open(unit=52,file=trim(Invar%output_prefix)//'dij-path-1.dat')
!FB    open(unit=53,file=trim(Invar%output_prefix)//'omega-path-1.dat')
!FB    if (Invar%enunit.eq.0) write(53,'(a)') '# Phonon frequencies in meV'
!FB    if (Invar%enunit.eq.1) write(53,'(a)') '# Phonon frequencies in cm-1'
!FB    if (Invar%enunit.eq.2) write(53,'(a)') '# Phonon frequencies in mHa'
!FB    if (Invar%enunit.eq.3) write(53,'(a)') '# Phonon frequencies in THz'
!FB    do iqpt=1,Qpt%nqpt
!FB      call tdep_write_dij(Eigen2nd_path,iqpt,Invar,Qpt%qpt_red(:,iqpt))
!FB    end do
!FB    close(51)
!FB    close(52)
!FB    close(53)
!FB  end if

! Compute the frequencies (2)
! =======================
  ABI_MALLOC(displ,(2*3*natom_unitcell*3*natom_unitcell,Qpt%nqpt)); displ(:,:)=zero
  Eigen2nd_path%eigenval=zero ; Eigen2nd_path%eigenvec=zero ; Eigen2nd_path%dynmat=zero
  do iqpt=1,Qpt%nqpt
    call ifc%fourq(Crystal,Qpt%qpt_red(:,iqpt),Eigen2nd_path%eigenval(:,iqpt),displ(:,iqpt),&
&                  out_eigvec=Eigen2nd_path%eigenvec(:,:,:,:,:,iqpt),&
&                  out_d2cart=Eigen2nd_path%dynmat  (:,:,:,:,:,iqpt))
  end do
  ABI_FREE(displ)
! Write the Dij, eigenvalues and eigenvectors (in ASCII and YAML)
! ======================================================
  if (MPIdata%iam_master) then
!FB    open(unit=51,file=trim(Invar%output_prefix)//'eigenvectors-path-2.dat')
!FB    open(unit=52,file=trim(Invar%output_prefix)//'dij-path-2.dat')
!FB    open(unit=53,file=trim(Invar%output_prefix)//'omega-path-2.dat')
    open(unit=51,file=trim(Invar%output_prefix)//'_eigenvectors.dat')
    open(unit=52,file=trim(Invar%output_prefix)//'_dij.dat')
    open(unit=53,file=trim(Invar%output_prefix)//'_omega.dat')
    if (Invar%enunit.eq.0) write(53,'(a)') '# Phonon frequencies in meV'
    if (Invar%enunit.eq.1) write(53,'(a)') '# Phonon frequencies in cm-1'
    if (Invar%enunit.eq.2) write(53,'(a)') '# Phonon frequencies in mHa'
    if (Invar%enunit.eq.3) write(53,'(a)') '# Phonon frequencies in THz'
    do iqpt=1,Qpt%nqpt
      call tdep_write_dij(Eigen2nd_path,iqpt,Invar,Qpt%qpt_red(:,iqpt))
    end do
    close(51)
    close(52)
    close(53)
    call tdep_write_yaml(Eigen2nd_path,Qpt,Invar%output_prefix)
  end if

!FB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!FB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ON THE MP GRID !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!FB!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!FB! Compute the spectrum and the dynamical matrix
!FB! =============================================
!FB  if (Qbz%nqbz.ne.Qbz%nqibz) then
!FB    write(Invar%stdlog,*) 'DIFFERENCES BETWEEN NQBZ \& NQIBZ :',Qbz%nqbz,Qbz%nqibz
!FB  end if
!FB! Compute the frequencies(1)
!FB! =======================
!FB  Eigen2nd_MP%eigenval=zero ; Eigen2nd_MP%eigenvec=zero ; Eigen2nd_MP%dynmat=zero
!FB  ABI_MALLOC(dij   ,(3*Invar%natom_unitcell,3*Invar%natom_unitcell))
!FB  ABI_MALLOC(eigenV,(3*Invar%natom_unitcell,3*Invar%natom_unitcell))
!FB  ABI_MALLOC(omega,(3*Invar%natom_unitcell))
!FB  do iq_ibz=1,Qbz%nqibz
!FB    omega=zero ; eigenV=zero ; dij=zero
!FB    call tdep_calc_dij(dij,eigenV,iq_ibz,Invar,omega,Phi2,Qbz%qibz_cart(:,iq_ibz),Rlatt_cart)
!FB    Eigen2nd_MP%eigenval(:,iq_ibz)= omega(:)
!FB    do iatom=1,Invar%natom_unitcell
!FB      do ii=1,3
!FB        do jatom=1,Invar%natom_unitcell
!FB          do jj=1,3
!FB            Eigen2nd_MP%eigenvec(1,ii,iatom,jj,jatom,iq_ibz)= real(eigenV(ii+3*(iatom-1),jj+3*(jatom-1)))
!FB            Eigen2nd_MP%eigenvec(2,ii,iatom,jj,jatom,iq_ibz)=aimag(eigenV(ii+3*(iatom-1),jj+3*(jatom-1)))
!FB            Eigen2nd_MP%dynmat(1,ii,iatom,jj,jatom,iq_ibz)  = real(dij(ii+3*(iatom-1),jj+3*(jatom-1)))
!FB            Eigen2nd_MP%dynmat(2,ii,iatom,jj,jatom,iq_ibz)  =aimag(dij(ii+3*(iatom-1),jj+3*(jatom-1)))
!FB          end do
!FB        end do
!FB      end do
!FB    end do
!FB  end do
!FB  ABI_FREE(dij)
!FB  ABI_FREE(eigenV)
!FB  ABI_FREE(omega)
!FB  if (MPIdata%iam_master) then
!FB    open(unit=51,file=trim(Invar%output_prefix)//'eigenvectors-MP-1.dat')
!FB    open(unit=52,file=trim(Invar%output_prefix)//'dij-MP-1.dat')
!FB    open(unit=53,file=trim(Invar%output_prefix)//'omega-MP-1.dat')
!FB    if (Invar%enunit.eq.0) write(53,'(a)') '# Phonon frequencies in meV'
!FB    if (Invar%enunit.eq.1) write(53,'(a)') '# Phonon frequencies in cm-1'
!FB    if (Invar%enunit.eq.2) write(53,'(a)') '# Phonon frequencies in mHa'
!FB    if (Invar%enunit.eq.3) write(53,'(a)') '# Phonon frequencies in THz'
!FB    do iq_ibz=1,Qbz%nqibz
!FB      call tdep_write_dij(Eigen2nd_MP,iq_ibz,Invar,Qbz%qibz(:,iq_ibz))
!FB    end do
!FB    close(51)
!FB    close(52)
!FB    close(53)
!FB  end if
!FB
! Compute the frequencies (2)
! =======================
  ABI_MALLOC(displ,(2*3*natom_unitcell*3*natom_unitcell,Qbz%nqibz)); displ(:,:)=zero
  Eigen2nd_MP%eigenval=zero ; Eigen2nd_MP%eigenvec=zero ; Eigen2nd_MP%dynmat=zero
  do iq_ibz=1,Qbz%nqibz
    call ifc%fourq(Crystal,Qbz%qibz(:,iq_ibz),Eigen2nd_MP%eigenval(:,iq_ibz),displ(:,iq_ibz),&
&                  out_eigvec=Eigen2nd_MP%eigenvec(:,:,:,:,:,iq_ibz),&
&                  out_d2cart=Eigen2nd_MP%dynmat  (:,:,:,:,:,iq_ibz))
  end do
  ABI_FREE(displ)
!FB  if (MPIdata%iam_master) then
!FB    open(unit=51,file=trim(Invar%output_prefix)//'eigenvectors-MP-2.dat')
!FB    open(unit=52,file=trim(Invar%output_prefix)//'dij-MP-2.dat')
!FB    open(unit=53,file=trim(Invar%output_prefix)//'omega-MP-2.dat')
!FB    if (Invar%enunit.eq.0) write(53,'(a)') '# Phonon frequencies in meV'
!FB    if (Invar%enunit.eq.1) write(53,'(a)') '# Phonon frequencies in cm-1'
!FB    if (Invar%enunit.eq.2) write(53,'(a)') '# Phonon frequencies in mHa'
!FB    if (Invar%enunit.eq.3) write(53,'(a)') '# Phonon frequencies in THz'
!FB    do iq_ibz=1,Qbz%nqibz
!FB      call tdep_write_dij(Eigen2nd_MP,iq_ibz,Invar,Qbz%qibz(:,iq_ibz))
!FB    end do
!FB    close(51)
!FB    close(52)
!FB    close(53)
!FB  end if

! Write the DDB
! =============
!BeginFB
! if (MPIdata%iam_master) call tdep_write_ddb(Crystal,DDB,Eigen2nd_MP,Invar,Lattice,MPIdata,Qbz,Sym)
!EndFB

end subroutine tdep_calc_phdos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tdep_calc_thermo(Invar,Lattice,MPIdata,PHdos,U0)

  implicit none

  double precision, intent(in) :: U0
  type(Input_type),intent(in) :: Invar
  type(Lattice_type), intent(inout) :: Lattice
  type(MPI_enreg_type), intent(in) :: MPIdata
  type(phdos_t),intent(in) :: PHdos

  integer :: iomega,itemp,iatom,itypat
  double precision :: k_B,wovert,heatcapa,entropy,internalE,freeE,expm2x,ln2shx,cothx,xx
  double precision :: Ftot,domega,MSD,Omega_m2,mass_amu,vdos

  write(Invar%stdout,'(a)')' '
  write(Invar%stdout,'(a)') ' #############################################################################'
  write(Invar%stdout,'(a)') ' ################# Thermodynamic quantities: Free energy,...##################'
  write(Invar%stdout,'(a)') ' #############################################################################'
  write(Invar%stdout,'(a)') ' See the thermo.dat file'

! The heat capacity, entropy, internal and free energies (direct calculation)
! ===========================================================================
  mass_amu=zero
  do iatom=1,Invar%natom_unitcell
    itypat=Invar%typat_unitcell(iatom)
    mass_amu=mass_amu+Invar%amu(itypat)
  end do
  mass_amu=mass_amu*amu_emass/real(Invar%natom_unitcell)

!FB  k_B=8.617343d-5 !in eV/K
  k_B=kb_HaK*Ha_eV
  domega=(Invar%dosdeltae*Ha_meV)
  wovert=1.d0/(2*Invar%temperature*k_B)
  heatcapa=0.d0 ; entropy=0.d0 ; internalE=0.d0 ; freeE=0.d0 ; MSD=0.d0 ; Omega_m2=0.d0 ; vdos=0.d0
  do iomega=1,PHdos%nomega
    xx=PHdos%omega(iomega)*Ha_eV
    if (xx.lt.tol8) cycle
    expm2x=exp(-2.d0*wovert*xx)
    ln2shx=wovert*xx+log(1.d0-expm2x)
    cothx=(1.d0+expm2x)/(1.d0-expm2x)
    heatcapa =heatcapa  + (wovert*xx/sinh(wovert*xx))**2*PHdos%phdos(iomega)*domega
    internalE=internalE + (wovert*xx)*cothx*PHdos%phdos(iomega)*domega
    entropy  =entropy   + ((wovert*xx)*cothx-ln2shx)*PHdos%phdos(iomega)*domega
    freeE    =freeE     + log(2*sinh(wovert*xx))*PHdos%phdos(iomega)*domega
    MSD      =MSD       + (cothx/PHdos%omega(iomega))*PHdos%phdos(iomega)*domega
    Omega_m2 =Omega_m2  + (1.d0/PHdos%omega(iomega))**2*PHdos%phdos(iomega)*domega
    vdos     =vdos      + PHdos%phdos(iomega)*domega
  end do
  heatcapa=heatcapa*3
  entropy=entropy*3
  internalE=internalE*3*k_B*Invar%temperature
  freeE=freeE*3*k_B*Invar%temperature
  MSD=MSD*3.d0/mass_amu/2.d0
  Omega_m2=Omega_m2*3.d0
  Ftot=U0*Ha_eV+freeE
  Lattice%HeatCapa_V=heatcapa
  if (.not.MPIdata%iam_master) return
! End of the calculation --> RETURN
! =================================

  open(unit=20,file=trim(Invar%output_prefix)//'_thermo.dat')
  write(20,'(a)')'============= Direct results (without any inter/extrapolation) =================='
  write(20,'(1x,a,f10.3)')'For present temperature (in Kelvin): T= ',Invar%temperature
  write(20,'(1x,a,f12.3)')'  The cold contribution (in eV/atom): U_0 =',U0*Ha_eV
  write(20,'(1x,a,f10.3)')'  The specific heat (in k_b/atom): C_v=',heatcapa
  write(20,'(1x,a,f10.3)')'  The vibrational entropy (in k_b/atom): S_vib =',entropy
  write(20,'(1x,a,f10.3)')'  The internal energy (in eV/atom): U_vib =',internalE
  write(20,'(1x,a,f10.3)')'  The vibrational contribution (in eV/atom): F_vib = U_vib -T.S_vib =',freeE
  write(20,'(1x,a,f12.3)')'  The harmonic free energy (in eV/atom) -->  F_tot^HA = U_0 + F_vib =',Ftot
  write(20,'(1x,a)')'  Useful quantities for melting :'
  write(20,'(1x,a,f10.3)')'     The mean square displacement (in a.u.): sqrt(<u^2>) =',(MSD)**0.5
  write(20,'(1x,a,f10.3)')'     The <Omega^(-2)> factor (in THz^(-2)) =',Omega_m2/(Ha_THz)**2
  write(20,'(1x,a,f10.3)')'     The Wigner-Seitz radius (in a.u.) : d_at =',(6*Lattice%ucvol/pi/real(Invar%natom_unitcell))**(1./3.)
  write(20,'(1x,a,f10.3)')'     The average mass / proton-electron mass ratio (in a.u.) =', mass_amu/amu_emass
  write(20,'(1x,a,f10.3)')'     The Lindemann constant : sqrt(<u^2>)/d_at =',(MSD)**0.5/&
&                               (6*Lattice%ucvol/pi/real(Invar%natom_unitcell))**(1./3.)
  write(20,'(1x,a,f10.3)')'     The integral of vDOS =',vdos
  write(20,'(a)')' '

! The free energy (extrapolation)
! ===============================
  write(20,'(a)')'============= Harmonic Approximation (HA) =================='
  write(20,'(1x,a)')'  Note that the following results come from an EXTRAPOLATION:'
  write(20,'(1x,a,i5,a)')'    1/ F_vib^HA(T) is computed for each T using vDOS(T=',int(Invar%temperature),')'
  write(20,'(1x,a)')'    2/ F_tot^HA(T) = F_vib^HA(T) + U_0'
  write(20,'(a)')'   T         F_vib^HA(T)   F_tot^HA(T)           C_v(T)  '&
&   //'       S_vib(T)        U_vib(T)        MSD(T)'
  do itemp=1,100
    wovert=1.d0/(2*real(itemp)*100*k_B)
    freeE=0.d0
    heatcapa=0.d0
    entropy=0.d0
    internalE=0.d0
    MSD=0.d0
    do iomega=1,PHdos%nomega
      xx=PHdos%omega(iomega)*Ha_eV
      if (xx.lt.tol8) cycle
      expm2x=exp(-2.d0*wovert*xx)
      ln2shx=wovert*xx+log(1.d0-expm2x)
      cothx=(1.d0+expm2x)/(1.d0-expm2x)
      heatcapa =heatcapa  + (wovert*xx/sinh(wovert*xx))**2*PHdos%phdos(iomega)*domega
      internalE=internalE + (wovert*xx)*cothx*PHdos%phdos(iomega)*domega
      entropy  =entropy   + ((wovert*xx)*cothx-ln2shx)*PHdos%phdos(iomega)*domega
      freeE    =freeE     + log(2*sinh(wovert*xx))*PHdos%phdos(iomega)*domega
      MSD      =MSD       + (cothx/PHdos%omega(iomega))*PHdos%phdos(iomega)*domega
    end do
    heatcapa=heatcapa*3
    entropy=entropy*3
    internalE=internalE*3*k_B*itemp*100
    freeE=freeE*3*k_B*itemp*100
    MSD=MSD*3.d0/mass_amu/2.d0
    Ftot=U0*Ha_eV+freeE
    write(20,'(1x,i5,6(1x,f15.3))') itemp*100,freeE,Ftot,heatcapa,entropy,internalE,(MSD)**0.5
  end do
  close(20)

end subroutine tdep_calc_thermo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine tdep_calc_elastic(Phi2,distance,Invar,Lattice)

  integer :: iatom,ii,jj,kk,ll,iatcell,itypat
  integer :: INFO,LWORK
  double precision :: BH,BR,BV,GR,GV,GH,Eaverage,Nuaverage,Laverage,Vp,Vs,Vphi
  double precision :: rho,E1,E2,E3,Nu12,Nu13,Nu23,Nu21,Nu31,Nu32,G23,G13,G12
  double precision :: mass_amu,bohr,V_Debye,T_Debye,A_U,A_B,A_G
  integer, allocatable :: IPIV(:)
  double precision, allocatable :: eigenvalues(:)
  double precision, allocatable :: WORK(:)
  double precision, allocatable :: Sij(:,:),Cij(:,:),aijkl(:,:,:,:),cijkl(:,:,:,:)
  type(Input_type), intent(in) :: Invar
  type(Lattice_type), intent(inout) :: Lattice
  double precision, intent(in) :: distance(Invar%natom,Invar%natom,4)
  double precision, intent(in) :: Phi2(3*Invar%natom,3*Invar%natom)

  write(Invar%stdout,'(a)')' '
  write(Invar%stdout,'(a)') ' #############################################################################'
  write(Invar%stdout,'(a)') ' ######################### Elastic constants #################################'
  write(Invar%stdout,'(a)') ' ################ Bulk and Shear modulus--Sound velocities ###################'
  write(Invar%stdout,'(a)') ' #############################################################################'

  bohr=0.5291772108e-10
! Define atomic mass average
  mass_amu=zero
  do iatom=1,Invar%natom_unitcell
    itypat=Invar%typat_unitcell(iatom)
    mass_amu=mass_amu+Invar%amu(itypat)
  end do
  mass_amu=mass_amu/real(Invar%natom_unitcell)

  rho=(mass_amu/1e3)*Invar%natom_unitcell/Lattice%ucvol/bohr**3/6.022e23

!==========================================================================================
!===================== Elastic constants ==================================================
!==========================================================================================
! New calculation of elastic constants using the formula (12.28 and 12.29 of
! Wallace, Statistical physics of crystals and liquids, Worl Scientific)
  ABI_MALLOC(aijkl,(3,3,3,3)); aijkl(:,:,:,:)=0.d0
  ABI_MALLOC(cijkl,(3,3,3,3)); cijkl(:,:,:,:)=0.d0
  do ii=1,3
    do jj=1,3
      do kk=1,3
        do ll=1,3
          do iatom=1,Invar%natom
            do iatcell=1,Invar%natom_unitcell
              aijkl(ii,jj,kk,ll)=aijkl(ii,jj,kk,ll)-Phi2(ii+(iatcell-1)*3,3*(iatom-1)+jj)&
&               *distance(iatcell,iatom,kk+1)*distance(iatcell,iatom,ll+1)/2.d0/Lattice%ucvol
            end do
          end do
        enddo
      enddo
    enddo
  enddo

  do ii=1,3
    do jj=1,3
      do kk=1,3
        do ll=1,3
          cijkl(ii,jj,kk,ll)=aijkl(ii,kk,jj,ll)+aijkl(jj,kk,ii,ll)-aijkl(ii,jj,kk,ll)
!FB          cijkl(ii,jj,kk,ll)=aijkl(ii,kk,jj,ll)+aijkl(jj,kk,ii,ll)-aijkl(jj,ii,kk,ll)
        enddo
      enddo
    enddo
  enddo
  ABI_FREE(aijkl)

  cijkl(:,:,:,:)=cijkl(:,:,:,:)*29421.033d0

  ABI_MALLOC(Cij,(6,6)) ; Cij(:,:)=0.d0
  Cij(1,1)=cijkl(1,1,1,1) ; Cij(1,2)=cijkl(1,1,2,2) ; Cij(1,3)=cijkl(1,1,3,3)
  Cij(1,4)=cijkl(1,1,2,3) ; Cij(1,5)=cijkl(1,1,1,3) ; Cij(1,6)=cijkl(1,1,1,2)
  Cij(2,1)=cijkl(2,2,1,1) ; Cij(2,2)=cijkl(2,2,2,2) ; Cij(2,3)=cijkl(2,2,3,3)
  Cij(2,4)=cijkl(2,2,2,3) ; Cij(2,5)=cijkl(2,2,1,3) ; Cij(2,6)=cijkl(2,2,1,2)
  Cij(3,1)=cijkl(3,3,1,1) ; Cij(3,2)=cijkl(3,3,2,2) ; Cij(3,3)=cijkl(3,3,3,3)
  Cij(3,4)=cijkl(3,3,2,3) ; Cij(3,5)=cijkl(3,3,1,3) ; Cij(3,6)=cijkl(3,3,1,2)
  Cij(4,1)=cijkl(2,3,1,1) ; Cij(4,2)=cijkl(2,3,2,2) ; Cij(4,3)=cijkl(2,3,3,3)
  Cij(4,4)=cijkl(2,3,2,3) ; Cij(4,5)=cijkl(2,3,1,3) ; Cij(4,6)=cijkl(2,3,1,2)
  Cij(5,1)=cijkl(1,3,1,1) ; Cij(5,2)=cijkl(1,3,2,2) ; Cij(5,3)=cijkl(1,3,3,3)
  Cij(5,4)=cijkl(1,3,2,3) ; Cij(5,5)=cijkl(1,3,1,3) ; Cij(5,6)=cijkl(1,3,1,2)
  Cij(6,1)=cijkl(1,2,1,1) ; Cij(6,2)=cijkl(1,2,2,2) ; Cij(6,3)=cijkl(1,2,3,3)
  Cij(6,4)=cijkl(1,2,2,3) ; Cij(6,5)=cijkl(1,2,1,3) ; Cij(6,6)=cijkl(1,2,1,2)
! Remove the rounding errors before writing (for non regression testing purposes)
  do ii=1,6
    do jj=1,6
      if (abs(Cij(ii,jj)).lt.tol8) Cij(ii,jj)=zero
    end do
  end do
  write(Invar%stdout,'(a)') ' '
  write(Invar%stdout,'(a)') ' ========== Using the formulation proposed by Wallace (using the IFC) ========='
  write(Invar%stdout,'(a)') ' Cijkl [in GPa]='
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | C11 C12 C13 C14 C15 C16 |   ',Cij(1,1),Cij(1,2),Cij(1,3),Cij(1,4),Cij(1,5),Cij(1,6)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | C21 C22 C23 C24 C25 C26 |   ',Cij(2,1),Cij(2,2),Cij(2,3),Cij(2,4),Cij(2,5),Cij(2,6)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | C31 C32 C33 C34 C35 C36 |   ',Cij(3,1),Cij(3,2),Cij(3,3),Cij(3,4),Cij(3,5),Cij(3,6)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | C41 C42 C43 C44 C45 C46 | = ',Cij(4,1),Cij(4,2),Cij(4,3),Cij(4,4),Cij(4,5),Cij(4,6)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | C51 C52 C53 C54 C55 C56 |   ',Cij(5,1),Cij(5,2),Cij(5,3),Cij(5,4),Cij(5,5),Cij(5,6)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | C61 C62 C63 C64 C65 C66 |   ',Cij(6,1),Cij(6,2),Cij(6,3),Cij(6,4),Cij(6,5),Cij(6,6)

! Compute the eigenvalues of the Cij matrix in order to find the Born-Huang
! stability criterion (see Wallace, Thermodynamics of crystal, p39)
  ABI_MALLOC(WORK,(1))
  ABI_MALLOC(Sij,(6,6)) ; Sij(:,:)=0.d0
  ABI_MALLOC(eigenvalues,(6)) ; eigenvalues(:)=0.d0
  Sij(:,:)=Cij(:,:)
  LWORK=-1
  call DSYEV('N','U',6,Sij,6,eigenvalues,WORK,LWORK,INFO)
  LWORK=WORK(1)
  ABI_FREE(WORK)
  ABI_MALLOC(WORK,(LWORK))
  call DSYEV('N','U',6,Sij,6,eigenvalues,WORK,LWORK,INFO)
  do ii=1,6
    if (eigenvalues(ii).lt.0.d0) then
      write(Invar%stdout,'(a)') ' WARNING :'
      write(Invar%stdout,'(a,i3,a,1x,f8.3,1x)') 'The eigenvalue number',ii,'is negative and equals to',eigenvalues(ii)
      write(Invar%stdout,'(a)') 'The Born-Huang stability criterion is not fulfilled'
    end if
  end do

  ABI_FREE(WORK)
  ABI_FREE(Sij)
  ABI_FREE(eigenvalues)

! For an anisotropic material
  write(Invar%stdout,'(a)') ' '
  write(Invar%stdout,'(a)') ' ========== For an Anisotropic Material ======================================='
  ABI_MALLOC(Sij,(6,6)) ; Sij(:,:)=0.d0
  ABI_MALLOC(IPIV,(6)); IPIV(:)=0
  ABI_MALLOC(WORK,(6)); WORK(:)=0.d0
  Sij(:,:)=Cij(:,:)
  call DGETRF(6,6,Sij,6,IPIV,INFO)
  call DGETRI(6,Sij,6,IPIV,WORK,6,INFO)
  ABI_FREE(IPIV)
  ABI_FREE(WORK)
  write(Invar%stdout,'(a)') ' Sijkl [in GPa-1]='
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | S11 S12 S13 S14 S15 S16 |   ',Sij(1,1),Sij(1,2),Sij(1,3),Sij(1,4),Sij(1,5),Sij(1,6)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | S21 S22 S23 S24 S25 S26 |   ',Sij(2,1),Sij(2,2),Sij(2,3),Sij(2,4),Sij(2,5),Sij(2,6)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | S31 S32 S33 S34 S35 S36 |   ',Sij(3,1),Sij(3,2),Sij(3,3),Sij(3,4),Sij(3,5),Sij(3,6)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | S41 S42 S43 S44 S45 S46 | = ',Sij(4,1),Sij(4,2),Sij(4,3),Sij(4,4),Sij(4,5),Sij(4,6)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | S51 S52 S53 S54 S55 S56 |   ',Sij(5,1),Sij(5,2),Sij(5,3),Sij(5,4),Sij(5,5),Sij(5,6)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | S61 S62 S63 S64 S65 S66 |   ',Sij(6,1),Sij(6,2),Sij(6,3),Sij(6,4),Sij(6,5),Sij(6,6)
  Lattice%Sij=Sij
  ABI_FREE(Sij)

! For an orthotropic material
  write(Invar%stdout,'(a)') ' '
  write(Invar%stdout,'(a)') ' ========== For an Orthotropic Material (see B. M. Lempriere (1968)) =========='

! Mean value of the off-diagonal elements
  Cij(1,2)=(Cij(1,2)+Cij(2,1))/2.d0 ; Cij(2,1)=Cij(1,2)
  Cij(1,3)=(Cij(1,3)+Cij(3,1))/2.d0 ; Cij(3,1)=Cij(1,3)
  Cij(2,3)=(Cij(2,3)+Cij(3,2))/2.d0 ; Cij(3,2)=Cij(2,3)

! Young's modulus
  E1=(Cij(1,1)*Cij(2,2)*Cij(3,3)+2.d0*Cij(2,3)*Cij(1,2)*Cij(1,3)-Cij(1,1)*Cij(2,3)**2-Cij(2,2)*Cij(1,3)**2-Cij(3,3)*Cij(1,2)**2)&
&   /(Cij(2,2)*Cij(3,3)-Cij(2,3)**2)
  E2=(Cij(1,1)*Cij(2,2)*Cij(3,3)+2.d0*Cij(2,3)*Cij(1,2)*Cij(1,3)-Cij(1,1)*Cij(2,3)**2-Cij(2,2)*Cij(1,3)**2-Cij(3,3)*Cij(1,2)**2)&
&   /(Cij(1,1)*Cij(3,3)-Cij(1,3)**2)
  E3=(Cij(1,1)*Cij(2,2)*Cij(3,3)+2.d0*Cij(2,3)*Cij(1,2)*Cij(1,3)-Cij(1,1)*Cij(2,3)**2-Cij(2,2)*Cij(1,3)**2-Cij(3,3)*Cij(1,2)**2)&
&   /(Cij(1,1)*Cij(2,2)-Cij(1,2)**2)
  write(Invar%stdout,'(a,3(f8.3,1x))') ' Young modulus E1, E2 and E3 [in GPa]=',E1,E2,E3

! Poisson Ratio
  Nu21=(Cij(1,2)*Cij(3,3)-Cij(1,3)*Cij(2,3))/(Cij(1,1)*Cij(3,3)-Cij(1,3)**2)
  Nu31=(Cij(1,3)*Cij(2,2)-Cij(1,2)*Cij(2,3))/(Cij(1,1)*Cij(2,2)-Cij(1,2)**2)
  Nu23=(Cij(1,1)*Cij(2,3)-Cij(1,2)*Cij(1,3))/(Cij(1,1)*Cij(3,3)-Cij(1,3)**2)
  Nu12=(Cij(1,2)*Cij(3,3)-Cij(1,3)*Cij(2,3))/(Cij(2,2)*Cij(3,3)-Cij(2,3)**2)
  Nu13=(Cij(2,2)*Cij(1,3)-Cij(1,2)*Cij(2,3))/(Cij(2,2)*Cij(3,3)-Cij(2,3)**2)
  Nu32=(Cij(1,1)*Cij(2,3)-Cij(1,2)*Cij(1,3))/(Cij(1,1)*Cij(2,2)-Cij(1,2)**2)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' Poisson ratio Nu21, Nu31, Nu23, Nu12, Nu13 and Nu32=',Nu21,Nu31,Nu23,Nu12,Nu13,Nu32

! Shear modulus
  G23=Cij(4,4) ; G13=Cij(5,5) ; G12=Cij(6,6)
  write(Invar%stdout,'(a,3(f8.3,1x))') ' Shear modulus G23, G13 and G12 [in GPa]=',G23,G13,G12

! Compliance matrix
  ABI_MALLOC(Sij,(6,6)) ; Sij(:,:)=0.d0
  Sij(1,1)= 1.d0/E1 ; Sij(1,2)=-Nu21/E2 ; Sij(1,3)=-Nu31/E3 ; Sij(1,4)=0.d0     ; Sij(1,5)=0.d0     ; Sij(1,6)=0.d0
  Sij(2,1)=-Nu12/E1 ; Sij(2,2)= 1.d0/E2 ; Sij(2,3)=-Nu32/E3 ; Sij(2,4)=0.d0     ; Sij(2,5)=0.d0     ; Sij(2,6)=0.d0
  Sij(3,1)=-Nu13/E1 ; Sij(3,2)=-Nu23/E2 ; Sij(3,3)= 1.d0/E3 ; Sij(3,4)=0.d0     ; Sij(3,5)=0.d0     ; Sij(3,6)=0.d0
  Sij(4,1)= 0.d0    ; Sij(4,2)= 0.d0    ; Sij(4,3)= 0.d0    ; Sij(4,4)=1.d0/G23 ; Sij(4,5)=0.d0     ; Sij(4,6)=0.d0
  Sij(5,1)= 0.d0    ; Sij(5,2)= 0.d0    ; Sij(5,3)= 0.d0    ; Sij(5,4)=0.d0     ; Sij(5,5)=1.d0/G13 ; Sij(5,6)=0.d0
  Sij(6,1)= 0.d0    ; Sij(6,2)= 0.d0    ; Sij(6,3)= 0.d0    ; Sij(6,4)=0.d0     ; Sij(6,5)=0.d0     ; Sij(6,6)=1.d0/G12
! Remove the rounding errors before writing (for non regression testing purposes)
  do ii=1,6
    do jj=1,6
      if (abs(Sij(ii,jj)).lt.tol8) Sij(ii,jj)=zero
    end do
  end do
  write(Invar%stdout,'(a)') ' Sijkl [in GPa-1]='
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | S11 S12 S13 S14 S15 S16 |   ',Sij(1,1),Sij(1,2),Sij(1,3),Sij(1,4),Sij(1,5),Sij(1,6)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | S21 S22 S23 S24 S25 S26 |   ',Sij(2,1),Sij(2,2),Sij(2,3),Sij(2,4),Sij(2,5),Sij(2,6)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | S31 S32 S33 S34 S35 S36 |   ',Sij(3,1),Sij(3,2),Sij(3,3),Sij(3,4),Sij(3,5),Sij(3,6)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | S41 S42 S43 S44 S45 S46 | = ',Sij(4,1),Sij(4,2),Sij(4,3),Sij(4,4),Sij(4,5),Sij(4,6)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | S51 S52 S53 S54 S55 S56 |   ',Sij(5,1),Sij(5,2),Sij(5,3),Sij(5,4),Sij(5,5),Sij(5,6)
  write(Invar%stdout,'(a,6(f8.3,1x))') ' | S61 S62 S63 S64 S65 S66 |   ',Sij(6,1),Sij(6,2),Sij(6,3),Sij(6,4),Sij(6,5),Sij(6,6)

!==========================================================================================
!===================== Bulk and Shear modulus--Sound velocities ===========================
!==========================================================================================
! Voigt notation
  write(Invar%stdout,'(a,f9.3)')' For density rho [in kg.m-3]=',rho
  write(Invar%stdout,'(a)')' '
  write(Invar%stdout,'(a)')' ========================= Voigt average (constant strain) ==================='
  BV=((Cij(1,1)+Cij(2,2)+Cij(3,3))+2.d0*(Cij(1,2)+Cij(1,3)+Cij(2,3)))/9.d0
  GV=((Cij(1,1)+Cij(2,2)+Cij(3,3))-     (Cij(1,2)+Cij(1,3)+Cij(2,3))+3.d0*(Cij(4,4)+Cij(5,5)+Cij(6,6)))/15.d0
  write(Invar%stdout,'(2(a,f9.3))')' ISOTHERMAL modulus [in GPa]: Bulk Kt=',BV,' and Shear G=',GV
  Eaverage=9.d0*BV*GV/(3*BV+GV)
  Nuaverage=0.5*(1.d0-(3.d0*GV)/(3.d0*BV+GV) )
  Laverage=(3.d0*BV-2.d0*GV)/3.d0
  write(Invar%stdout,'(3(a,f9.3))')' Average of Young modulus E [in GPa]=',Eaverage,' Lame modulus Lambda [in GPa]=',Laverage,&
&   ' and Poisson ratio Nu=',Nuaverage
  Vp=dsqrt(1.d9*(BV+4.d0*GV/3.d0)/rho)
  Vs=dsqrt(1.d9*GV/rho)
  Vphi=dsqrt(1.d9*BV/rho)
  write(Invar%stdout,'(3(a,f9.3,1x))')' Velocities [in m.s-1]: compressional Vp=',Vp,' shear Vs=',Vs,' and bulk Vphi=',Vphi
  V_Debye=(1./3.*(1./Vp**3+2./Vs**3))**(-1./3.)
  T_Debye= (V_Debye/1.d3/Bohr_Ang/1.d-13*Time_Sec)*(six*pi**2*Invar%natom_unitcell/Lattice%ucvol)**(1./3.)*Ha_K
  write(Invar%stdout,'(2(a,f9.3,1x))')' Debye velocity [in m.s-1]=',V_Debye,' and temperature [in K]=',T_Debye

! Reuss notation
  write(Invar%stdout,'(a)')' '
  write(Invar%stdout,'(a)')' ========================= Reuss average (constant stress) ==================='
  BR=1.d0/(Sij(1,1)+Sij(2,2)+Sij(3,3)+2.d0*(Sij(1,2)+Sij(1,3)+Sij(2,3)))
  GR=15.d0/(4.d0*(Sij(1,1)+Sij(2,2)+Sij(3,3))-4.d0*(Sij(1,2)+Sij(1,3)+Sij(2,3))+3.d0*(Sij(4,4)+Sij(5,5)+Sij(6,6)))
  write(Invar%stdout,'(2(a,f9.3))')' ISOTHERMAL modulus [in GPa]: Bulk Kt=',BR,' and Shear G=',GR
  Eaverage=9.d0*BR*GR/(3*BR+GR)
  Nuaverage=0.5*(1.d0-(3.d0*GR)/(3.d0*BR+GR) )
  Laverage=(3.d0*BR-2.d0*GR)/3.d0
  write(Invar%stdout,'(3(a,f9.3))')' Average of Young modulus E [in GPa]=',Eaverage,' Lame modulus Lambda [in GPa]=',Laverage,&
&   ' and Poisson ratio Nu=',Nuaverage
  Vp=dsqrt(1.d9*(BR+4.d0*GR/3.d0)/rho)
  Vs=dsqrt(1.d9*GR/rho)
  Vphi=dsqrt(1.d9*BR/rho)
  write(Invar%stdout,'(3(a,f9.3,1x))')' Velocities [in m.s-1]: compressional Vp=',Vp,' shear Vs=',Vs,' and bulk Vphi=',Vphi
  V_Debye=(1./3.*(1./Vp**3+2./Vs**3))**(-1./3.)
  T_Debye= (V_Debye/1.d3/Bohr_Ang/1.d-13*Time_Sec)*(six*pi**2*Invar%natom_unitcell/Lattice%ucvol)**(1./3.)*Ha_K
  write(Invar%stdout,'(2(a,f9.3,1x))')' Debye velocity [in m.s-1]=',V_Debye,' and temperature [in K]=',T_Debye

! Voigt-Reuss-Hill notation
  write(Invar%stdout,'(a)')' '
  write(Invar%stdout,'(a)')' ============================== Hill average ================================='
  BH=(BR+BV)/2.d0
  GH=(GR+GV)/2.d0
  write(Invar%stdout,'(2(a,f9.3))')' ISOTHERMAL modulus [in GPa]: Bulk Kt=',BH,' and Shear G=',GH
  Eaverage=9.d0*BH*GH/(3*BH+GH)
  Nuaverage=0.5*(1.d0-(3.d0*GH)/(3.d0*BH+GH) )
  Laverage=(3.d0*BH-2.d0*GH)/3.d0
  write(Invar%stdout,'(3(a,f9.3))')' Average of Young modulus E [in GPa]=',Eaverage,' Lame modulus Lambda [in GPa]=',Laverage,&
&   ' and Poisson ratio Nu=',Nuaverage
  Vp=dsqrt(1.d9*(BH+4.d0*GH/3.d0)/rho)
  Vs=dsqrt(1.d9*GH/rho)
  Vphi=dsqrt(1.d9*BH/rho)
  write(Invar%stdout,'(3(a,f9.3,1x))')' Velocities [in m.s-1]: compressional Vp=',Vp,' shear Vs=',Vs,' and bulk Vphi=',Vphi
  V_Debye=(1./3.*(1./Vp**3+2./Vs**3))**(-1./3.)
  T_Debye= (V_Debye/1.d3/Bohr_Ang/1.d-13*Time_Sec)*(six*pi**2*Invar%natom_unitcell/Lattice%ucvol)**(1./3.)*Ha_K
  write(Invar%stdout,'(2(a,f9.3,1x))')' Debye velocity [in m.s-1]=',V_Debye,' and temperature [in K]=',T_Debye

! Compute the elastic anisotropy
  A_U=5.d0*GV/GR+BV/BR -6.d0
  A_B=(BV-BR)/(BV+BR)
  A_G=(GV-GR)/(GV+GR)
  write(Invar%stdout,'(a)')' '
  write(Invar%stdout,'(a)')' ========================= Elastic anisotropy ================================='
  write(Invar%stdout,'(a,f9.3)')' Elastic anisotropy index : A_U= 5*G_V/G_R + K_V/K_R - 6 =',A_U
  write(Invar%stdout,'(a,f9.3)')' Bulk anisotropy ratio : A_B= (B_V-B_R)/(B_V+B_R) =',A_B
  write(Invar%stdout,'(a,f9.3)')' Shear anisotropy ratio : A_G= (G_V-G_R)/(G_V+G_R) =',A_G


! Store the RVH value of the elastic moduli (will be useful for the Gruneisen)
  Lattice%BulkModulus_T=BH
  Lattice%Shear=GH
  Lattice%Density=rho

  ABI_FREE(cijkl)
  ABI_FREE(Cij)
  ABI_FREE(Sij)
end subroutine tdep_calc_elastic
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module m_tdep_phdos
