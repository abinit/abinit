!{\src2tex{textfont=tt}}
!!****f* m_phonon_effective_potential/ddb_to_effective_potential
!!
!! NAME
!! ddb_to_effective_potential
!!
!! FUNCTION
!!  Transfert ddb into effective potential structure.
!!  Also calculate the IFC
!!
!! INPUTS
!! crystal  = number of atoms in primitive cell
!! ddb  = number of type of atoms
!! inp  = input of epigene
!! OUTPUT
!! effective_potantial = effective_potential structure to be initialized
!!
!! PARENTS
!!    epigene
!!
!! CHILDREN
!!    asria_calc,gtdyn9,ifc_init,ddb_internalstr,gtblk9,dfpt_phfrq,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ddb_to_effective_potential(crystal,ddb, effective_potential,inp)

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_dynmat
 use m_xmpi

 use m_ddb
 use m_ifc
 use m_crystal,  only : crystal_t
 use m_dynmat,  only  : asrif9,gtdyn9,cell9,canat9
 use m_epigene_dataset, only : epigene_dataset_type
 use m_effective_potential, only : effective_potential_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ddb_to_effective_potential'
 use interfaces_14_hidewrite
 use interfaces_77_ddb
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(ddb_type),intent(inout) :: ddb
 type(effective_potential_type), intent(inout) :: effective_potential
 type(crystal_t),intent(in) :: crystal
 type(epigene_dataset_type),intent(in) :: inp
!arrays

!Local variables-------------------------------
!scalar
 integer :: chneut,ia,ib,iblok,idir1,idir2,ii,ipert1,iphl1
 integer :: ipert2,irpt,ivarA,ivarB,msize,mpert,natom,nblok,rftyp,selectz
!arrays
 integer :: rfelfd(4),rfphon(4),rfstrs(4)
 real(dp):: dielt(3,3),elast_clamped(6,6),fact
 real(dp):: qphnrm(3),qphon(3,3)
 real(dp),allocatable :: blkval(:,:,:,:,:,:),d2asr(:,:,:,:,:)
 real(dp),allocatable :: d2cart(:,:,:,:,:),displ(:)
 real(dp),allocatable :: eigval(:,:),eigvec(:,:,:,:,:),instrain(:,:)
 real(dp),allocatable :: phfrq(:),zeff(:,:,:)
 character(len=500) :: message
 type(ifc_type) :: ifc

! *************************************************************************

!Initialisation of usefull values  
  natom = ddb%natom
  nblok = ddb%nblok
  mpert=natom+6
  msize=3*mpert*3*mpert;

!Tranfert the ddb into usable array (ipert and idir format like in abinit)
  ABI_ALLOCATE(blkval,(2,3,mpert,3,mpert,nblok))
  blkval = zero
  blkval = reshape(ddb%val,(/2,3,mpert,3,mpert,nblok/))

!**********************************************************************
! Transfert basics values 
!**********************************************************************
  effective_potential%natom  = crystal%natom
  effective_potential%ucvol  = crystal%ucvol
  effective_potential%rprimd = crystal%rprimd
  effective_potential%acell  = ddb%acell
  effective_potential%amu    = ddb%amu
  effective_potential%typat  = crystal%typat
  effective_potential%xcart  = crystal%xcart
  effective_potential%znucl  = crystal%znucl

!**********************************************************************
! Transfert energy from input file
!**********************************************************************
  write(message, '(2a,(80a),6a)') ch10,('=',ii=1,80),ch10,ch10,&
&     ' Extraction of the energy of the structure (unit: Hartree)',ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')
  if(inp%energy_reference==zero)then
   write(message,'(6a)')ch10,&
&    ' Warning : Energy of the reference structure is not specify in',&
&    ' the input file.',ch10,' Energy will set to zero):',ch10
    call wrtout(std_out,message,'COLL')
    call wrtout(ab_out,message,'COLL')

  else
    effective_potential%energy = inp%energy_reference
    write(message,'(a,es25.12)') ' Energy = ',&
&                  effective_potential%energy
    call wrtout(std_out,message,'COLL')
    call wrtout(ab_out,message,'COLL')
  end if

!**********************************************************************
! Dielectric Tensor and Effective Charges
!**********************************************************************
  ABI_ALLOCATE(zeff,(3,3,natom))
  rftyp   = 1 ! Blocks obtained by a non-stationary formulation.
  chneut  = 1 ! The ASR for effective charges is imposed
  selectz = 0 ! No selection of some parts of the effective charge tensor
  iblok = ddb_get_dielt_zeff(ddb,crystal,rftyp,chneut,selectz,dielt,zeff)
  if (iblok /=0) then
    effective_potential%epsilon_inf = dielt
    effective_potential%zeff = zeff
  else
    effective_potential%epsilon_inf(1,1) = one 
    effective_potential%epsilon_inf(2,2) = one 
    effective_potential%epsilon_inf(3,3) = one 
    effective_potential%zeff = zero
  end if

!**********************************************************************
! Look after the blok no. that contains the stress tensor
!**********************************************************************
  write(message, '(a,a,(80a),a,a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&   ' Extraction of the stress tensor (unit: GPa) and forces (unit: Ha/bohr)'
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  effective_potential%internal_stress = zero
  effective_potential%forces = zero

  qphon(:,1)=zero
  qphnrm(1)=zero
  rfphon(1:2)=0
  rfelfd(1:2)=0
  rfstrs(1:2)=0
  rftyp=4

  call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

  if (iblok /=0) then
!  firts give the corect stress values store in hartree
!  diagonal parts
   effective_potential%internal_stress(1)=blkval(1,1,natom+3,1,1,iblok)
   effective_potential%internal_stress(2)=blkval(1,2,natom+3,1,1,iblok)
   effective_potential%internal_stress(3)=blkval(1,3,natom+3,1,1,iblok)
!  the shear parts
   effective_potential%internal_stress(4)=blkval(1,1,natom+4,1,1,iblok)
   effective_potential%internal_stress(5)=blkval(1,2,natom+4,1,1,iblok)
   effective_potential%internal_stress(6)=blkval(1,3,natom+4,1,1,iblok)

!  Get forces
   effective_potential%forces(:,1:natom) = blkval(1,:,1:natom,1,1,iblok)

   write(message, '(3a)' )ch10,&
&   ' Cartesian components of forces (hartree/bohr)',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   do ii = 1, natom
   write(message, '(I4,a,3(e16.8))' ) &
&   ii,'   ',effective_potential%forces(:,ii)

   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
     
   end do

   write(message, '(a,a)' )ch10,&
&   ' Cartesian components of stress tensor (hartree/bohr^3)'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(1 1)=',effective_potential%internal_stress(1),&
&   '  sigma(3 2)=',effective_potential%internal_stress(4)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(2 2)=',effective_potential%internal_stress(2),&
&   '  sigma(3 1)=',effective_potential%internal_stress(5)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(3 3)=',effective_potential%internal_stress(3),&
&   '  sigma(2 1)=',effective_potential%internal_stress(6)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a)' ) ' '
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

 else
   
    write(message,'(2a)')ch10,&
&    ' Warning : Stress Tensor and forces are set to zero (not available in the DDB)'
    call wrtout(std_out,message,'COLL')
    call wrtout(ab_out,message,'COLL')
   
  end if

!**********************************************************************
! Elastic tensors at Gamma Point
!**********************************************************************
  write(message, '(a,a,(80a),a,a,a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&   ' Extraction of the clamped elastic tensor (unit:10^2GPa)',ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

! look after the blok no.iblok that contains the elastic tensor
  qphon(:,1)=zero
  qphnrm(1)=zero
  rfphon(1:2)=0
  rfelfd(1:2)=0
  rfstrs(1:2)=3 ! Need uniaxial  both stresses and  shear stresses
  rftyp=1 ! Blocks obtained by a non-stationary formulation.
! for both diagonal and shear parts
  call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)

  if (iblok /=0) then
!   extraction of the elastic constants from the blkvals (GPa)
    do ivarA=1,6
      do ivarB=1,6
!       because the elastic constant is 6*6,
!       so we should judge if the idir is larger than 3
!       or not
        if(ivarA>3) then
          idir1=ivarA-3
          ipert1=natom+4  !for the shear modulus
        else if(ivarA<=3) then
          idir1=ivarA
          ipert1=natom+3  !for the diagonal part
        end if
        if(ivarB>3) then
          idir2=ivarB-3
          ipert2=natom+4  !for the shear modulus
        else if(ivarB<=3) then
          idir2=ivarB
          ipert2=natom+3  !for the diagonal part
        end if
        elast_clamped(ivarA,ivarB) = blkval(1,idir1,ipert1,idir2,ipert2,iblok)/crystal%ucvol
      end do
    end do
    fact=HaBohr3_GPa
    do ivarA=1,6
      write(message,'(6f12.7)')elast_clamped(ivarA,1)*fact/100.00_dp,&
&                              elast_clamped(ivarA,2)*fact/100.00_dp,&
&                              elast_clamped(ivarA,3)*fact/100.00_dp,&
&                              elast_clamped(ivarA,4)*fact/100.00_dp,&
&                              elast_clamped(ivarA,5)*fact/100.00_dp,&
&                              elast_clamped(ivarA,6)*fact/100.00_dp
    call wrtout(std_out,message,'COLL')
    call wrtout(ab_out,message,'COLL')
    end do
    
!   Set the clamped tensor into the effective potentiel
    effective_potential%elastic_constants = elast_clamped

  else
    
    write(message,'(3a)')ch10,&
&    ' Warning : Elastic Tensor is set to zero (not available in the DDB)'
    call wrtout(std_out,message,'COLL')
    call wrtout(ab_out,message,'COLL')

!   Set the clamped tensor to zero into the effective potentiel (not available in the DDB)
    effective_potential%elastic_constants = zero
  end if

!**********************************************************************
!   Acoustic Sum Rule
!***************************************************************************
! ASR-correction (d2asr) has to be determined here from the Dynamical matrix at Gamma.
  ABI_CALLOC(d2asr,(2,3,natom,3,natom))

  write(message, '(a,a,(80a),a,a,a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of acoustic sum rule',ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

! Find the Gamma block in the DDB (no need for E-field entries)
  qphon(:,1)=zero
  qphnrm(1)=zero
  rfphon(1:2)=1
  rfelfd(:)=0
  rfstrs(:)=0
  rftyp=inp%rfmeth

  call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)
  
  d2asr = zero
  if (iblok /=0) then
    call asria_calc(inp%asr,d2asr,ddb%val(:,:,iblok),ddb%mpert,ddb%natom)
  end if


!**********************************************************************
! Interatomic Forces Calculation
!**********************************************************************
! ifc to be calculated for interpolation
  write(message, '(a,a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of the interatomic forces ',ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  call ifc_init(ifc,crystal,ddb,inp%brav,inp%asr,inp%symdynmat,inp%dipdip,inp%rfmeth,&
&   inp%ngqpt(1:3),inp%nqshft,inp%q1shft,dielt,effective_potential%zeff,inp%nsphere,inp%rifcsph,&
&   inp%prtsrlr,inp%enunit,prtfreq=.True.)

!Reorder cell from canonical coordinates to reduced coordinates (for epigene)
 call cell9(ifc%atmfrc,inp%brav,ifc%cell,ddb%gprim,effective_potential%natom,ifc%nrpt,&
&           ifc%rcan,ifc%rpt,ifc%wghatm,crystal%xred)

! Apply weight on each R point
  do irpt=1,ifc%nrpt
    do ia=1,effective_potential%natom 
      do ib=1,effective_potential%natom 
        ifc%atmfrc(:,:,ia,:,ib,irpt) = ifc%atmfrc(:,:,ia,:,ib,irpt)*ifc%wghatm(ia,ib,irpt) 
      end do
    end do
  end do

  effective_potential%ifcs%atmfrc = ifc%atmfrc    !depend of dipdip (see help)

  if (inp%dipdip == 1) then
    effective_potential%ifcs%short_atmfrc = ifc%atmfrc ! if dipdip==1 atmfrc = short range
  else
    effective_potential%ifcs%short_atmfrc = zero
  end if
  effective_potential%ifcs%ewald_atmfrc = zero       ! fill with 0
  effective_potential%ifcs%cell   = ifc%cell         ! indexation of cells
  effective_potential%ifcs%nrpt   = ifc%nrpt         ! new number of cell
  effective_potential%ifcs%wghatm = ifc%wghatm       ! weight of r point

!**********************************************************************
! Dynamical matrix calculation for each qpoint
!**********************************************************************

 !  ABI_ALLOCATE(d2cart,(2,3,mpert,3,mpert))
 !  ABI_ALLOCATE(displ,(2*3*natom*3*natom))
 !  ABI_ALLOCATE(eigval,(3,natom))
 !  ABI_ALLOCATE(eigvec,(2,3,natom,3,natom))
 !  ABI_ALLOCATE(phfrq,(3*natom))

 !  write(message,'(a,(80a),3a)')ch10,('=',ii=1,80),ch10,ch10,&
 ! &   ' Calculation of dynamical matrix '
 !  call wrtout(ab_out,message,'COLL')
 !  call wrtout(std_out,message,'COLL')

! Transfer value in effective_potential structure
!  effective_potential%qph1l = inp%qph1l

!  do iphl1=1,inp%nph1l

!    ! Initialisation of the phonon wavevector
!     qphon(:,1)=inp%qph1l(:,iphl1)
!     if (inp%nph1l /= 0) qphnrm(1) = inp%qnrml1(iphl1)

!     ! Get d2cart using the interatomic forces and the
!     ! long-range coulomb interaction through Ewald summation
!     call gtdyn9(ddb%acell,ifc%atmfrc,ifc%dielt,ifc%dipdip,ifc%dyewq0,d2cart,crystal%gmet,&
! &     ddb%gprim,mpert,natom,ifc%nrpt,qphnrm(1),qphon(:,1),crystal%rmet,ddb%rprim,ifc%rpt,&
! &     ifc%trans,crystal%ucvol,ifc%wghatm,crystal%xred,zeff)

!     ! Calculation of the eigenvectors and eigenvalues of the dynamical matrix
!     call dfpt_phfrq(ddb%amu,displ,d2cart,eigval,eigvec,crystal%indsym,&
! &     mpert,crystal%nsym,natom,crystal%nsym,crystal%ntypat,phfrq,qphnrm(1),qphon,&
! &     crystal%rprimd,inp%symdynmat,crystal%symrel,crystal%symafm,crystal%typat,crystal%ucvol)

!     ! Write the phonon frequencies
!     call dfpt_prtph(displ,inp%eivec,inp%enunit,ab_out,natom,phfrq,qphnrm(1),qphon)

!     effective_potential%dynmat(:,:,:,:,:,iphl1) = d2cart(:,:,:natom,:,:natom)
!     effective_potential%phfrq(:,iphl1) = phfrq(:) * Ha_cmm1
    
!   end do

!   ABI_DEALLOCATE(d2cart)
!   ABI_DEALLOCATE(displ)
!   ABI_DEALLOCATE(eigval)
!   ABI_DEALLOCATE(eigvec)
!   ABI_DEALLOCATE(phfrq)

!**********************************************************************
! Internal strain tensors at Gamma point
!**********************************************************************
  write(message, '(a,a,(80a),a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&   ' Calculation of the internal-strain  tensor'
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')
  ABI_ALLOCATE(instrain,(3*natom,6))
! looking after the no. of blok that contains the internal strain tensor
  qphon(:,1)=zero
  qphnrm(1)=zero
  rfphon(1:2)=0
  rfelfd(1:2)=0
  rfstrs(1:2)=3
  rftyp=1
  call gtblk9(ddb,iblok,qphon,qphnrm,rfphon,rfelfd,rfstrs,rftyp)
  if (iblok /=0) then

!    then print the internal stain tensor
    call ddb_internalstr(inp%asr,ddb%val,d2asr,iblok,instrain,ab_out,mpert,natom,nblok)

    do ipert1=1,6
      do ipert2=1,natom
        do idir2=1,3
          ii=3*(ipert2-1)+idir2
          effective_potential%internal_strain(ipert1,ipert2,idir2) = instrain(ii,ipert1)
        end do
      end do
    end do
  else
    write(message,'(3a)')ch10,&
&    ' Warning : Internal strain is set to zero (not available in the DDB)'
    call wrtout(std_out,message,'COLL')
    call wrtout(ab_out,message,'COLL')
  end if
!   DEALLOCATION OF ARRAYS
  ABI_DEALLOCATE(blkval)
  ABI_DEALLOCATE(zeff)
  ABI_DEALLOCATE(instrain)

!-------------------------------------------------------------------------------------
  write(message,'(a)')ch10
  call wrtout(std_out,message,'COLL')
  call wrtout(ab_out,message,'COLL')

end subroutine ddb_to_effective_potential
!!***
