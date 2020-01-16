!!****m* ABINIT/m_nucprop
!! NAME
!!  m_nucprop
!!
!! FUNCTION
!!  routines used to compute properties at the nuclear sites, including
!!  electric field gradient and Fermi contact
!!
!! COPYRIGHT
!! Copyright (C) 1998-2019 ABINIT group (MT, JWZ)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_nucprop

  use defs_basis
  use m_abicore
  use m_errors

  use defs_abitypes, only : MPI_type
  use m_mpinfo,   only : ptabs_fourdp
  use m_xmpi, only : xmpi_comm_self, xmpi_sum
  use m_geometry,       only : xred2xcart
  use m_linalg_interfaces, only: dsyev
  use m_paw_an,      only : paw_an_type
  use m_pawang,      only : pawang_type
  use m_pawrad,     only : pawrad_type
  use m_pawtab,     only : pawtab_type
  use m_pawrhoij,   only : pawrhoij_type
  use m_paw_nmr,    only : make_efg_onsite,make_fc_paw
  use m_paral_atom, only : get_my_atmtab,free_my_atmtab
  use m_special_funcs,  only : abi_derfc
  use m_symtk,          only : matr3inv, matpointsym
  use m_fft,           only : fourdp

  implicit none

  private
!!***

  public :: calc_efg
  public :: calc_fc
  public :: make_efg_ion
  public :: make_efg_el
!!***

contains

!!****f* ABINIT/calc_efg
!! NAME
!! calc_efg
!!
!! FUNCTION
!! calculation and output of electric field gradient tensor at each atomic site
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  nfft=number of points on fft grid
!!  ngfft(18)=details of fft
!!  nspden=number of spin densities
!!  nsym=number of symmetries in space group
!!  ntypat=number of atom types
!!  ptcharge(ntypat)=user input charges on atoms to make simple point charge calc
!!  paw_an(my_natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  prtefg=1 to print summary output, 2 for detailed output
!!  quadmom(ntypat)=quadrupole moments in barns of different atomic nuclei
!!  rhor(nfft,nspden)=electron density on grid (strictly $\tilde{n}+\hat{n}$)
!!  rprimd(3,3)=matrix relating cartesian coordinates to crystal coordinates
!!  symrel(3,3,nsym)=symmetry operators in terms of action on primitive translations
!!  tnons(3,nsym)=nonsymmorphic translations
!!  typat(natom)=type (integer) for each atom
!!  ucvol=unit cell volume in Bohr^3
!!  usepaw=1 if we are using PAW formalism, 0 else
!!  xred(3,natom)=vectors locating each atom in the unit cell, in crystal coords
!!  zion(ntypat)=net core charge on each type of atom
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      dsyev,free_my_atmtab,get_my_atmtab,make_efg_el,make_efg_ion
!!      make_efg_onsite,wrtout
!!
!! SOURCE

  subroutine calc_efg(mpi_enreg,my_natom,natom,nfft,ngfft,nspden,nsym,ntypat,&
                      paw_an,pawang,pawrad,pawrhoij,pawtab,&
                      ptcharge,prtefg,quadmom,rhor,rprimd,symrel,tnons,typat,ucvol,usepaw,xred,zion,&
                      mpi_atmtab,comm_atom) ! optional arguments (parallelism)

    !Arguments ------------------------------------
    !scalars
    integer,intent(in) :: my_natom,natom,nfft,nspden,nsym,ntypat,prtefg,usepaw
    integer,optional,intent(in) :: comm_atom
    real(dp),intent(in) :: ucvol
    type(MPI_type),intent(in) :: mpi_enreg
    type(pawang_type),intent(in) :: pawang
    !arrays
    integer,intent(in) :: ngfft(18),symrel(3,3,nsym),typat(natom)
    integer,optional,target,intent(in) :: mpi_atmtab(:)
    real(dp),intent(in) :: ptcharge(ntypat)
    real(dp),intent(in) :: quadmom(ntypat),rhor(nfft,nspden),rprimd(3,3)
    real(dp),intent(in) :: tnons(3,nsym),zion(ntypat)
    real(dp),intent(inout) :: xred(3,natom)
    type(paw_an_type),intent(in) :: paw_an(my_natom)
    type(pawrad_type),intent(in) :: pawrad(ntypat)
    type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
    type(pawtab_type),intent(in) :: pawtab(ntypat)

    !Local variables-------------------------------
    !scalars
    integer :: INFO,LDA,LWORK,N,iatom,my_comm_atom
    logical :: my_atmtab_allocated,paral_atom
    real(dp) :: cq,eta,vxx,vyy,vzz
    character(len=500) :: message
    !arrays
    integer,pointer :: my_atmtab(:)
    real(dp) :: eigval(3),matr(3,3),work(8)
    real(dp),allocatable :: efg(:,:,:),efg_el(:,:,:),efg_ion(:,:,:),efg_paw(:,:,:)
    real(dp),allocatable :: efg_point_charge(:,:,:)

    ! ************************************************************************

    !Compatibility tests
    if (usepaw /= 1) then
       message = ' usepaw /= 1 but EFG calculation requires PAW '
       MSG_ERROR(message)
    end if

    !Set up parallelism over atoms
    paral_atom=(present(comm_atom).and.(my_natom/=natom))
    nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
    my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
    call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

    ABI_ALLOCATE(efg,(3,3,natom))
    ABI_ALLOCATE(efg_el,(3,3,natom))
    ABI_ALLOCATE(efg_ion,(3,3,natom))
    ABI_ALLOCATE(efg_paw,(3,3,natom))
    ABI_ALLOCATE(efg_point_charge,(3,3,natom))
    efg_el(:,:,:) = zero
    efg_ion(:,:,:) = zero
    efg_paw(:,:,:) = zero
    efg_point_charge(:,:,:) = zero

    call make_efg_el(efg_el,mpi_enreg,natom,nfft,ngfft,nspden,nsym,rhor,rprimd,symrel,tnons,xred)

    call make_efg_ion(efg_ion,natom,nsym,ntypat,rprimd,symrel,tnons,typat,ucvol,xred,zion)

    if (paral_atom) then
       call make_efg_onsite(efg_paw,my_natom,natom,nsym,ntypat,paw_an,pawang,pawrhoij,pawrad,pawtab,&
            &   rprimd,symrel,tnons,xred,comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
    else
       call make_efg_onsite(efg_paw,my_natom,natom,nsym,ntypat,paw_an,pawang,pawrhoij,pawrad,pawtab,&
            &   rprimd,symrel,tnons,xred)
    end if

    !calculate efg due to pure point charges, as input in variable ptcharge(ntypat)
    !note here all atoms of the same type will have the same valence; in the future this
    !could be made more flexible by having ptcharge(natom) but that will require a slightly
    !different version than the existing make_efg_ion routine
    if(prtefg > 2) then
       call make_efg_ion(efg_point_charge,natom,nsym,ntypat,rprimd,symrel,tnons,typat,ucvol,xred,ptcharge)
    end if

    efg(:,:,:) = efg_el(:,:,:) + efg_ion(:,:,:) + efg_paw(:,:,:)

    write(message,'(a,a,a)' ) ch10,' Electric Field Gradient Calculation ',ch10
    call wrtout(ab_out,message,'COLL')

    LDA=3; LWORK=8;N=3 ! these parameters are needed for the LAPACK dsyev routine
    do iatom = 1, natom
       matr(:,:) = efg(:,:,iatom)
       call dsyev('V','U',N,matr,LDA,eigval,work,LWORK,INFO) ! get eigenvalues and eigenvectors
       if (eigval(3) > abs(eigval(1)) ) then ! In NMR, the convention is that whatever component is
          !    largest in magnitude is called Vzz, next comes Vxx, then Vyy
          vzz = eigval(3)
          vxx = eigval(1)
          vyy = eigval(2)
       else
          vzz = eigval(1)
          vxx = eigval(3)
          vyy = eigval(2)
       end if
       if (abs(quadmom(typat(iatom))) > tol8 ) then ! only relevant when quadmom > 0 for a given atom
          !    cq = (eQ)*Vzz/h, where Q is the electric quadrupole moment and Vzz is the largest in magnitude
          !    principal component of the EFG tensor. Q is input in quadmom in barns, and Vzz is computed in atomic
          !    units. The factor 2349647.81 Ha^{-1}Bohr^2 fm^{-2} sec^-1 converts from atomic units to frequency (see
          !    http://www.ismar.org/ISMARpedia/index.php/Nuclear_Quadrupole_Resonance for discussion); we divide by
          !    10^6 to convert to MHz from Hz and multiply by 100 to convert from fm^2 to Barns.
          cq = vzz*quadmom(typat(iatom))*2349647.81/1.0E4
          if(abs(cq) > tol6 )then ! if Cq is non-zero, eta is meaningful, otherwise it s numerical noise
             eta = abs(vxx - vyy)/abs(vzz)
          else
             eta=zero
          end if
       else
          cq =zero
          eta =zero
       end if
       !  we always write Cq and eta, these are the NMR observables
       write(message,'(a,i3,a,i3,a,f13.6,a,f13.6)') ' Atom ',iatom,', typat ',typat(iatom),': Cq = ',cq,' MHz     eta = ',eta
       call wrtout(ab_out,message,'COLL')
       if (prtefg > 1) then ! print detailed results on component EFG's
          write(message,'(a,a,f13.6,a,a,3f13.6)')ch10,'      efg eigval : ',eigval(1),ch10,&
               &     '-         eigvec : ',matr(1,1),matr(2,1),matr(3,1)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,f13.6,a,a,3f13.6)')'      efg eigval : ',eigval(2),ch10,&
               &     '-         eigvec : ',matr(1,2),matr(2,2),matr(3,2)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,f13.6,a,a,3f13.6)')'      efg eigval : ',eigval(3),ch10,&
               &     '-         eigvec : ',matr(1,3),matr(2,3),matr(3,3)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,a,3f13.6)')ch10,'      total efg : ',efg(1,1,iatom),efg(1,2,iatom),efg(1,3,iatom)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,3f13.6)')'      total efg : ',efg(2,1,iatom),efg(2,2,iatom),efg(2,3,iatom)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,3f13.6,a)')'      total efg : ',efg(3,1,iatom),efg(3,2,iatom),efg(3,3,iatom),ch10
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,a,3f13.6)')ch10,'      efg_el : ',efg_el(1,1,iatom),efg_el(1,2,iatom),efg_el(1,3,iatom)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,3f13.6)')'      efg_el : ',efg_el(2,1,iatom),efg_el(2,2,iatom),efg_el(2,3,iatom)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,3f13.6,a)')'      efg_el : ',efg_el(3,1,iatom),efg_el(3,2,iatom),efg_el(3,3,iatom),ch10
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,3f13.6)')'      efg_ion : ',efg_ion(1,1,iatom),efg_ion(1,2,iatom),efg_ion(1,3,iatom)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,3f13.6)')'      efg_ion : ',efg_ion(2,1,iatom),efg_ion(2,2,iatom),efg_ion(2,3,iatom)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,3f13.6,a)')'      efg_ion : ',efg_ion(3,1,iatom),efg_ion(3,2,iatom),efg_ion(3,3,iatom),ch10
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,3f13.6)')'      efg_paw : ',efg_paw(1,1,iatom),efg_paw(1,2,iatom),efg_paw(1,3,iatom)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,3f13.6)')'      efg_paw : ',efg_paw(2,1,iatom),efg_paw(2,2,iatom),efg_paw(2,3,iatom)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,3f13.6,a)')'      efg_paw : ',efg_paw(3,1,iatom),efg_paw(3,2,iatom),efg_paw(3,3,iatom),ch10
          call wrtout(ab_out,message,'COLL')
       end if
       if (prtefg > 2) then ! write output of pure pointcharge calculation
          matr(:,:) = efg_point_charge(:,:,iatom)
          call dsyev('V','U',N,matr,LDA,eigval,work,LWORK,INFO) ! get eigenvalues and eigenvectors
          if (eigval(3) > abs(eigval(1)) ) then ! In NMR, the convention is that whatever component is
             !      largest in magnitude is called Vzz, next comes Vxx, then Vyy
             vzz = eigval(3)
             vxx = eigval(1)
             vyy = eigval(2)
          else
             vzz = eigval(1)
             vxx = eigval(3)
             vyy = eigval(2)
          end if
          if (abs(quadmom(typat(iatom))) > tol8 ) then ! only relevant when quadmom > 0 for a given atom
             !      cq = e2Qq/h, where Vzz = eq and quadmom = Q; the other factors convert from atomic units to MHz
             cq = vzz*quadmom(typat(iatom))*2349647.81/1.0E4
             if(abs(cq) > tol6 )then ! if Cq is non-zero, eta is meaningful, otherwise it s numerical noise
                eta = abs(vxx - vyy)/abs(vzz)
             else
                eta=zero
             end if
          else
             cq =zero
             eta =zero
          end if
          !    we always write Cq and eta, these are the NMR observables
          write(message,'(a,i3,a,i3,a,f13.6,a,f13.6)') ' Atom ',iatom,', typat ',typat(iatom),&
               &     ': Point charge Cq = ',cq,' MHz     eta = ',eta
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,a,f13.6,a,a,3f13.6)')ch10,'      point charge efg eigval : ',eigval(1),ch10,&
               &     '-         eigvec : ',matr(1,1),matr(2,1),matr(3,1)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,f13.6,a,a,3f13.6)')'      point charge efg eigval : ',eigval(2),ch10,&
               &     '-         eigvec : ',matr(1,2),matr(2,2),matr(3,2)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,f13.6,a,a,3f13.6)')'      point charge efg eigval : ',eigval(3),ch10,&
               &     '-         eigvec : ',matr(1,3),matr(2,3),matr(3,3)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,a,3f13.6)')ch10,'      point charge efg : ',efg_point_charge(1,1,iatom),&
               &     efg_point_charge(1,2,iatom),efg_point_charge(1,3,iatom)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,3f13.6)')'      point charge efg : ',efg_point_charge(2,1,iatom),&
               &     efg_point_charge(2,2,iatom),efg_point_charge(2,3,iatom)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,3f13.6,a)')'      point charge efg : ',efg_point_charge(3,1,iatom),&
               &     efg_point_charge(3,2,iatom),efg_point_charge(3,3,iatom),ch10
          call wrtout(ab_out,message,'COLL')
       end if
    end do
    write(message,'(3a)')ch10,ch10,ch10
    call wrtout(ab_out,message,'COLL')

    ABI_DEALLOCATE(efg)
    ABI_DEALLOCATE(efg_el)
    ABI_DEALLOCATE(efg_ion)
    ABI_DEALLOCATE(efg_paw)
    ABI_DEALLOCATE(efg_point_charge)

    !Destroy atom table used for parallelism
    call free_my_atmtab(my_atmtab,my_atmtab_allocated)

    !DEBUG
    !write(std_out,*)' calc_efg : exit '
    !stop
    !ENDDEBUG

  end subroutine calc_efg
!!***

!!***
!!****f* ABINIT/calc_fc
!! NAME
!! calc_fc
!!
!! FUNCTION
!! calculation and output of Fermi-contact term at each atomic site
!!
!! INPUTS
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  nspden=number of spin density components
!!  ntypat=number of atom types
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  typat(natom)=type (integer) for each atom
!!  usepaw=1 if PAW is activated
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,make_fc_paw,wrtout
!!
!! SOURCE

  subroutine calc_fc(my_natom,natom,nspden,ntypat,pawrad,pawrhoij,pawtab,typat,usepaw,&
       &                  mpi_atmtab,comm_atom) ! optional arguments (parallelism)

    !Arguments ------------------------------------
    !scalars
    integer,intent(in) :: my_natom,natom,nspden,ntypat,usepaw
    integer,optional,intent(in) :: comm_atom
    !arrays
    integer,intent(in) :: typat(natom)
    integer,optional,target,intent(in) :: mpi_atmtab(:)
    type(pawrad_type),intent(in) :: pawrad(ntypat)
    type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
    type(pawtab_type),intent(in) :: pawtab(ntypat)

    !Local variables-------------------------------
    !scalars
    integer :: iatom,my_comm_atom
    logical :: my_atmtab_allocated,paral_atom
    character(len=500) :: message
    !arrays
    integer,pointer :: my_atmtab(:)
    real(dp),allocatable :: fc(:,:)

!***********************************************************************

    !Compatibility tests
    if (usepaw /= 1) then
       message = ' usepaw /= 1 but Fermi-contact calculation requires PAW '
       MSG_ERROR(message)
    end if

    !Set up parallelism over atoms
    paral_atom=(present(comm_atom).and.(my_natom/=natom))
    nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
    my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
    call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

    !Initialization
    ABI_ALLOCATE(fc,(nspden,natom))

    !Computation
    if (paral_atom) then
       call make_fc_paw(fc,my_natom,natom,nspden,ntypat,pawrhoij,pawrad,pawtab,&
            &   comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
    else
       call make_fc_paw(fc,my_natom,natom,nspden,ntypat,pawrhoij,pawrad,pawtab)
    end if

    !Printing
    write(message,'(a,a,a)' ) ch10,' Fermi-contact Term Calculation ',ch10
    call wrtout(ab_out,message,'COLL')

    do iatom = 1, natom
       if (nspden == 2) then
          write(message,'(a,i3,a,i3,a,f12.4)') ' Atom ',iatom,', typat ',typat(iatom),': FC total = ',&
               &     fc(1,iatom)+fc(2,iatom)
          call wrtout(ab_out,message,'COLL')
          write(message,'(a,i3,a,i3,a,f12.4)') ' Atom ',iatom,', typat ',typat(iatom),': FC up - down = ',&
               &     fc(1,iatom)-fc(2,iatom)
          call wrtout(ab_out,message,'COLL')
       else
          write(message,'(a,i3,a,i3,a,f12.4)') ' Atom ',iatom,', typat ',typat(iatom),': FC = ',&
               &     fc(1,iatom)
          call wrtout(ab_out,message,'COLL')
       end if
    end do

    write(message,'(3a)')ch10,ch10,ch10
    call wrtout(ab_out,message,'COLL')

    !Memory deallocation
    ABI_DEALLOCATE(fc)

    !Destroy atom table used for parallelism
    call free_my_atmtab(my_atmtab,my_atmtab_allocated)

  end subroutine calc_fc
!!***

!!****f* ABINIT/make_efg_ion
!! NAME
!! make_efg_ion
!!
!! FUNCTION
!! compute the electric field gradient due to ionic cores
!!
!! INPUTS
!! natom, number of atoms in the unit cell
!! nsym=number of symmetries in space group
!! ntypat, the number of types of atoms in the unit cell
!! rprimd(3,3), the matrix giving the transformation from crystal to cartesian coordinates
!! symrel(3,3,nsym)=symmetry operators in terms of action on primitive translations
!! tnons(3,nsym) = nonsymmorphic translations
!! typat(natom), the type of each atom in the unit cell
!! ucvol, the volume of the unit cell in atomic units
!! xred(3,natom) the location of each atom in the cell in crystallographic coordinates
!! zion(ntypat) the net charge on each type of atom
!!
!! OUTPUT
!! efg(3,3,natom), the 3x3 efg tensors at each atomic site
!!
!! SIDE EFFECTS
!!
!! NOTES
!! This routine computes the electric field gradient, specifically the components
!! $\partial^2 V/\partial x_\alpha \partial x_\beta$ of the potential generated by the ionic cores,
!! at each atomic site in the unit cell.
!! Key references:
!! Profeta, Mauri, and Pickard, ``Accurate first principles prediction of $^{17}$O NMR parameters in
!! SiO$_2$: Assignment of the zeolite ferrierite spectrum'', J. Am. Chem. Soc. 125, 541--548 (2003) [[cite:Profeta2003]];
!! A. Honma, ``Dipolar lattice-sums with applications to the exciton bands of anthracene crystal and
!! the crystal field due to point charges'', J. Phys. Soc. Jpn. 42, 1129--1135 (1977) [[cite:Honma1977]];
!! and Kresse and Joubert, ``From ultrasoft pseudopotentials to the projector augmented wave method'',
!! Phys. Rev. B. 59, 1758--1775 (1999) [[cite:Kresse1999]]. In Kresse and Joubert's notation, the ionic cores are $n_{Zc}$;
!! these charges are given by the net core charges on the pseudoatoms. Due to otherwise slow convergence,
!! the sum over atoms is carried out by an Ewald method as detailed in the Honma reference, specifically
!! his Eq. 4.8.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine make_efg_ion(efg,natom,nsym,ntypat,rprimd,symrel,tnons,typat,ucvol,xred,zion)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: natom,nsym,ntypat
  real(dp) :: ucvol
  !arrays
  integer,intent(in) :: symrel(3,3,nsym),typat(natom)
  real(dp),intent(in) :: rprimd(3,3),tnons(3,nsym)
  real(dp),intent(in) :: zion(ntypat)
  real(dp),intent(inout) :: xred(3,natom)
  real(dp),intent(out) :: efg(3,3,natom)
  !Local variables-------------------------------
  !scalars
  integer :: iatom,ishell,ii,jatom,jj,nshell,sx,sy,sz
  real(dp) :: cph,dampfac,derfc_karg,derivs,gsq,karg
  real(dp) :: lenrho,phase,qk,rlkcut,trace,xi0
  real(dp) :: glkcut
  !arrays
  real(dp) :: cvec(3),gvec(3),gpl(3),gprimd(3,3)
  real(dp) :: rhok(3),rhored(3),rpl(3)
  real(dp),allocatable :: efg_g(:,:,:),efg_r(:,:,:)
  real(dp),allocatable :: xcart(:,:)

  ! ************************************************************************

  !DEBUG
  !write(std_out,*)' make_efg_ion : enter'
  !ENDDEBUG

  ABI_ALLOCATE(efg_g,(3,3,natom))
  ABI_ALLOCATE(efg_r,(3,3,natom))
  ABI_ALLOCATE(xcart,(3,natom))
  efg(:,:,:) = zero ! final efg tensor
  efg_g(:,:,:) = zero ! part of tensor accumulated in G space
  efg_r(:,:,:) = zero ! part of tensor accumulated in R space

  call xred2xcart(natom,rprimd,xcart,xred) ! get atomic locations in cartesian coords

  do ii = 1, 3 ! generate the lengths of the unit cell edges in atomic units
     rpl(ii) = sqrt(rprimd(1,ii)**2+rprimd(2,ii)**2+rprimd(3,ii)**2)
  end do
  xi0 = sqrt(pi/(maxval(rpl)*minval(rpl))) ! this estimate for xi0 is from Honma's paper

  call matr3inv(rprimd,gprimd) ! gprimd holds the inverse transpose of rprimd
  !remember ordering: rprimd( (x_comp,y_comp,z_comp), (edge 1, edge 2, edge 3) )
  !while gprimd( (edge 1, edge 2, edge 3),(x_comp, y_comp, z_comp) )
  do ii = 1, 3 ! generate the lengths of the reciprocal cell edges
     gpl(ii) = sqrt(gprimd(ii,1)**2+gprimd(ii,2)**2+gprimd(ii,3)**2)
  end do

  !go out enough shells such that g**2/4*xi0**2 is of order 30
  nshell = int(anint(sqrt(30.0)*xi0/(pi*minval(gpl))))
  glkcut = (0.95*nshell*two*pi*minval(gpl))**2

  do ishell = 0, nshell ! loop over shells
     do sx = -ishell, ishell
        do sy = -ishell, ishell
           do sz = -ishell, ishell
              if ( .not. (sx==0 .and. sy==0 .and. sz==0) ) then ! avoid origin
                 !          constrain to be on shell surface, not interior
                 if ( abs(sx)==ishell .or. abs(sy)==ishell .or. abs(sz)==ishell ) then
                    cvec(1)=sx;cvec(2)=sy;cvec(3)=sz
                    !            make the g vector in cartesian coords
                    gvec(:) = zero
                    do ii = 1, 3
                       do jj = 1, 3
                          gvec(ii) = gvec(ii) + gprimd(ii,jj)*cvec(jj)*two*pi
                       end do
                    end do
                    gsq = dot_product(gvec,gvec)
                    if(gsq < glkcut) then
                       dampfac = exp(-gsq/(4.0*xi0*xi0)) ! see Honma eq. 4.8
                       do iatom = 1, natom
                          do jatom = 1, natom
                             qk = zion(typat(jatom)) ! charge on neighbor atom
                             rhok = xcart(:,jatom)-xcart(:,iatom)
                             phase = dot_product(gvec,rhok)
                             cph = cos(phase)
                             do ii = 1, 3
                                do jj = 1, 3
                                   derivs = -3.0*gvec(ii)*gvec(jj)/gsq
                                   if (ii == jj) derivs = 1.0 + derivs
                                   efg_g(ii,jj,iatom) = efg_g(ii,jj,iatom) + &
                                        &                       qk*cph*derivs*dampfac
                                end do ! end loop over jj
                             end do ! end loop over ii
                          end do ! end loop over jatom
                       end do ! end loop over iatom
                    end if ! constrain to gsq < glkcut
                 end if ! end selection on shell edge
              end if ! end avoidance of origin
           end do ! end loop over sz
        end do ! end loop over sy
     end do ! end loop over sx
  end do ! end loop over ishell

  !sum in real space begins here

  !go out enough shells such that (r*xi0)**2 is of order 30
  nshell = int(anint(sqrt(30.)/(minval(rpl)*xi0)))
  rlkcut = nshell*minval(rpl)*0.95
!
  !go out enough shells so that rlkcut is of order 30 bohr
  !nshell=int(anint(30.0/minval(rpl)))
  !rlkcut = 0.95*nshell*minval(rpl)

  do ishell = 0, nshell ! total set of cells to loop over
     do sx = -ishell, ishell ! loop over all cells in each dimension
        do sy = -ishell, ishell
           do sz = -ishell, ishell
              !        constrain to shell surface, not interior
              if ( abs(sx)==ishell .or. abs(sy)==ishell .or. abs(sz)==ishell ) then
                 do jatom = 1, natom ! loop over atoms in shell cell
                    do iatom = 1, natom ! loop over atoms in central unit cell
                       if (.NOT. (jatom == iatom .AND. sx == 0 .AND. sy == 0 .AND. sz == 0)) then ! avoid self term
                          qk = zion(typat(jatom)) ! charge on each neighbor atom
                          !                ! rhored is the vector in crystal coords from neighbor to target
                          rhored(1) = xred(1,jatom) + sx - xred(1,iatom)
                          rhored(2) = xred(2,jatom) + sy - xred(2,iatom)
                          rhored(3) = xred(3,jatom) + sz - xred(3,iatom)
                          !                !  rhok is rhored in cartesian coords
                          rhok(1) = rprimd(1,1)*rhored(1)+rprimd(1,2)*rhored(2)+rprimd(1,3)*rhored(3)
                          rhok(2) = rprimd(2,1)*rhored(1)+rprimd(2,2)*rhored(2)+rprimd(2,3)*rhored(3)
                          rhok(3) = rprimd(3,1)*rhored(1)+rprimd(3,2)*rhored(2)+rprimd(3,3)*rhored(3)
                          trace = dot_product(rhok,rhok)
                          lenrho = sqrt(trace)
                          if (lenrho < rlkcut) then ! this restriction is critical as it ensures
                             !                  ! that we sum over a sphere of atoms in real space
                             !                  ! no matter what shape the unit cell has
                             karg = xi0*lenrho
                             derfc_karg = abi_derfc(karg)
                             !                  see Honma eq. 2.10 for derivation of the following damping factor
                             dampfac = (1.0+3.0/(2.0*karg*karg))*exp(-karg*karg)+3.0*sqrt(pi)*derfc_karg/(4.0*karg**3)
                             do ii = 1, 3 ! loop over tensor elements
                                do jj = 1, 3 ! loop over tensor elements
                                   derivs = -3.0*rhok(ii)*rhok(jj)/trace
                                   if(ii == jj) derivs = derivs + 1.0 ! see Honma eq 4.8 re: sign
                                   !                      accumulate real space tensor element,
                                   !                      weighted by charge of neighbor and Ewald damping factor
                                   efg_r(ii,jj,iatom) = efg_r(ii,jj,iatom) + qk*derivs*dampfac
                                end do ! end loop over jj in efg(ii,jj,iatom)
                             end do ! end loop over ii in efg(ii,jj,iatom)
                          end if ! end if statement restricting to a sphere of radius rlkcut
                       end if ! end if statement avoiding the self atom term
                    end do ! end loop over i atoms in cell
                 end do ! end loop over j atoms in cell
              end if ! end selection on outer shell of cells only
           end do ! end loop over sz cells
        end do ! end loop over sy cells
     end do ! end loop over sx cells
  end do ! end loop over shells

  !now combine the g-space and r-space parts, properly weighted (see Honma)
  do iatom = 1, natom
     do ii = 1, 3
        do jj = 1, 3
           efg(ii,jj,iatom) = four_pi*efg_g(ii,jj,iatom)/(three*ucvol)-&
                &       four*xi0**3*efg_r(ii,jj,iatom)/(three*sqrt(pi))
           !      note extra factor of two: compare Honma eq. 4.6
        end do
     end do
  end do

  ! symmetrize tensor at each atomic site using point symmetry operations
  do iatom = 1, natom
     call matpointsym(iatom,efg(:,:,iatom),natom,nsym,rprimd,symrel,tnons,xred)
  end do

  ABI_DEALLOCATE(efg_g)
  ABI_DEALLOCATE(efg_r)
  ABI_DEALLOCATE(xcart)

  !DEBUG
  !write(std_out,*)' make_efg_ion : exit '
  !stop
  !ENDDEBUG

end subroutine make_efg_ion
!!***

!!****f* ABINIT/make_efg_el
!! NAME
!! make_efg_el
!!
!! FUNCTION
!! compute the electric field gradient due to electron density
!!
!! INPUTS
!! mpi_enreg=information about MPI parallelization
!! natom, number of atoms in unit cell
!! nfft,ngfft(18), number of FFT points and details of FFT
!! nspden, number of spin components
!! nsym=number of symmetries in space group
!! rhor(nfft,nspden), valence electron density, here $\tilde{n} + \hat{n}$
!! rprimd(3,3), conversion from crystal coordinates to cartesian coordinates
!! symrel(3,3,nsym)=symmetry operators in terms of action on primitive translations
!! tnons(3,nsym) = nonsymmorphic translations
!! xred(3,natom), location of atoms in crystal coordinates.
!!
!! OUTPUT
!! efg(3,3,natom), the 3x3 efg tensor at each atomic site due to rhor
!!
!! NOTES
!! This routine computes the electric field gradient, specifically the components
!! $\partial^2 V/\partial x_\alpha \partial x_\beta$ of the potential generated by the valence
!! electrons, at each atomic site in the unit cell. Key references: Kresse and Joubert, ``From
!! ultrasoft pseudopotentials to the projector augmented wave method'', Phys. Rev. B. 59, 1758--1775 (1999) [[cite:Kresse1999]],
!! and Profeta, Mauri, and Pickard, ``Accurate first principles prediction of $^{17}$O NMR parameters in
!! SiO$_2$: Assignment of the zeolite ferrierite spectrum'', J. Am. Chem. Soc. 125, 541--548 (2003) [[cite:Profeta2003]]. This
!! routine computes the second derivatives of the potential generated by $\tilde{n}$ (see Kresse and Joubert
!! for notation, Fourier-transforming the density, doing the sum in G space, and then transforming back at
!! each atomic site. The final formula is
!! \begin{displaymath}
!! \frac{\partial^2 V}{\partial x_\alpha\partial x_\beta} = -4\pi^2\sum_G (G_\alpha G_\beta - \delta_{\alpha,\beta}G^2/3)
!! \left(\frac{\tilde{n}(G)}{\pi G^2}\right)e^{2\pi i G\cdot R}
!! \end{displaymath}
!!
!!
!! PARENTS
!!      calc_efg
!!
!! CHILDREN
!!      fourdp,matpointsym,matr3inv,ptabs_fourdp,xmpi_sum,xred2xcart
!!
!! SOURCE

subroutine make_efg_el(efg,mpi_enreg,natom,nfft,ngfft,nspden,nsym,rhor,rprimd,symrel,tnons,xred)

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: natom,nfft,nspden,nsym
  type(MPI_type),intent(in) :: mpi_enreg
  !arrays
  integer,intent(in) :: ngfft(18),symrel(3,3,nsym)
  real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3),tnons(3,nsym),xred(3,natom)
  real(dp),intent(out) :: efg(3,3,natom)

  !Local variables-------------------------------
  !scalars
  integer :: cplex,fftdir,fofg_index,iatom,i1,i2,i2_local,i23,i3,id1,id2,id3
  integer :: ierr,ig,ig2,ig3,ii,ii1,ing,jj
  integer :: me_fft,n1,n2,n3,nproc_fft,tim_fourdp
  real(dp) :: cph,derivs,phase,sph,trace
  ! type(MPI_type) :: mpi_enreg_seq
  !arrays
  integer :: id(3)
  integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
  integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
  real(dp) :: gprimd(3,3),gred(3),gvec(3),ratom(3)
  real(dp),allocatable :: fofg(:,:),fofr(:),gq(:,:),xcart(:,:)

  ! ************************************************************************

  !DEBUG
  !write(std_out,*)' make_efg_el : enter'
  !ENDDEBUG

  ABI_ALLOCATE(fofg,(2,nfft))
  ABI_ALLOCATE(fofr,(nfft))
  ABI_ALLOCATE(xcart,(3,natom))

  efg(:,:,:) = zero
  call xred2xcart(natom,rprimd,xcart,xred) ! get atomic locations in cartesian coords
  call matr3inv(rprimd,gprimd)

  tim_fourdp = 0 ! timing code, not using
  fftdir = -1 ! FT from R to G
  cplex = 1 ! fofr is real
  !here we are only interested in the total charge density including nhat, which is rhor(:,1)
  !regardless of the value of nspden. This may change in the future depending on
  !developments with noncollinear magnetization and so forth. Such a change will
  !require an additional loop over nspden.
  !Multiply by -1 to convert the electron particle density to the charge density
  fofr(:) = -rhor(:,1)

  ! Get the distrib associated with this fft_grid  See hartre.F90 for another example where
  ! this is done
  n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
  nproc_fft = mpi_enreg%nproc_fft; me_fft = mpi_enreg%me_fft
  call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

  call fourdp(cplex,fofg,fofr,fftdir,mpi_enreg,nfft,1,ngfft,tim_fourdp) ! construct charge density in G space

  ! the following loops over G vectors has been copied from hartre.F90 in order to be compatible with
  ! possible FFT parallelism

  ! In order to speed the routine, precompute the components of g
  ! Also check if the booked space was large enough...
  ABI_ALLOCATE(gq,(3,max(n1,n2,n3)))
  do ii=1,3
     id(ii)=ngfft(ii)/2+2
     do ing=1,ngfft(ii)
        ig=ing-(ing/id(ii))*ngfft(ii)-1
        gq(ii,ing)=ig
     end do
  end do
  id1=n1/2+2;id2=n2/2+2;id3=n3/2+2

  ! Triple loop on each dimension
  do i3=1,n3
     ig3=i3-(i3/id3)*n3-1
     gred(3) = gq(3,i3)

     do i2=1,n2
        ig2=i2-(i2/id2)*n2-1
        if (fftn2_distrib(i2) == me_fft) then

           gred(2) = gq(2,i2)
           i2_local = ffti2_local(i2)
           i23=n1*(i2_local-1 +(n2/nproc_fft)*(i3-1))
           ! Do the test that eliminates the Gamma point outside of the inner loop
           ii1=1
           if(i23==0 .and. ig2==0 .and. ig3==0) ii1=2

           ! Final inner loop on the first dimension (note the lower limit)
           do i1=ii1,n1
              !         gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
              gred(1) = gq(1,i1)
              gvec(1:3) = MATMUL(gprimd,gred)
              fofg_index=i1+i23
              trace = dot_product(gvec,gvec)
              do ii = 1, 3 ! sum over components of efg tensor
                 do jj = 1, 3 ! sum over components of efg tensor
                    derivs = gvec(ii)*gvec(jj) ! This term is $G_\alpha G_\beta$
                    if (ii == jj) derivs = derivs - trace/three
                    do iatom = 1, natom ! sum over atoms in unit cell
                       ratom(:) = xcart(:,iatom) ! extract location of atom iatom
                       phase = two_pi*dot_product(gvec,ratom) ! argument of $e^{2\pi i G\cdot R}$
                       cph = cos(phase)
                       sph = sin(phase)
                       efg(ii,jj,iatom) = efg(ii,jj,iatom) - &
                            &               four_pi*derivs*(fofg(1,fofg_index)*cph-fofg(2,fofg_index)*sph)/trace ! real part of efg tensor
                    end do ! end loop over atoms in cell
                 end do ! end loop over jj in V_ij
              end do ! end loop over ii in V_ij
           end do ! End loop on i1
        end if
     end do ! End loop on i2
  end do ! End loop on i3

  call xmpi_sum(efg,mpi_enreg%comm_fft,ierr)

  ! symmetrize tensor at each atomic site using point symmetry operations
  do iatom = 1, natom
     call matpointsym(iatom,efg(:,:,iatom),natom,nsym,rprimd,symrel,tnons,xred)
  end do

  ABI_DEALLOCATE(fofg)
  ABI_DEALLOCATE(fofr)
  ABI_DEALLOCATE(xcart)
  ABI_DEALLOCATE(gq)

  !DEBUG
  !write(std_out,*)' make_efg_el : exit '
  !stop
  !ENDDEBUG

end subroutine make_efg_el
!!***

end module m_nucprop
!!***
