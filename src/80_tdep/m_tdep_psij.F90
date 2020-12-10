
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_psij

  use defs_basis
  use m_errors
  use m_abicore
  use m_numeric_tools
  use m_xmpi
  use m_symkpt

  use m_crystal,          only : crystal_t
  use m_ddb,              only : ddb_type
  use m_ifc,              only : ifc_type
  use m_tdep_qpt,         only : Qpoints_type
  use m_tdep_readwrite,   only : Input_Variables_type
  use m_tdep_latt,        only : Lattice_Variables_type
  use m_tdep_shell,       only : Shell_Variables_type
  use m_tdep_sym,         only : Symetries_Variables_type
  use m_tdep_phij,        only : Eigen_Variables_type,Eigen_Variables_type,tdep_init_eigen2nd,tdep_destroy_eigen2nd
  use m_tdep_utils,       only : Coeff_Moore_type

  implicit none

  public :: tdep_calc_psijfcoeff
  public :: tdep_calc_psijtot
  public :: tdep_build_psij333
  public :: tdep_calc_gruneisen
  public :: tdep_calc_alpha_gamma

contains

!====================================================================================================
subroutine tdep_calc_psijfcoeff(CoeffMoore,InVar,proj,Shell3at,Sym,ucart)

  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell3at
  type(Coeff_Moore_type), intent(inout) :: CoeffMoore
  double precision, intent(in) :: ucart(3,InVar%natom,InVar%nstep)
  double precision, intent(in) :: proj(27,27,Shell3at%nshell)

  integer :: ishell,ncoeff,ncoeff_prev,istep,iatom,jatom,katom
  integer :: icoeff,isym,iatshell,trans
  integer :: mu,nu,xi,alpha,beta,gama
  double precision :: terme,temp
  double precision :: udiff_ki(3),udiff_ji(3),SSSu(3,27)
  double precision, allocatable :: SSS_ref(:,:,:,:,:,:)

! For each couple of atoms, transform the Psij (3x3x3) ifc matrix using the symetry operation (S)
! Note: iatom=1 is excluded in order to take into account the atomic sum rule (see below)
  ABI_MALLOC(SSS_ref,(3,27,3,3,Sym%nsym,6)); SSS_ref(:,:,:,:,:,:)=zero
  do isym=1,Sym%nsym
    do mu=1,3
      do alpha=1,3
        do nu=1,3
          do beta=1,3
            do xi=1,3
              do gama=1,3
                temp=Sym%S_ref(mu,alpha,isym,1)*Sym%S_ref(nu,beta,isym,1)*Sym%S_ref(xi,gama,isym,1)
                SSS_ref(mu,gama+(beta-1)*3+(alpha-1)*9,nu,xi,isym,1)=temp !\Psi_efg
                SSS_ref(mu,beta+(gama-1)*3+(alpha-1)*9,nu,xi,isym,2)=temp !\Psi_egf
                SSS_ref(mu,gama+(alpha-1)*3+(beta-1)*9,nu,xi,isym,3)=temp !\Psi_feg
                SSS_ref(mu,alpha+(gama-1)*3+(beta-1)*9,nu,xi,isym,4)=temp !\Psi_fge
                SSS_ref(mu,beta+(alpha-1)*3+(gama-1)*9,nu,xi,isym,5)=temp !\Psi_gef
                SSS_ref(mu,alpha+(beta-1)*3+(gama-1)*9,nu,xi,isym,6)=temp !\Psi_gfe
              end do
            end do
          end do
        end do
      end do
    end do
  end do

  write(InVar%stdout,*) ' Compute the coefficients (at the 3rd order) used in the Moore-Penrose...'
  do ishell=1,Shell3at%nshell
    do iatom=1,InVar%natom
      if (Shell3at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
      do iatshell=1,Shell3at%neighbours(iatom,ishell)%n_interactions
        jatom=Shell3at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
        katom=Shell3at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
        if (iatom==jatom.or.iatom==katom) cycle
!FB        if (iatom==katom) cycle
        isym =Shell3at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
        trans=Shell3at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
        ncoeff     =Shell3at%ncoeff(ishell)
        ncoeff_prev=Shell3at%ncoeff_prev(ishell)+CoeffMoore%ncoeff2nd+CoeffMoore%ncoeff1st

        do istep=1,InVar%nstep

!         In order to impose the acoustic sum rule we use :
!         (u_k^\gamma-u_i^\gamma).u_j^\beta = udiff^\gamma.u_j^\beta
          udiff_ki(1)=ucart(1,katom,istep)-ucart(1,iatom,istep)
          udiff_ki(2)=ucart(2,katom,istep)-ucart(2,iatom,istep)
          udiff_ki(3)=ucart(3,katom,istep)-ucart(3,iatom,istep)
          udiff_ji(1)=ucart(1,jatom,istep)-ucart(1,iatom,istep)
          udiff_ji(2)=ucart(2,jatom,istep)-ucart(2,iatom,istep)
          udiff_ji(3)=ucart(3,jatom,istep)-ucart(3,iatom,istep)

!         F_i^{\mu}(t)=\sum_{\alpha\beta,jk,\nu\xi} S^{\mu\alpha}.S^{\nu\beta}.S^{\xi\gamma}.\Psi_{ijk}^{\alpha\beta\gamma}.udiff^\xi(t).u_j^\nu
          SSSu(:,:)=zero
          do nu=1,3
            do xi=1,3
!FB              SSSu(:,:)=SSSu(:,:)+SSS_ref(:,:,nu,xi,isym,trans)*udiff(xi)*ucart(nu,jatom,istep)
              SSSu(:,:)=SSSu(:,:)+SSS_ref(:,:,nu,xi,isym,trans)*udiff_ki(xi)*udiff_ji(nu)
!FB              SSSu(:,:)=SSSu(:,:)+SSS_ref(:,:,nu,xi,isym,trans)*(ucart(nu,jatom,istep)*ucart(xi,katom,istep)+ucart(nu,iatom,istep)*ucart(xi,iatom,istep))
!FB              SSSu(:,:)=SSSu(:,:)+SSS_ref(:,:,nu,xi,isym,trans)*udiff_ki(xi)*ucart(nu,jatom,istep)
!FB              SSSu(:,:)=SSSu(:,:)+SSS_ref(:,:,nu,xi,isym,trans)*ucart(nu,jatom,istep)*ucart(xi,katom,istep)
            end do
          end do
          do mu=1,3
            do icoeff=1,ncoeff
              terme=sum(SSSu(mu,:)*proj(:,icoeff,ishell))
              CoeffMoore%fcoeff(mu+3*(iatom-1)+3*InVar%natom*(istep-1),icoeff+ncoeff_prev)= &
&             CoeffMoore%fcoeff(mu+3*(iatom-1)+3*InVar%natom*(istep-1),icoeff+ncoeff_prev)+terme
            end do
          end do

        end do !istep
      end do !iatshell
    end do !iatom
  end do !ishell
  write(InVar%stdout,*) ' ------- achieved'
  ABI_FREE(SSS_ref)

end subroutine tdep_calc_psijfcoeff

!=====================================================================================================
subroutine tdep_calc_psijtot(distance,InVar,ntotcoeff,proj,Psij_coeff,Psij_ref,Shell3at,Sym)

  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell3at
  integer,intent(in) :: ntotcoeff
  double precision, intent(in) :: distance(InVar%natom,InVar%natom,4),proj(27,27,Shell3at%nshell)
  double precision, intent(in) :: Psij_coeff(ntotcoeff,1)
  double precision, intent(inout) :: Psij_ref(3,3,3,Shell3at%nshell)

  integer :: ishell,jshell,isym,jatom,katom,eatom,fatom,gatom,ncoeff,ncoeff_prev
  integer :: iatref,jatref,katref,iatshell,nshell,trans
  integer :: ii,jj,kk,kappa !,alpha,beta,gama
  double precision, allocatable :: Psij_333(:,:,:)

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '#### For each shell, list of coefficients (IFC), number of neighbours... ####'
  write(InVar%stdout,*) '#############################################################################'

  nshell=Shell3at%nshell
!==========================================================================================
!======== 1/ Build the Psi_ijk ============================================================
!==========================================================================================
  ABI_MALLOC(Psij_333,(3,3,3)) ; Psij_333(:,:,:)=0.d0
  do ishell=1,nshell
!   Build the 3x3x3 IFC per shell
    ncoeff     =Shell3at%ncoeff(ishell)
    ncoeff_prev=Shell3at%ncoeff_prev(ishell)
    kappa=0
    do ii=1,3
      do jj=1,3
        do kk=1,3
          kappa=kappa+1
          Psij_ref(ii,jj,kk,ishell)=sum(proj(kappa,1:ncoeff,ishell)*Psij_coeff(ncoeff_prev+1:ncoeff_prev+ncoeff,1))
        end do
      end do
    end do
  end do
  do eatom=1,Invar%natom
    do ishell=1,nshell
!     Build the 3x3x3 IFC of an atom in this shell
      if (Shell3at%neighbours(eatom,ishell)%n_interactions.eq.0) cycle
      do iatshell=1,Shell3at%neighbours(eatom,ishell)%n_interactions
        fatom=Shell3at%neighbours(eatom,ishell)%atomj_in_shell(iatshell)
        gatom=Shell3at%neighbours(eatom,ishell)%atomk_in_shell(iatshell)
        isym =Shell3at%neighbours(eatom,ishell)%sym_in_shell(iatshell)
        trans=Shell3at%neighbours(eatom,ishell)%transpose_in_shell(iatshell)
!       Skip the onsite terms
        if ((eatom.eq.fatom).and.(eatom.eq.gatom)) cycle
        call tdep_build_psij333(isym,InVar,Psij_ref(:,:,:,ishell),Psij_333,Sym,trans)
!       Compute the onsite terms : Psi_iii=sum_{k.ne.i j.ne.i} Psi_ijk
        if (fatom.ne.eatom.and.gatom.ne.eatom) then
          do jshell=1,Shell3at%nshell
            iatref=Shell3at%iatref(jshell)
            jatref=Shell3at%jatref(jshell)
            katref=Shell3at%katref(jshell)
            if ((iatref.eq.jatref).and.(jatref.eq.katref).and.(katref.eq.eatom)) then
              Psij_ref(:,:,:,jshell)=Psij_ref(:,:,:,jshell)+Psij_333(:,:,:)
            end if
          end do !jshell
        end if
      end do !iatshell
    end do !ishell
  end do !eatom

!==========================================================================================
!======== 2/ Write the Psi_ijk in output ==================================================
!==========================================================================================
! Remove the rounding errors before writing (for non regression testing purposes)
  do ishell=1,Shell3at%nshell
    do ii=1,3
      do jj=1,3
        do kk=1,3
          if (abs(Psij_ref(ii,jj,kk,ishell)).lt.tol8) Psij_ref(ii,jj,kk,ishell)=zero
        end do
      end do
    end do
  end do

! Write the IFCs in the data.out file (with others specifications:
! number of atoms in a shell, Trace...)
  do ishell=1,Shell3at%nshell
    iatref=Shell3at%iatref(ishell)
    if (Shell3at%neighbours(iatref,ishell)%n_interactions.ne.0) then
      jatref=Shell3at%jatref(ishell)
      katref=Shell3at%katref(ishell)
      write(InVar%stdout,'(a,i4,a,i4,a)') ' ======== NEW SHELL (ishell=',ishell,&
&           '): There are',Shell3at%neighbours(iatref,ishell)%n_interactions,' atoms on this shell'
      do iatshell=1,Shell3at%neighbours(iatref,ishell)%n_interactions
        jatom=Shell3at%neighbours(iatref,ishell)%atomj_in_shell(iatshell)
        katom=Shell3at%neighbours(iatref,ishell)%atomk_in_shell(iatshell)
        isym =Shell3at%neighbours(iatref,ishell)%sym_in_shell(iatshell)
        trans=Shell3at%neighbours(iatref,ishell)%transpose_in_shell(iatshell)
!FB        if ((iatref.ne.jatom).or.(jatom.ne.katom)) then
          call tdep_build_psij333(isym,InVar,Psij_ref(:,:,:,ishell),Psij_333,Sym,trans)
!FB	else
!FB	  Psij_333(:,:,:)=Psij_ref(:,:,:,ishell)
!FB	end if
        write(InVar%stdout,'(a,i4,a,i4)') '  For iatcell=',iatref,' ,with type=',mod(iatref-1,InVar%natom_unitcell)+1
        write(InVar%stdout,'(a,i4,a,i4)') '  For jatom  =',jatom ,' ,with type=',mod(jatom -1,InVar%natom_unitcell)+1
        write(InVar%stdout,'(a,i4,a,i4)') '  For katom  =',katom ,' ,with type=',mod(katom -1,InVar%natom_unitcell)+1
        write(InVar%stdout,'(a)') '  \Psi^{\alpha\beta x}='
        do ii=1,3
          write(InVar%stdout,'(2x,3(f9.6,1x))') (Psij_333(ii,jj,1),jj=1,3)
        end do
        write(InVar%stdout,'(a)') '  \Psi^{\alpha\beta y}='
        do ii=1,3
          write(InVar%stdout,'(2x,3(f9.6,1x))') (Psij_333(ii,jj,2),jj=1,3)
        end do
        write(InVar%stdout,'(a)') '  \Psi^{\alpha\beta z}='
        do ii=1,3
          write(InVar%stdout,'(2x,3(f9.6,1x))') (Psij_333(ii,jj,3),jj=1,3)
        end do
        write(InVar%stdout,'(a,3(f9.6,1x))') '  (i,j) vector components:', (distance(iatref,jatom,jj+1),jj=1,3)
        write(InVar%stdout,'(a,3(f9.6,1x))') '  (j,k) vector components:', (distance(jatom ,katom,jj+1),jj=1,3)
        write(InVar%stdout,'(a,3(f9.6,1x))') '  (k,i) vector components:', (distance(katom,iatref,jj+1),jj=1,3)
        write(InVar%stdout,*) ' '
      end do !iatshell
    end if !n_interactions
  end do !ishell
  ABI_FREE(Psij_333)

end subroutine tdep_calc_psijtot

!=====================================================================================================
subroutine tdep_calc_gruneisen(distance,Eigen2nd,Gruneisen,iqpt,InVar,Lattice,Psij_ref,qpt_cart,Rlatt_cart,Shell3at,Sym)

  type(Symetries_Variables_type),intent(in) :: Sym
  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(in) :: Shell3at
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Eigen_Variables_type),intent(in) :: Eigen2nd
  integer,intent(in) :: iqpt
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision,intent(in) :: Psij_ref(3,3,3,Shell3at%nshell)
  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)
  double precision,intent(in) :: qpt_cart(3)
  double complex  ,intent(out) :: Gruneisen(3*InVar%natom_unitcell)

  integer :: ii,jj,kk,iatcell,jatom,katom,itypat,jtypat,natom_unitcell
  integer :: imode,nmode,ncomp,jat_mod,jatcell,iatshell,ishell,isym,trans
  double precision :: phase
  double precision, allocatable :: Psij_333(:,:,:)
  double complex, allocatable :: eigen_prod(:,:,:,:,:),Grun_shell(:,:)

  ABI_UNUSED(Lattice%brav)

! Define quantities
  natom_unitcell=InVar%natom_unitcell
  nmode=3*InVar%natom_unitcell
  ncomp=3

! Calculation of the eigenvectors product
  ABI_MALLOC(eigen_prod,(ncomp,ncomp,natom_unitcell,natom_unitcell,nmode)) ; eigen_prod(:,:,:,:,:)=zero
  do iatcell=1,natom_unitcell
    do jatcell=1,natom_unitcell
      do imode=1,nmode
        do ii=1,3
          do jj=1,3
            eigen_prod(ii,jj,iatcell,jatcell,imode)=conjg(Eigen2nd%eigenvec((iatcell-1)*3+ii,imode,iqpt))*&
&                                                         Eigen2nd%eigenvec((jatcell-1)*3+jj,imode,iqpt)
          end do !ii
        end do !jj
      end do !imode
    end do !jatcell
  end do !iatcell

! Calculation of the Gruneisen
  ABI_MALLOC(Grun_shell,(3*InVar%natom_unitcell,Shell3at%nshell)); Grun_shell(:,:)=czero
  ABI_MALLOC(Psij_333,(3,3,3)) ; Psij_333(:,:,:)=0.d0
  Gruneisen(:)=zero
  do ishell=1,Shell3at%nshell
    do iatcell=1,InVar%natom_unitcell
      itypat=InVar%typat_unitcell(iatcell)
      do iatshell=1,Shell3at%neighbours(iatcell,ishell)%n_interactions
        jatom=Shell3at%neighbours(iatcell,ishell)%atomj_in_shell(iatshell)
        katom=Shell3at%neighbours(iatcell,ishell)%atomk_in_shell(iatshell)
        isym =Shell3at%neighbours(iatcell,ishell)%sym_in_shell(iatshell)
        trans=Shell3at%neighbours(iatcell,ishell)%transpose_in_shell(iatshell)
        jat_mod=mod(jatom+InVar%natom_unitcell-1,InVar%natom_unitcell)+1
        jtypat=InVar%typat_unitcell(jat_mod)
        phase=0.d0
        do jj=1,3
          phase=phase+2*pi*Rlatt_cart(jj,iatcell,jatom)*qpt_cart(jj)
        end do
        call tdep_build_psij333(isym,InVar,Psij_ref(:,:,:,ishell),Psij_333,Sym,trans)
        do ii=1,3
          do jj=1,3
            do kk=1,3
              do imode=1,nmode
                Grun_shell(imode,ishell)=Grun_shell(imode,ishell)-dcmplx(Psij_333(ii,jj,kk),0.d0)*&
&                                eigen_prod(ii,jj,iatcell,jat_mod,imode)*&
&                                exp(dcmplx(0.d0,phase))*dcmplx(distance(iatcell,katom,kk+1),0.d0)/&
&                                dcmplx(dsqrt(InVar%amu(itypat)*InVar%amu(jtypat))*amu_emass,0.d0)/&
&                                dcmplx(6*Eigen2nd%eigenval(imode,iqpt)**2,0.d0)
              end do !imode
            end do !kk
          end do !jj
        end do !ii
      end do !iatshell
    end do !iatcell
  end do !ishell
  ABI_FREE(Psij_333)
  ABI_FREE(eigen_prod)
  do ishell=1,Shell3at%nshell
    Gruneisen(:)=Gruneisen(:)+Grun_shell(:,ishell)
  end do
  ABI_FREE(Grun_shell)

end subroutine tdep_calc_gruneisen

!=====================================================================================================
subroutine tdep_build_psij333(isym,InVar,Psij_ref,Psij_333,Sym,trans)

  type(Symetries_Variables_type),intent(in) :: Sym
  type(Input_Variables_type),intent(in) :: InVar
  double precision, intent(in) :: Psij_ref(3,3,3)
  double precision, intent(out) :: Psij_333(3,3,3)
  integer,intent(in) :: isym,trans

  integer :: alpha,beta,gama
  integer :: ii,jj,kk,ee,ff,gg,mu,nu,xi
  double precision :: Psij_tmp(3,3,3)

  ABI_UNUSED(invar%natom)

! Transform in the new basis wrt S_ref
  Psij_333(:,:,:)=zero
  do mu=1,3
    do alpha=1,3
      do nu=1,3
        do beta=1,3
          do xi=1,3
            do gama=1,3
              Psij_333(mu,nu,xi)=Psij_333(mu,nu,xi)+&
&             Sym%S_ref(mu,alpha,isym,1)*Sym%S_ref(nu,beta,isym,1)*Sym%S_ref(xi,gama,isym,1)*Psij_ref(alpha,beta,gama)
            end do
          end do
        end do
      end do
    end do
  end do

! Take into account the 6 allowed permutations
  Psij_tmp(:,:,:)=Psij_333(:,:,:)
  if ((trans.lt.1).or.(trans.gt.6)) then
    ABI_BUG('This value of the symmetry index is not permitted')
  end if
  do ii=1,3
    do jj=1,3
      do kk=1,3
        if (trans==1) then ; ee=ii ; ff=jj ; gg=kk ; endif !\Psi_efg
        if (trans==2) then ; ee=ii ; ff=kk ; gg=jj ; endif !\Psi_egf
        if (trans==3) then ; ee=jj ; ff=ii ; gg=kk ; endif !\Psi_feg
        if (trans==4) then ; ee=jj ; ff=kk ; gg=ii ; endif !\Psi_fge
        if (trans==5) then ; ee=kk ; ff=ii ; gg=jj ; endif !\Psi_gef
        if (trans==6) then ; ee=kk ; ff=jj ; gg=ii ; endif !\Psi_gfe
        Psij_333(ee,ff,gg)=Psij_tmp(ii,jj,kk)
      end do
    end do
  end do

end subroutine tdep_build_psij333

!=====================================================================================================
subroutine tdep_calc_alpha_gamma(Crystal,distance,DDB,Ifc,InVar,Lattice,Psij_ref,Rlatt_cart,Shell3at,Sym)

  type(crystal_t),intent(in) :: Crystal
  type(ddb_type),intent(in) :: DDB
  type(ifc_type),intent(in) :: Ifc
  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Shell_Variables_type),intent(in) :: Shell3at
  type(Symetries_Variables_type),intent(in) :: Sym
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision,intent(in) :: Psij_ref(3,3,3,Shell3at%nshell)
  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)

  integer :: iq_ibz,nqbz,nqibz,ii,iatcell,jatcell,jj,nmode,ntemp,itemp
  integer, allocatable :: ibz2bz(:), bz2ibz_smap(:,:)
  double precision :: k_B,wovert,xx,C_v,Gama,alpha_v,P_th,E_th,P_th1,P_th2
  double precision, allocatable :: qbz(:,:),wtq(:),wtq_folded(:),wtqibz(:),qibz(:,:),qibz_cart(:,:)
  double precision, allocatable :: omega(:,:),displ(:,:),out_eigvec(:,:,:,:,:,:),heatcapa(:),grun_thermo(:),u_vib(:)
  double precision, allocatable :: heatcapa_HA(:,:),grun_thermo_HA(:,:),p_thermo1(:),p_thermo_HA(:,:),u_vib_HA(:,:)
  double precision, allocatable :: p_thermo2(:),tmp(:)
  double complex :: Gruneisen(3*InVar%natom_unitcell)
  type(Eigen_Variables_type) :: Eigen2nd_tmp

! Initialize useful quantities stored in the DDB file
  nqbz=DDB%nblok
  ABI_MALLOC(qbz,(3,nqbz)); qbz(:,:)=zero
  qbz(:,:)=DDB%qpt(1:3,:)

! Reduce the number of such points by symmetrization.
  ABI_MALLOC(ibz2bz,(nqbz))
  ABI_MALLOC(wtq,(nqbz))
  ABI_MALLOC(bz2ibz_smap, (6, nqbz))
  ABI_MALLOC(wtq_folded,(nqbz))
  wtq(:)=one/nqbz         ! Weights sum up to one

!FB  write(InVar%stdlog,*) 'nqbz = ', nqbz
  call symkpt(0,Crystal%gmet,ibz2bz,InVar%stdlog,qbz,nqbz,nqibz,Sym%nsym,Crystal%symrec,1,wtq,wtq_folded, &
    bz2ibz_smap, xmpi_comm_self)
!FB  write(InVar%stdlog,*) 'nqibz = ', nqibz

  ABI_FREE(bz2ibz_smap)

  ABI_MALLOC(wtqibz   ,(nqibz))
  ABI_MALLOC(qibz     ,(3,nqibz))
  ABI_MALLOC(qibz_cart,(3,nqibz))
  do iq_ibz=1,nqibz
    wtqibz(iq_ibz)=wtq_folded(ibz2bz(iq_ibz))
    qibz(:,iq_ibz)=qbz(:,ibz2bz(iq_ibz))
  end do
  ABI_FREE(wtq_folded)

! Loop over irreducible q-points
! =======================
  ntemp=1000
  nmode=3*InVar%natom_unitcell
  k_B=kb_HaK*Ha_eV
  ABI_MALLOC(heatcapa,(nmode))   ; heatcapa   (:)=0.d0
  ABI_MALLOC(grun_thermo,(nmode)); grun_thermo(:)=0.d0
  ABI_MALLOC(p_thermo1,(nmode))  ; p_thermo1  (:)=0.d0
  ABI_MALLOC(p_thermo2,(ntemp))  ; p_thermo2  (:)=0.d0
  ABI_MALLOC(u_vib,(nmode))      ; u_vib      (:)=0.d0
  ABI_MALLOC(heatcapa_HA,(nmode,ntemp))   ; heatcapa_HA   (:,:)=0.d0
  ABI_MALLOC(grun_thermo_HA,(nmode,ntemp)); grun_thermo_HA(:,:)=0.d0
  ABI_MALLOC(p_thermo_HA,(nmode,ntemp))   ; p_thermo_HA   (:,:)=0.d0
  ABI_MALLOC(u_vib_HA,(nmode,ntemp))      ; u_vib_HA      (:,:)=0.d0
  ABI_MALLOC(displ,(2*3*InVar%natom_unitcell*3*InVar%natom_unitcell,nqibz)); displ(:,:)=zero
  ABI_MALLOC(out_eigvec,(2,3,InVar%natom_unitcell,3,InVar%natom_unitcell,nqibz)); out_eigvec(:,:,:,:,:,:)=zero
  ABI_MALLOC(omega,(3*InVar%natom_unitcell,nqibz)); omega(:,:)=zero
  call tdep_init_eigen2nd(Eigen2nd_tmp,InVar%natom_unitcell,nqibz)
  do iq_ibz=1,nqibz
!   Compute the frequencies
!   =======================
    call ifc%fourq(Crystal,qibz(:,iq_ibz),omega(:,iq_ibz),displ(:,iq_ibz),out_eigvec=out_eigvec(:,:,:,:,:,iq_ibz))
    qibz_cart(:,iq_ibz)=matmul(Crystal%gprimd,qibz(:,iq_ibz))

!   Compute the gruneisen for this q-point
!   ======================================
    Eigen2nd_tmp%eigenval(:,iq_ibz)  =omega(:,iq_ibz)
    do iatcell=1,InVar%natom_unitcell
      do jatcell=1,InVar%natom_unitcell
        do ii=1,3
          do jj=1,3
            Eigen2nd_tmp%eigenvec((iatcell-1)*3+ii,(jatcell-1)*3+jj,iq_ibz)=dcmplx(out_eigvec(1,ii,iatcell,jj,jatcell,iq_ibz),&
&                                                                                  out_eigvec(2,ii,iatcell,jj,jatcell,iq_ibz))
          end do
        end do
      end do
    end do
    call tdep_calc_gruneisen(distance,Eigen2nd_tmp,Gruneisen,iq_ibz,InVar,Lattice,Psij_ref,qibz_cart(:,iq_ibz),&
&     Rlatt_cart,Shell3at,Sym)

!   Compute the heat capacity and thermodynamical gruneisen parameter at present temperature
!   ========================================================================================
    wovert=1.d0/(2*InVar%temperature*k_B)
    do ii=1,nmode
      xx=omega(ii,iq_ibz)*Ha_eV
      if (xx.le.0) cycle
      heatcapa(ii)   =heatcapa(ii)   +(wovert*xx/sinh(wovert*xx))**2*wtqibz(iq_ibz)
      grun_thermo(ii)=grun_thermo(ii)+(wovert*xx/sinh(wovert*xx))**2*wtqibz(iq_ibz)*Gruneisen(ii) ! Gamma= sum_i Gamma_i*C_Vi / C_V
      p_thermo1(ii)  =p_thermo1(ii)  +  (xx/2.d0/tanh(wovert*xx))   *wtqibz(iq_ibz)*Gruneisen(ii) ! P    = sum_i Gamma_i*E_i / V
      u_vib(ii)      =u_vib(ii)      +  (xx/2.d0/tanh(wovert*xx))   *wtqibz(iq_ibz)
    end do
!   Compute the heat capacity and thermodynamical gruneisen parameter as a function of temperature
!   ==============================================================================================
    do itemp=1,ntemp
      wovert=1.d0/(2*real(itemp)*10*k_B)
      do ii=1,nmode
        xx=omega(ii,iq_ibz)*Ha_eV
        if (xx.le.0) cycle
        heatcapa_HA(ii,itemp)   =heatcapa_HA(ii,itemp)   +(wovert*xx/sinh(wovert*xx))**2*wtqibz(iq_ibz)
        grun_thermo_HA(ii,itemp)=grun_thermo_HA(ii,itemp)+(wovert*xx/sinh(wovert*xx))**2*wtqibz(iq_ibz)*Gruneisen(ii)
        p_thermo_HA(ii,itemp)   =p_thermo_HA(ii,itemp)   +  (xx/2.d0/tanh(wovert*xx))   *wtqibz(iq_ibz)*Gruneisen(ii)
        u_vib_HA(ii,itemp)      =u_vib_HA(ii,itemp)      +  (xx/2.d0/tanh(wovert*xx))   *wtqibz(iq_ibz)
      end do
    end do
  end do
! Compute the pressure as the integral of Gamma*C_v overt T
  ABI_MALLOC(tmp, (ntemp))
  tmp(:)=0.d0
  do itemp=1,ntemp
!FB    tmp(itemp)=sum(heatcapa_HA(:,itemp))
    tmp(itemp)=sum(grun_thermo_HA(:,itemp))
  end do
  call simpson_int(ntemp,10.d0,tmp,p_thermo2)
  ABI_FREE(tmp)

  open(unit=20,file=trim(InVar%output_prefix)//'thermo3.dat')
  write(20,'(a)')'#   T(K)    C_v(k_B/fu)        Gamma     alpha_v*10^6(K^-1)   E_th(eV)                       P_th_(GPa)'
  write(20,'(a,72x,a)')'#',' ----------------------------------------------'
  write(20,'(a,72x,a)')'#','  {sum G_i.U_iV}  {int G.C_v/V dT}    {G.U/V}'
  do itemp=1,ntemp
    C_v    =sum(heatcapa_HA   (:,itemp))
    E_th   =sum(u_vib_HA    (:,itemp))
    Gama   =sum(grun_thermo_HA(:,itemp))/C_v
    alpha_v=Gama*C_v*kb_HaK*Ha_J/(Lattice%BulkModulus*1.d9)/(Lattice%ucvol*Bohr_Ang**3*1.d-30)
    P_th1  =sum(p_thermo_HA (:,itemp))/(Lattice%ucvol*Bohr_Ang**3*1.d-30)*e_Cb/10**9
    P_th2  =      p_thermo2(itemp)*k_B/(Lattice%ucvol*Bohr_Ang**3*1.d-30)*e_Cb/10**9
    P_th   =                 Gama*E_th/(Lattice%ucvol*Bohr_Ang**3*1.d-30)*e_Cb/10**9
    write(20,'(1x,i5,7(1x,f15.5))') itemp*10,C_v,Gama,alpha_v*1.d6,E_th,P_th1,P_th2,P_th
  end do
  close(20)
  ABI_FREE(heatcapa_HA)
  ABI_FREE(grun_thermo_HA)
  ABI_FREE(p_thermo_HA)
  ABI_FREE(u_vib_HA)

  call tdep_destroy_eigen2nd(Eigen2nd_tmp)
  ABI_FREE(displ)
  ABI_FREE(out_eigvec)
  ABI_FREE(omega)
  ABI_FREE(ibz2bz)
  ABI_FREE(wtq)
  ABI_FREE(qbz)
  ABI_FREE(wtqibz)
  ABI_FREE(qibz)
  ABI_FREE(qibz_cart)

! Summary
! =======
  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '####### Gruneisen parameter, Thermal expansion, Thermal pressure... #########'
  write(InVar%stdout,*) '#############################################################################'
!FB  heatcapa   (:)=heatcapa   (:)/InVar%natom_unitcell
!FB  write(InVar%stdlog,*) 'Summary:'
!FB  write(InVar%stdlog,'(a,1x,100(f15.6,1x))') 'C_v(per_mode \& per unit cell)=',(heatcapa(ii)   ,ii=1,nmode)
!FB  write(InVar%stdlog,'(a,1x,100(f15.6,1x))') 'Gamma(per mode)=',(grun_thermo(ii),ii=1,nmode)
  C_v    =sum(heatcapa(:))
  E_th   =sum(u_vib   (:))
  Gama   =sum(grun_thermo(:))/C_v
  alpha_v=Gama*C_v*kb_HaK*Ha_J/(Lattice%BulkModulus*1.d9)/(Lattice%ucvol*Bohr_Ang**3*1.d-30)
  P_th1  =                   sum(p_thermo1(:))    /(Lattice%ucvol*Bohr_Ang**3*1.d-30)*e_Cb/10**9
  P_th2  =p_thermo2(int(InVar%temperature/10))*k_B/(Lattice%ucvol*Bohr_Ang**3*1.d-30)*e_Cb/10**9
  P_th   =                           Gama*E_th    /(Lattice%ucvol*Bohr_Ang**3*1.d-30)*e_Cb/10**9
  write(InVar%stdout,'(a,1x,f15.5)') ' Specific heat (k_B/f.u.):             C_v=',C_v
  write(InVar%stdout,'(a,1x,f15.5)') ' Gruneisen parameter :               Gamma=',Gama
  write(InVar%stdout,'(a,1x,f15.5)') ' Isothermal Bulk Modulus (GPa):        B_T=',Lattice%BulkModulus
  write(InVar%stdout,'(a,1x,f15.5)') ' Volume (bohr^3 per unit cell):          V=',Lattice%ucvol
  write(InVar%stdout,'(a,1x,f15.5)') ' Thermal energy (eV) :                E_th=',E_th
  write(InVar%stdout,'(a,1x,f15.5)') ' Thermal expansion (K^{-1}*10^6) : alpha_v=',alpha_v*1.d6
  write(InVar%stdout,'(a)') ' Thermal pressure (in GPa) : '
  write(InVar%stdout,'(a,1x,f15.5)') '  - with    intrinsic effects and with    ZPE : P_th=sum_i Gamma_i*E_i/V   =',P_th1
  write(InVar%stdout,'(a,1x,f15.5)') '  - with    intrinsic effects and without ZPE : P_th=integ{Gamma*C_v/V dT} =',P_th2
  write(InVar%stdout,'(a,1x,f15.5)') '  - without intrinsic effects and without ZPE : P_th=Gamma*E_th/V          =',P_th
  ABI_FREE(heatcapa)
  ABI_FREE(grun_thermo)
  ABI_FREE(p_thermo1)
  ABI_FREE(p_thermo2)
  ABI_FREE(u_vib)

!FB  0.473294364993209*(2.38255605878933*1.38065e-23)/1.20512e11/3.0/((30.6135754000*0.529177e-10)**3/216)

end subroutine tdep_calc_alpha_gamma

!=====================================================================================================
subroutine tdep_write_gruneisen(distance,Eigen2nd,InVar,Lattice,Psij_ref,Qpt,Rlatt_cart,Shell3at,Sym)

  type(Symetries_Variables_type),intent(in) :: Sym
  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(in) :: Shell3at
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Eigen_Variables_type),intent(in) :: Eigen2nd
  type(Qpoints_type),intent(in) :: Qpt
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision,intent(in) :: Psij_ref(3,3,3,Shell3at%nshell)
  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)

  integer :: iqpt,nmode,ii !,jj
  double precision :: qpt_cart(3)
  double complex, allocatable :: Gruneisen(:)

  nmode=3*InVar%natom_unitcell
  ABI_MALLOC(Gruneisen,(3*InVar%natom_unitcell)); Gruneisen(:)=zero
  open(unit=53,file=trim(InVar%output_prefix)//'gruneisen.dat')
  do iqpt=1,Qpt%nqpt
    qpt_cart(:)=Qpt%qpt_cart(:,iqpt)
    if ((sum(abs(Qpt%qpt_red(:,iqpt)))).lt.tol8) cycle  ! G point
    if (abs(sum(Qpt%qpt_red(:,iqpt)**2)-1.d0).lt.tol8) cycle ! Gp point
    call tdep_calc_gruneisen(distance,Eigen2nd,Gruneisen,iqpt,InVar,Lattice,Psij_ref,qpt_cart,Rlatt_cart,Shell3at,Sym)
!   Write the Gruneisen
    if (sum(abs(dimag(Gruneisen(:)))).gt.tol8) then
      write(53,'(i5,1x,100(e15.6,1x))') iqpt,(real(Gruneisen(ii)),ii=1,nmode),(dimag(Gruneisen(ii)),ii=1,nmode)
      ABI_BUG('The imaginary part of the Gruneisen is not equal to zero')
    else
!FB      write(53,'(i5,1x,500(e15.6,1x))') iqpt,(real(Gruneisen(ii)),ii=1,nmode),((real(Grun_shell(ii,jj)),ii=1,nmode),jj=1,Shell3at%nshell)
      write(53,'(i5,1x,500(e15.6,1x))') iqpt,(real(Gruneisen(ii)),ii=1,nmode)
    end if
  end do
  close(53)
  ABI_FREE(Gruneisen)

end subroutine tdep_write_gruneisen

!=====================================================================================================
end module m_tdep_psij
