
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_phi3

  use defs_basis
  use m_errors
  use m_abicore
  use m_numeric_tools
  use m_linalg_interfaces
  use m_xmpi
  use m_io_tools
  use m_crystal,          only : crystal_t
  use m_ddb,              only : ddb_type
  use m_ifc,              only : ifc_type
  use m_tdep_abitypes,    only : Qbz_type
  use m_htetra
  use m_kpts,             only : kpts_ibz_from_kptrlatt, tetra_from_kptrlatt
  use m_tdep_qpt,         only : Qpoints_type
  use m_tdep_readwrite,   only : Input_Variables_type, MPI_enreg_type
  use m_tdep_latt,        only : Lattice_Variables_type
  use m_tdep_shell,       only : Shell_Variables_type
  use m_tdep_sym,         only : Symetries_Variables_type
  use m_tdep_phi2,        only : Eigen_Variables_type,tdep_init_eigen2nd,tdep_destroy_eigen2nd
  use m_tdep_utils,       only : Coeff_Moore_type, Constraints_Variables_type

  implicit none

  public :: tdep_calc_ftot3
  public :: tdep_calc_phi3fcoeff
  public :: tdep_calc_phi3ref
  public :: tdep_write_phi3
  public :: tdep_build_phi3_333
  public :: tdep_calc_gruneisen
  public :: tdep_write_gruneisen
  public :: tdep_calc_alpha_gamma
  public :: tdep_calc_lifetime1
!  public :: tdep_calc_lifetime2

contains

!====================================================================================================
 subroutine tdep_calc_ftot3(Forces_TDEP,InVar,MPIdata,Phi3_ref,Phi3UiUjUk,Shell3at,ucart,Sym) 

  implicit none 

  type(Input_Variables_type),intent(in) :: InVar
  type(MPI_enreg_type), intent(in) :: MPIdata
  type(Shell_Variables_type),intent(in) :: Shell3at
  type(Symetries_Variables_type),intent(in) :: Sym
  double precision, intent(in)  :: ucart(3,InVar%natom,InVar%my_nstep)
  double precision, intent(in)  :: Phi3_ref(3,3,3,Shell3at%nshell)
  double precision, intent(out) :: Phi3UiUjUk(InVar%my_nstep)
  double precision, intent(inout) :: Forces_TDEP(3*InVar%natom*InVar%my_nstep)

  integer :: iatom,jatom,katom,isym,itrans,ishell,iatshell
  integer :: ii,jj,kk,istep
  double precision, allocatable :: Phi3_333(:,:,:)
  double precision, allocatable :: ucart_blas(:)
  double precision, allocatable :: ftot3(:,:)

  ABI_MALLOC(Phi3_333,(3,3,3)) ; Phi3_333(:,:,:)=0.d0
  ABI_MALLOC(ftot3,(3*InVar%natom,InVar%my_nstep)); ftot3(:,:)=0.d0
  do iatom=1,InVar%natom
    do ishell=1,Shell3at%nshell
!     Build the 3x3x3 IFC of an atom in this shell    
      if (Shell3at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
      do iatshell=1,Shell3at%neighbours(iatom,ishell)%n_interactions
        jatom=Shell3at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
        katom=Shell3at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
        isym =Shell3at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
        itrans=Shell3at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
        call tdep_build_phi3_333(isym,InVar,Phi3_ref(:,:,:,ishell),Phi3_333,Sym,itrans) 
!       Calculation of the force components (third order)
        do istep=1,InVar%my_nstep
          do ii=1,3
            do jj=1,3
              do kk=1,3
                ftot3(3*(iatom-1)+ii,istep)=ftot3(3*(iatom-1)+ii,istep)+&
&                    Phi3_333(ii,jj,kk)*ucart(jj,jatom,istep)*ucart(kk,katom,istep)
              end do !kk 
            end do !jj 
          end do !ii 
        end do !istep
      end do !iatshell
    end do !ishell  
  end do !iatom  
  ABI_FREE(Phi3_333)
  ftot3(:,:)=ftot3(:,:)/2.d0

  ABI_MALLOC(ucart_blas,(3*InVar%natom)); ucart_blas(:)=0.d0
  do istep=1,InVar%my_nstep
    ucart_blas(:)=0.d0 
    do jatom=1,InVar%natom
      do jj=1,3
        ucart_blas(3*(jatom-1)+jj)=ucart(jj,jatom,istep)
      end do
    end do
    call DGEMM('T','N',1,1,3*InVar%natom,1./3.d0,ftot3(:,istep),3*InVar%natom,ucart_blas,&
&              3*InVar%natom,0.d0,Phi3UiUjUk(istep),3*InVar%natom)
    Forces_TDEP(3*InVar%natom*(istep-1)+1:3*InVar%natom*istep)=&
&   Forces_TDEP(3*InVar%natom*(istep-1)+1:3*InVar%natom*istep)-ftot3(:,istep)
  end do
  ABI_FREE(ucart_blas)
  ABI_FREE(ftot3)

 end subroutine tdep_calc_ftot3

!====================================================================================================
subroutine tdep_calc_phi3fcoeff(CoeffMoore,InVar,proj,Shell3at,Sym,ucart)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell3at
  type(Coeff_Moore_type), intent(inout) :: CoeffMoore
  double precision, intent(in) :: ucart(3,InVar%natom,InVar%my_nstep)
  double precision, intent(in) :: proj(27,27,Shell3at%nshell)

  integer :: ishell,ncoeff,ncoeff_prev,istep,iatom,jatom,katom
  integer :: icoeff,isym,itrans,iatshell
  integer :: mu,nu,xi,alpha,beta,gama,iindex
  double precision :: temp
  double precision :: udiff_ki(3),udiff_ji(3)
  double precision, allocatable :: SSS_proj(:,:,:,:)
  type(Constraints_Variables_type) :: Const

  ABI_MALLOC(Const%Sprod,(Sym%nsym,6))
  do isym=1,Sym%nsym
    do itrans=1,6
      ABI_MALLOC(Const%Sprod(isym,itrans)%SSS,(3,27,3,3)); Const%Sprod(isym,itrans)%SSS(:,:,:,:)=zero
    end do  
  end do  

! For each couple of atoms, transform the Phi3 (3x3x3) ifc matrix using the symetry operation (S)
! Note: iatom=1 is excluded in order to take into account the atomic sum rule (see below)
  do isym=1,Sym%nsym
    do mu=1,3
      do alpha=1,3
        do nu=1,3
          do beta=1,3
            do xi=1,3
              do gama=1,3
                temp=Sym%S_ref(mu,alpha,isym,1)*Sym%S_ref(nu,beta,isym,1)*Sym%S_ref(xi,gama,isym,1)
                Const%Sprod(isym,1)%SSS(mu,gama+(beta-1)*3+(alpha-1)*9,nu,xi)=temp !\Psi_efg
                Const%Sprod(isym,2)%SSS(mu,gama+(beta-1)*3+(alpha-1)*9,xi,nu)=temp !\Psi_egf
                Const%Sprod(isym,3)%SSS(nu,gama+(beta-1)*3+(alpha-1)*9,mu,xi)=temp !\Psi_feg
                Const%Sprod(isym,4)%SSS(nu,gama+(beta-1)*3+(alpha-1)*9,xi,mu)=temp !\Psi_fge
                Const%Sprod(isym,5)%SSS(xi,gama+(beta-1)*3+(alpha-1)*9,mu,nu)=temp !\Psi_gef
                Const%Sprod(isym,6)%SSS(xi,gama+(beta-1)*3+(alpha-1)*9,nu,mu)=temp !\Psi_gfe
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
!FB        if (iatom==jatom.or.iatom==katom) cycle
        isym =Shell3at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
        itrans=Shell3at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
        ncoeff     =Shell3at%ncoeff(ishell)
        ncoeff_prev=Shell3at%ncoeff_prev(ishell)+CoeffMoore%ncoeff2nd+CoeffMoore%ncoeff1st
  
        ABI_MALLOC(SSS_proj,(3,3,3,ncoeff)) ; SSS_proj(:,:,:,:)=zero
        do mu=1,3
          do nu=1,3
            do xi=1,3
              do icoeff=1,ncoeff
                SSS_proj(mu,nu,xi,icoeff)=DDOT(27,Const%Sprod(isym,itrans)%SSS(mu,:,nu,xi),1,proj(:,icoeff,ishell),1)
              end do
            end do
          end do
        end do
        do istep=1,InVar%my_nstep
          iindex=3*(iatom-1)+3*InVar%natom*(istep-1)
!         In order to impose the acoustic sum rule we use : 
!FB          udiff_ji(:)=ucart(:,jatom,istep)-ucart(:,iatom,istep)
!FB          udiff_ki(:)=ucart(:,katom,istep)-ucart(:,iatom,istep)
          udiff_ji(:)=ucart(:,jatom,istep)
          udiff_ki(:)=ucart(:,katom,istep)
!         F_i^{\mu}(t)=\sum_{\alpha\beta\gamma,jk,\nu\xi} S^{\mu\alpha}.S^{\nu\beta}.S^{\xi\gamma}.\Psi_{ijk}^{\alpha\beta\gamma}.udiff_k^\xi(t).udiff_j^\nu(t)
          do nu=1,3
            do xi=1,3
              CoeffMoore%fcoeff(iindex+1:iindex+3,ncoeff_prev+1:ncoeff_prev+ncoeff)= &
&             CoeffMoore%fcoeff(iindex+1:iindex+3,ncoeff_prev+1:ncoeff_prev+ncoeff)+&
&             SSS_proj(1:3,nu,xi,1:ncoeff)*udiff_ji(nu)*udiff_ki(xi)/2.d0
            end do  
          end do  
        end do !istep
        ABI_FREE(SSS_proj)
      end do !iatshell
    end do !iatom
  end do !ishell
  write(InVar%stdout,*) ' ------- achieved'
  do isym=1,Sym%nsym
    do itrans=1,6
      ABI_FREE(Const%Sprod(isym,itrans)%SSS)
    end do
  end do
  ABI_FREE(Const%Sprod)

end subroutine tdep_calc_phi3fcoeff

!=====================================================================================================
subroutine tdep_calc_phi3ref(InVar,ntotcoeff,proj,Phi3_coeff,Phi3_ref,Shell3at)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(in) :: Shell3at
  integer,intent(in) :: ntotcoeff
  double precision, intent(in) :: proj(27,27,Shell3at%nshell)
  double precision, intent(in) :: Phi3_coeff(ntotcoeff,1)
  double precision, intent(inout) :: Phi3_ref(3,3,3,Shell3at%nshell)

  integer :: ishell,ncoeff,ncoeff_prev
  integer :: ii,jj,kk,kappa

  do ishell=1,Shell3at%nshell
!   Build the 3x3x3 IFC per shell    
    ncoeff     =Shell3at%ncoeff(ishell)
    ncoeff_prev=Shell3at%ncoeff_prev(ishell)
    kappa=0  
    do ii=1,3
      do jj=1,3
        do kk=1,3
          kappa=kappa+1
          Phi3_ref(ii,jj,kk,ishell)=sum(proj(kappa,1:ncoeff,ishell)*Phi3_coeff(ncoeff_prev+1:ncoeff_prev+ncoeff,1))
        end do  
      end do  
    end do  
!   Remove the rounding errors before writing (for non regression testing purposes)
    do ii=1,3
      do jj=1,3
        do kk=1,3
          if (abs(Phi3_ref(ii,jj,kk,ishell)).lt.tol8) Phi3_ref(ii,jj,kk,ishell)=zero
        end do
      end do
    end do  
  end do  

end subroutine tdep_calc_phi3ref

!=====================================================================================================
subroutine tdep_write_phi3(distance,InVar,Phi3_ref,Shell3at,Sym)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell3at
  double precision, intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision, intent(in) :: Phi3_ref(3,3,3,Shell3at%nshell)

  integer :: ishell,jshell,isym,jatom,katom
  integer :: iatref,jatref,katref,iatshell,itrans
  integer :: ii,jj
  double precision, allocatable :: Phi3_333(:,:,:)

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '#### For each shell, list of coefficients (IFC), number of neighbours... ####'
  write(InVar%stdout,*) '#############################################################################'

! Write the IFCs in the data.out file (with others specifications: 
! number of atoms in a shell, Trace...)
  ABI_MALLOC(Phi3_333,(3,3,3)) ; Phi3_333(:,:,:)=0.d0
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
        itrans=Shell3at%neighbours(iatref,ishell)%transpose_in_shell(iatshell)
        call tdep_build_phi3_333(isym,InVar,Phi3_ref(:,:,:,ishell),Phi3_333,Sym,itrans) 
        write(InVar%stdout,'(a,i4,a,i4)') '  For iatcell=',iatref,' ,with type=',mod(iatref-1,InVar%natom_unitcell)+1
        write(InVar%stdout,'(a,i4,a,i4)') '  For jatom  =',jatom ,' ,with type=',mod(jatom -1,InVar%natom_unitcell)+1
        write(InVar%stdout,'(a,i4,a,i4)') '  For katom  =',katom ,' ,with type=',mod(katom -1,InVar%natom_unitcell)+1
        write(InVar%stdout,'(a)') '  \Psi^{\alpha\beta x}='
        do ii=1,3
          write(InVar%stdout,'(2x,3(f9.6,1x))') (Phi3_333(ii,jj,1),jj=1,3)
        end do
        write(InVar%stdout,'(a)') '  \Psi^{\alpha\beta y}='
        do ii=1,3
          write(InVar%stdout,'(2x,3(f9.6,1x))') (Phi3_333(ii,jj,2),jj=1,3)
        end do
        write(InVar%stdout,'(a)') '  \Psi^{\alpha\beta z}='
        do ii=1,3
          write(InVar%stdout,'(2x,3(f9.6,1x))') (Phi3_333(ii,jj,3),jj=1,3)
        end do
        write(InVar%stdout,'(a,3(f9.6,1x))') '  (i,j) vector components:', (distance(iatref,jatom,jj+1),jj=1,3)
        write(InVar%stdout,'(a,3(f9.6,1x))') '  (j,k) vector components:', (distance(jatom ,katom,jj+1),jj=1,3)
        write(InVar%stdout,'(a,3(f9.6,1x))') '  (k,i) vector components:', (distance(katom,iatref,jj+1),jj=1,3)
        write(InVar%stdout,*) ' '
      end do !iatshell
    end if !n_interactions 
  end do !ishell 
  ABI_FREE(Phi3_333)

end subroutine tdep_write_phi3

!=====================================================================================================
subroutine tdep_calc_gruneisen(distance,Eigen2nd,Gruneisen,iqpt,InVar,Lattice,Phi3_ref,qpt_cart,Rlatt_cart,Shell3at,Sym)

  implicit none

  type(Symetries_Variables_type),intent(in) :: Sym
  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(in) :: Shell3at
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Eigen_Variables_type),intent(in) :: Eigen2nd
  integer,intent(in) :: iqpt
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision,intent(in) :: Phi3_ref(3,3,3,Shell3at%nshell)
  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)
  double precision,intent(in) :: qpt_cart(3)
  double complex  ,intent(out) :: Gruneisen(3*InVar%natom_unitcell,3,3)

  integer :: ii,jj,kk,iatcell,jatom,katom,itypat,jtypat,natom_unitcell,ll
  integer :: imode,nmode,jat_mod,jatcell,iatshell,ishell,isym,itrans
  double precision :: phase
  double precision :: F_strain(3,3)
  double precision, allocatable :: Phi3_333(:,:,:)
  double complex, allocatable :: mass_mat(:,:),eigen_prod(:,:,:,:,:),eigenvec_loc(:,:),omega(:)

! Define quantities
  natom_unitcell=InVar%natom_unitcell
  nmode=3*InVar%natom_unitcell
  do ii=1,3
    do jj=1,3
      if (ii.eq.jj) then
!FB        F_strain(ii,jj)= 1./dsqrt(6.d0)
!FB        F_strain(ii,jj)=1./dsqrt(15.d0)
        F_strain(ii,jj)= 1.d0
!FB        F_strain(ii,jj)= 0.8
      else 
!FB        F_strain(ii,jj)= 1./dsqrt(6.d0)/2.d0
!FB        F_strain(ii,jj)=2./dsqrt(15.d0)
        F_strain(ii,jj)=0.d0
!FB        F_strain(ii,jj)= 0.3
      end if
    end do
  end do

! Calculation of the eigenvectors product (and phonon frequencies)
  ABI_MALLOC(eigen_prod,  (3,3,natom_unitcell,natom_unitcell,nmode)) ; eigen_prod(:,:,:,:,:)=czero
  ABI_MALLOC(eigenvec_loc,(3*natom_unitcell,3*natom_unitcell))       ; eigenvec_loc(:,:)    =czero
  do iatcell=1,natom_unitcell
    do jatcell=1,natom_unitcell
      do ii=1,3
        do jj=1,3
          eigenvec_loc((iatcell-1)*3+ii,(jatcell-1)*3+jj)=dcmplx(Eigen2nd%eigenvec(1,ii,iatcell,jj,jatcell,iqpt),&
&                                                                Eigen2nd%eigenvec(2,ii,iatcell,jj,jatcell,iqpt))
        end do !ii 
      end do !jj 
    end do !jatcell 
  end do !iatcell
  ABI_MALLOC(omega,(3*natom_unitcell)); omega(:)=czero
  do iatcell=1,natom_unitcell
    do jatcell=1,natom_unitcell
      do ii=1,3
        do jj=1,3
          do imode=1,nmode
            eigen_prod(ii,jj,iatcell,jatcell,imode)=conjg(eigenvec_loc((iatcell-1)*3+ii,imode))*&
&                                                         eigenvec_loc((jatcell-1)*3+jj,imode)
            omega(imode)=omega(imode)+eigen_prod(ii,jj,iatcell,jatcell,imode)*&
&                              dcmplx(Eigen2nd%dynmat(1,ii,iatcell,jj,jatcell,iqpt),&
&                                     Eigen2nd%dynmat(2,ii,iatcell,jj,jatcell,iqpt))
!FB            write(InVar%stdlog,'(a,5(1x,i4),2(1x,f15.5))') 'imode,iatcell,jatcell,ii,jj=',imode,iatcell,jatcell,ii,jj,eigen_prod(ii,jj,iatcell,jatcell,imode)
          end do !imode 
        end do !ii 
      end do !jj 
    end do !jatcell 
  end do !iatcell

  do imode=1,nmode
    if (abs(aimag(omega(imode))).gt.tol8) then
      write(InVar%stdlog,'(a,1x,e15.8,1x,a,i4)') '>>>> WARNING : Imaginary part of the phonon frequency is not zero (',aimag(omega(imode)),') for mode :',imode 
!FB      stop
    end if
  end do  
!FB  write(InVar%stdlog,'(a,1(1x,i4),6(1x,f15.5))') 'iqpt,Grun=',iqpt,(dsqrt(real(omega(imode))),imode=1,6)
  ABI_FREE(omega)


! Calculation of the Gruneisen
  ABI_MALLOC(Phi3_333,(3,3,3)) ; Phi3_333(:,:,:)=0.d0
  Gruneisen(:,:,:)=czero
  do ishell=1,Shell3at%nshell
    do iatcell=1,InVar%natom_unitcell
      itypat=InVar%typat_unitcell(iatcell)
      do iatshell=1,Shell3at%neighbours(iatcell,ishell)%n_interactions
        jatom=Shell3at%neighbours(iatcell,ishell)%atomj_in_shell(iatshell)
        katom=Shell3at%neighbours(iatcell,ishell)%atomk_in_shell(iatshell)
        isym =Shell3at%neighbours(iatcell,ishell)%sym_in_shell(iatshell)
        itrans=Shell3at%neighbours(iatcell,ishell)%transpose_in_shell(iatshell)
        jat_mod=mod(jatom+InVar%natom_unitcell-1,InVar%natom_unitcell)+1
        jtypat=InVar%typat_unitcell(jat_mod)
        phase=0.d0
        do jj=1,3
          phase=phase+2*pi*Rlatt_cart(jj,iatcell,jatom)*qpt_cart(jj)
        end do
        call tdep_build_phi3_333(isym,InVar,Phi3_ref(:,:,:,ishell),Phi3_333,Sym,itrans) 
        do ii=1,3
          do jj=1,3
            do kk=1,3
             do ll=1,3
                do imode=1,nmode
                  Gruneisen(imode,kk,ll)=Gruneisen(imode,kk,ll)-dcmplx(Phi3_333(ii,jj,kk),0.d0)*&
&                                  eigen_prod(ii,jj,iatcell,jat_mod,imode)*&
&                                  exp(dcmplx(0.d0,phase))*dcmplx(F_strain(kk,ll),0.d0)*dcmplx(distance(iatcell,katom,ll+1),0.d0)/&
&                                  dcmplx(dsqrt(InVar%amu(itypat)*InVar%amu(jtypat))*amu_emass,0.d0)/&
&                                  dcmplx(6*Eigen2nd%eigenval(imode,iqpt)**2,0.d0)
!FB                  do ll=1,3
!FB                    Gruneisen(imode,kk,ll)=Gruneisen(imode,kk,ll)-dcmplx(Phi3_333(ii,jj,kk),0.d0)*&
!FB&                                    eigen_prod(ii,jj,iatcell,jat_mod,imode)*&
!FB&                                    exp(dcmplx(0.d0,phase))*F_strain(kk,ll)*dcmplx(distance(iatcell,katom,ll+1),0.d0)/&
!FB&                                    dcmplx(dsqrt(InVar%amu(itypat)*InVar%amu(jtypat))*amu_emass,0.d0)/&
!FB&                                    dcmplx(6*Eigen2nd%eigenval(imode,iqpt)**2,0.d0)
!FB!FB                  write(InVar%stdlog,'(a,4(1x,i4),2(1x,f15.5))') 'ii,jj,kk,imode,Grun=',ii,jj,kk,imode,Gruneisen(imode,kk,kk)
!FB                  end do !ll
                end do !imode
              end do !ll
            end do !kk
          end do !jj
        end do !ii 
      end do !iatshell 
    end do !iatcell 
  end do !ishell

!FB  Gruneisen(:,:,:)=Gruneisen(:,:,:)*dcmplx(dsqrt(6.d0),0.d0)
!FB  do imode=1,nmode
!FB    write(InVar%stdlog,'(a,1(1x,i4),6(1x,f15.5))') 'imode,Grun=',imode,(Gruneisen(imode,kk,kk),kk=1,3)
!FB  end do  
  ABI_FREE(Phi3_333)
  ABI_FREE(eigen_prod)
  ABI_FREE(eigenvec_loc)

end subroutine tdep_calc_gruneisen

!=====================================================================================================
subroutine tdep_build_phi3_333(isym,InVar,Phi3_ref,Phi3_333,Sym,itrans) 

  implicit none

  type(Symetries_Variables_type),intent(in) :: Sym
  type(Input_Variables_type),intent(in) :: InVar
  double precision, intent(in) :: Phi3_ref(3,3,3)
  double precision, intent(out) :: Phi3_333(3,3,3)
  integer,intent(in) :: isym,itrans

  integer :: alpha,beta,gama
  integer :: ii,jj,kk,ee,ff,gg,mu,nu,xi
  double precision :: Phi3_tmp(3,3,3)


! Transform in the new basis wrt S_ref
  Phi3_333(:,:,:)=zero
  do mu=1,3
    do alpha=1,3
      do nu=1,3
        do beta=1,3
          do xi=1,3
            do gama=1,3
              Phi3_333(mu,nu,xi)=Phi3_333(mu,nu,xi)+&
&             Sym%S_ref(mu,alpha,isym,1)*Sym%S_ref(nu,beta,isym,1)*Sym%S_ref(xi,gama,isym,1)*Phi3_ref(alpha,beta,gama)
            end do
          end do
        end do
      end do
    end do
  end do
  
! Take into account the 6 allowed permutations
  Phi3_tmp(:,:,:)=Phi3_333(:,:,:)
  if ((itrans.lt.1).or.(itrans.gt.6)) then  
    MSG_BUG('This value of the symmetry index is not permitted')
  end if
  do ii=1,3
    do jj=1,3
      do kk=1,3
        if (itrans==1) then ; ee=ii ; ff=jj ; gg=kk ; endif !\Psi_efg
        if (itrans==2) then ; ee=ii ; ff=kk ; gg=jj ; endif !\Psi_egf
        if (itrans==3) then ; ee=jj ; ff=ii ; gg=kk ; endif !\Psi_feg
        if (itrans==4) then ; ee=jj ; ff=kk ; gg=ii ; endif !\Psi_fge
        if (itrans==5) then ; ee=kk ; ff=ii ; gg=jj ; endif !\Psi_gef
        if (itrans==6) then ; ee=kk ; ff=jj ; gg=ii ; endif !\Psi_gfe
        Phi3_333(ee,ff,gg)=Phi3_tmp(ii,jj,kk)
      end do
    end do  
  end do          

end subroutine tdep_build_phi3_333

!=====================================================================================================
subroutine tdep_calc_alpha_gamma(distance,DDB,Eigen2nd,InVar,Lattice,MPIdata,Phi3_ref,Qbz,Rlatt_cart,Shell3at,Sym)

  implicit none

  type(ddb_type),intent(in) :: DDB
  type(Eigen_Variables_type),intent(in) :: Eigen2nd
  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(inout) :: Lattice
  type(Shell_Variables_type),intent(in) :: Shell3at
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Qbz_type),intent(in) :: Qbz
  type(MPI_enreg_type), intent(in) :: MPIdata
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision,intent(in) :: Phi3_ref(3,3,3,Shell3at%nshell)
  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)

  integer :: iq_ibz,ii,iatcell,jatcell,jj,nmode,ntemp,itemp,kk,imode,ll
  integer :: alpha,beta
  double precision :: k_B,wovert,xx,C_v,Gama,alpha_v,P_th,E_th,P_th1,P_th2,Vp,Vs,tmp1
  double precision, allocatable :: heatcapa(:),grun_thermo(:,:,:),u_vib(:)
  double precision, allocatable :: heatcapa_HA(:,:),grun_thermo_HA(:,:,:,:),p_thermo1(:,:,:)
  double precision, allocatable :: p_thermo_HA(:,:,:,:),u_vib_HA(:,:)
  double precision, allocatable :: p_thermo2(:),tmp(:)
  double precision :: Gama_tensor(3,3),alpha_v_tensor(3,3),P_th_tensor(3,3),P_th1_tensor(3,3)
  double complex, allocatable :: Gruneisen(:,:,:),Grun_mean(:)
  character(len=500) :: message
  integer :: ierr

  ierr = 0

  ABI_MALLOC(Gruneisen,(3*InVar%natom_unitcell,3,3)); Gruneisen(:,:,:)=czero
  ABI_MALLOC(Grun_mean,(3*InVar%natom_unitcell))    ; Grun_mean(:)    =czero

  ntemp=1000
  nmode=3*InVar%natom_unitcell
  k_B=kb_HaK*Ha_eV
  ABI_MALLOC(heatcapa,(nmode))       ; heatcapa   (:)    =0.d0 
  ABI_MALLOC(grun_thermo,(nmode,3,3)); grun_thermo(:,:,:)=0.d0 
  ABI_MALLOC(p_thermo1,(nmode,3,3))  ; p_thermo1  (:,:,:)=0.d0 
  ABI_MALLOC(p_thermo2,(ntemp))      ; p_thermo2  (:)    =0.d0 
  ABI_MALLOC(u_vib,(nmode))          ; u_vib      (:)    =0.d0 
  ABI_MALLOC(heatcapa_HA,(nmode,ntemp))       ; heatcapa_HA   (:,:)    =0.d0 
  ABI_MALLOC(grun_thermo_HA,(nmode,ntemp,3,3)); grun_thermo_HA(:,:,:,:)=0.d0 
  ABI_MALLOC(p_thermo_HA,(nmode,ntemp,3,3))   ; p_thermo_HA   (:,:,:,:)=0.d0 
  ABI_MALLOC(u_vib_HA,(nmode,ntemp))          ; u_vib_HA      (:,:)    =0.d0 
! Loop over irreducible q-points
! =======================
  do iq_ibz=1,Qbz%nqibz

!   Compute the gruneisen for this q-point
!   ======================================
    if ((sum(abs(Qbz%qibz_cart(:,iq_ibz)))).lt.tol8) cycle  ! G point
    call tdep_calc_gruneisen(distance,Eigen2nd,Gruneisen,iq_ibz,InVar,Lattice,&
&                            Phi3_ref,Qbz%qibz_cart(:,iq_ibz),Rlatt_cart,Shell3at,Sym)
    Grun_mean(:)=czero
    do ii=1,3*InVar%natom_unitcell
      do jj=1,3
        do kk=1,3
          Grun_mean(ii)=Grun_mean(ii)+Gruneisen(ii,jj,kk)
        end do  
      end do  
    end do
    if (sum(abs(aimag(Grun_mean(:)))).gt.tol8) then
      write(message,'(i5,1x,100(e15.6,1x))') iq_ibz,( real(Grun_mean(ii)),ii=1,nmode)
      MSG_ERROR_NOSTOP(message,ierr)
      write(message,'(i5,1x,100(e15.6,1x))') iq_ibz,(aimag(Grun_mean(ii)),ii=1,nmode)
      MSG_ERROR_NOSTOP(message,ierr)
      MSG_ERROR('tdep_calc_alpha_gamma : The imaginary part of the Gruneisen is not equal to zero')
    end if
!   If the Gruneisen is isotropic    
!FB    do ii=1,3
!FB      tmp1=0.d0
!FB      do jj=1,3
!FB        tmp1=tmp1 + real(Gruneisen(ii,jj,jj))
!FB      end do
!FB      do jj=1,3
!FB        Gruneisen(ii,jj,jj)=dcmplx(tmp1/3.d0,0.d0)
!FB      end do
!FB      write(InVar%stdlog,'(2(i5,1x),100(e15.6,1x))') iq_ibz,ii,Gruneisen(ii,jj,jj)
!FB    end do  

!   Compute the heat capacity and thermodynamical gruneisen parameter at present temperature
!   ========================================================================================
    wovert=1.d0/(2*InVar%temperature*k_B)
    do ii=1,nmode
      xx=Eigen2nd%eigenval(ii,iq_ibz)*Ha_eV
      if (xx.le.0) cycle
      heatcapa(ii)       =heatcapa(ii)       +(wovert*xx/sinh(wovert*xx))**2*Qbz%wtqibz(iq_ibz)
      grun_thermo(ii,:,:)=grun_thermo(ii,:,:)+(wovert*xx/sinh(wovert*xx))**2*Qbz%wtqibz(iq_ibz)*real(Gruneisen(ii,:,:))*3.d0 ! Gamma= sum_i Gamma_i*C_Vi / C_V
      p_thermo1(ii,:,:)  =p_thermo1(ii,:,:)  +  (xx/2.d0/tanh(wovert*xx))   *Qbz%wtqibz(iq_ibz)*real(Gruneisen(ii,:,:))*3.d0 ! P    = sum_i Gamma_i*E_i / V
      u_vib(ii)          =u_vib(ii)          +  (xx/2.d0/tanh(wovert*xx))   *Qbz%wtqibz(iq_ibz)
    end do  
!   Compute the heat capacity and thermodynamical gruneisen parameter as a function of temperature
!   ==============================================================================================
    do itemp=1,ntemp
      wovert=1.d0/(2*real(itemp)*10*k_B)
      do ii=1,nmode
        xx=Eigen2nd%eigenval(ii,iq_ibz)*Ha_eV
        if (xx.le.0) cycle
        heatcapa_HA(ii,itemp)       =heatcapa_HA(ii,itemp)       +(wovert*xx/sinh(wovert*xx))**2*Qbz%wtqibz(iq_ibz)
        grun_thermo_HA(ii,itemp,:,:)=grun_thermo_HA(ii,itemp,:,:)+(wovert*xx/sinh(wovert*xx))**2*Qbz%wtqibz(iq_ibz)*real(Gruneisen(ii,:,:))*3.d0
        p_thermo_HA(ii,itemp,:,:)   =p_thermo_HA(ii,itemp,:,:)   +  (xx/2.d0/tanh(wovert*xx))   *Qbz%wtqibz(iq_ibz)*real(Gruneisen(ii,:,:))*3.d0
        u_vib_HA(ii,itemp)          =u_vib_HA(ii,itemp)          +  (xx/2.d0/tanh(wovert*xx))   *Qbz%wtqibz(iq_ibz)
      end do  
    end do  
  end do
  ABI_FREE(Gruneisen)
  ABI_FREE(Grun_mean)
! Compute the pressure as the integral of Gamma*C_v overt T  
  allocate(tmp(ntemp)); tmp(:)=0.d0
  do itemp=1,ntemp
    do ii=1,3
      do jj=1,3
        tmp(itemp)=tmp(itemp)+sum(grun_thermo_HA(:,itemp,ii,jj))/3.d0
      end do  
    end do  
  end do  
  call simpson_int(ntemp,10.d0,tmp,p_thermo2)
  deallocate(tmp)

  if (MPIdata%iam_master) then
    open(unit=20,file=trim(InVar%output_prefix)//'thermo3.dat')
    open(unit=21,file=trim(InVar%output_prefix)//'alpha_gamma.dat')
    write(20,'(a)')'#   T(K)    C_v(k_B/fu)        Gamma     alpha_v*10^6(K^-1)   E_th(eV)                       P_th_(GPa)'
    write(20,'(a,72x,a)')'#',' ----------------------------------------------'
    write(20,'(a,72x,a)')'#','  {sum G_i.U_iV}  {int G.C_v/V dT}    {G.U/V}'
    write(21,'(a)')'#   T(K)      Gamma_11         Gamma_22        Gamma_33      alpha_11        alpha_22        alpha_33        alpha_12        alpha_13        alpha_23'
    write(20,'(a)')'# ---------------------------------------------------------------------------------------------------------------------------------------------------'
  end if  
  do itemp=1,ntemp
    C_v    =sum(heatcapa_HA   (:,itemp))
    E_th   =sum(u_vib_HA    (:,itemp))
    alpha_v_tensor(:,:)=zero
    Gama_tensor   (:,:)=zero
    P_th1_tensor  (:,:)=zero
    P_th_tensor   (:,:)=zero
    do imode=1,nmode
      Gama_tensor (:,:)=Gama_tensor (:,:)+grun_thermo_HA(imode,itemp,:,:)/C_v
      P_th1_tensor(:,:)=P_th1_tensor(:,:)+p_thermo_HA   (imode,itemp,:,:)/(Lattice%ucvol*Bohr_Ang**3*1.d-30)*e_Cb/10**9
    end do
!   alpha=sum{I=1,6;J=1,6} S(I,J) Gama_tensor (J)     
    do ii=1,3
      do jj=1,3
        if (ii.eq.1.and.jj.eq.1) alpha=1
        if (ii.eq.2.and.jj.eq.2) alpha=2
        if (ii.eq.3.and.jj.eq.3) alpha=3
        if (ii.eq.2.and.jj.eq.3) alpha=4
        if (ii.eq.1.and.jj.eq.3) alpha=5
        if (ii.eq.1.and.jj.eq.2) alpha=6
        if (ii.eq.3.and.jj.eq.2) alpha=4
        if (ii.eq.3.and.jj.eq.1) alpha=5
        if (ii.eq.2.and.jj.eq.1) alpha=6
        do kk=1,3
          do ll=1,3
            if (kk.eq.1.and.ll.eq.1) beta=1
            if (kk.eq.2.and.ll.eq.2) beta=2
            if (kk.eq.3.and.ll.eq.3) beta=3
            if (kk.eq.2.and.ll.eq.3) beta=4
            if (kk.eq.1.and.ll.eq.3) beta=5
            if (kk.eq.1.and.ll.eq.2) beta=6
            if (kk.eq.3.and.ll.eq.2) beta=4
            if (kk.eq.3.and.ll.eq.1) beta=5
            if (kk.eq.2.and.ll.eq.1) beta=6
            alpha_v_tensor(ii,jj)=alpha_v_tensor(ii,jj)+Lattice%Sij(alpha,beta)*Gama_tensor(kk,ll)*C_v*kb_HaK*Ha_J/1.d9/(Lattice%ucvol*Bohr_Ang**3*1.d-30)
          end do
        end do
        P_th_tensor(ii,jj)   =Gama_tensor(ii,jj)*E_th/(Lattice%ucvol*Bohr_Ang**3*1.d-30)*e_Cb/10**9
      end do
    end do
    P_th2  =      p_thermo2(itemp)*k_B/(Lattice%ucvol*Bohr_Ang**3*1.d-30)*e_Cb/10**9
    if (MPIdata%iam_master) then
      write(20,'(1x,i5,7(1x,f15.5))') itemp*10,C_v,(Gama_tensor(1,1)+Gama_tensor(2,2)+Gama_tensor(3,3))/3.d0,&
&                                                  (alpha_v_tensor(1,1)+alpha_v_tensor(2,2)+alpha_v_tensor(3,3))*1.d6,&
&                                             E_th,(P_th1_tensor(1,1)+P_th1_tensor(2,2)+P_th1_tensor(3,3))/3.d0,&
&                                            P_th2,(P_th_tensor(1,1)+P_th_tensor(2,2)+P_th_tensor(3,3))/3.d0
      write(21,'(1x,i5,9(1x,f15.5))') itemp*10,Gama_tensor(1,1),Gama_tensor(2,2),Gama_tensor(3,3),&
&                                     alpha_v_tensor(1,1)*1.d6,alpha_v_tensor(2,2)*1.d6,alpha_v_tensor(3,3)*1.d6,&
&                                     alpha_v_tensor(1,2)*1.d6,alpha_v_tensor(1,3)*1.d6,alpha_v_tensor(2,3)*1.d6
    end if
  end do  
  if (MPIdata%iam_master) then
    close(20)
    close(21)
  end if  
  ABI_FREE(heatcapa_HA)
  ABI_FREE(grun_thermo_HA)
  ABI_FREE(p_thermo_HA)
  ABI_FREE(u_vib_HA)

! Summary
! =======
  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '####### Gruneisen parameter, Thermal expansion, Thermal pressure... #########'
  write(InVar%stdout,*) '#############################################################################'
  C_v    =sum(heatcapa(:))
  E_th   =sum(u_vib   (:))
  alpha_v_tensor(:,:)=zero
  Gama_tensor   (:,:)=zero
  P_th1_tensor  (:,:)=zero
  P_th_tensor   (:,:)=zero
  do imode=1,nmode
    Gama_tensor (:,:)=Gama_tensor (:,:)+grun_thermo(imode,:,:)/C_v
    P_th1_tensor(:,:)=P_th1_tensor(:,:)+p_thermo1  (imode,:,:)/(Lattice%ucvol*Bohr_Ang**3*1.d-30)*e_Cb/10**9
  end do
! alpha=sum{I=1,6;J=1,6} S(I,J) Gama_tensor (J)     
  do ii=1,3
    do jj=1,3
      if (ii.eq.1.and.jj.eq.1) alpha=1
      if (ii.eq.2.and.jj.eq.2) alpha=2
      if (ii.eq.3.and.jj.eq.3) alpha=3
      if (ii.eq.2.and.jj.eq.3) alpha=4
      if (ii.eq.1.and.jj.eq.3) alpha=5
      if (ii.eq.1.and.jj.eq.2) alpha=6
      if (ii.eq.3.and.jj.eq.2) alpha=4
      if (ii.eq.3.and.jj.eq.1) alpha=5
      if (ii.eq.2.and.jj.eq.1) alpha=6
      do kk=1,3
        do ll=1,3
          if (kk.eq.1.and.ll.eq.1) beta=1
          if (kk.eq.2.and.ll.eq.2) beta=2
          if (kk.eq.3.and.ll.eq.3) beta=3
          if (kk.eq.2.and.ll.eq.3) beta=4
          if (kk.eq.1.and.ll.eq.3) beta=5
          if (kk.eq.1.and.ll.eq.2) beta=6
          if (kk.eq.3.and.ll.eq.2) beta=4
          if (kk.eq.3.and.ll.eq.1) beta=5
          if (kk.eq.2.and.ll.eq.1) beta=6
          alpha_v_tensor(ii,jj)=alpha_v_tensor(ii,jj)+Lattice%Sij(alpha,beta)*Gama_tensor(kk,ll)*C_v*kb_HaK*Ha_J/1.d9/(Lattice%ucvol*Bohr_Ang**3*1.d-30)
        end do
      end do
      P_th_tensor(ii,jj)=Gama_tensor(ii,jj)*E_th/(Lattice%ucvol*Bohr_Ang**3*1.d-30)*e_Cb/10**9
    end do  
  end do  
  P_th2  =p_thermo2(int(InVar%temperature/10))*k_B/(Lattice%ucvol*Bohr_Ang**3*1.d-30)*e_Cb/10**9
  Lattice%BulkModulus_S=Lattice%BulkModulus_T*(1.+(alpha_v_tensor(1,1)+alpha_v_tensor(2,2)+alpha_v_tensor(3,3))*&
&                                                    (Gama_tensor(1,1)+   Gama_tensor(2,2)+   Gama_tensor(3,3))/3.d0*InVar%temperature)
  Lattice%HeatCapa_P=C_v*Lattice%BulkModulus_S/Lattice%BulkModulus_T
  Vp=dsqrt(1.d9*(Lattice%BulkModulus_S+4.d0*Lattice%Shear/3.d0)/Lattice%Density)
  Vs=dsqrt(1.d9*Lattice%Shear/Lattice%Density)
  write(InVar%stdout,'(a)') ' See the gruneisen.dat, alpha_gamma.dat and thermo3.dat files'
  write(InVar%stdout,'(a)') ' '
  write(InVar%stdout,'(a,1x,f15.5)') ' Gruneisen parameter :                 Gamma=',(Gama_tensor(1,1)+Gama_tensor(2,2)+Gama_tensor(3,3))/3.d0
  write(InVar%stdout,'(a,3(1x,f8.3),a)') '                The Gruneisen matrix is  |',Gama_tensor(1,1),Gama_tensor(1,2),Gama_tensor(1,3),'|'
  write(InVar%stdout,'(a,3(1x,f8.3),a)') '                                         |',Gama_tensor(2,1),Gama_tensor(2,2),Gama_tensor(2,3),'|'
  write(InVar%stdout,'(a,3(1x,f8.3),a)') '                                         |',Gama_tensor(3,1),Gama_tensor(3,2),Gama_tensor(3,3),'|'
  write(InVar%stdout,'(a,1x,f15.5)') ' Thermal expansion (K^{-1}*10^6) :   alpha_v=',(alpha_v_tensor(1,1)+alpha_v_tensor(2,2)+alpha_v_tensor(3,3))*1.d6
  write(InVar%stdout,'(a,3(1x,f8.3),a)') '        The thermal expansion matrix is  |',alpha_v_tensor(1,1)*1.d6,alpha_v_tensor(1,2)*1.d6,alpha_v_tensor(1,3)*1.d6,'|'
  write(InVar%stdout,'(a,3(1x,f8.3),a)') '                                         |',alpha_v_tensor(2,1)*1.d6,alpha_v_tensor(2,2)*1.d6,alpha_v_tensor(2,3)*1.d6,'|'
  write(InVar%stdout,'(a,3(1x,f8.3),a)') '                                         |',alpha_v_tensor(3,1)*1.d6,alpha_v_tensor(3,2)*1.d6,alpha_v_tensor(3,3)*1.d6,'|'
  write(InVar%stdout,'(a)') ' Thermal pressure (in GPa) : '
  write(InVar%stdout,'(a)') '    ------- w   intrinsic effects and w   ZPE --------'
  write(InVar%stdout,'(a,1x,f15.5)') '                 P_th=sum_i Gamma_i*E_i/V   =',(P_th1_tensor(1,1)+P_th1_tensor(2,2)+P_th1_tensor(3,3))/3.d0
  write(InVar%stdout,'(a,3(1x,f8.3),a)') '           The thermal pressure matrix is|',P_th1_tensor(1,1),P_th1_tensor(1,2),P_th1_tensor(1,3),'|'
  write(InVar%stdout,'(a,3(1x,f8.3),a)') '                                         |',P_th1_tensor(2,1),P_th1_tensor(2,2),P_th1_tensor(2,3),'|'
  write(InVar%stdout,'(a,3(1x,f8.3),a)') '                                         |',P_th1_tensor(3,1),P_th1_tensor(3,2),P_th1_tensor(3,3),'|'
  write(InVar%stdout,'(a)') '    ------- w   intrinsic effects and w/o ZPE --------'
  write(InVar%stdout,'(a,1x,f15.5)') '                 P_th=integ{Gamma*C_v/V dT} =',P_th2
  write(InVar%stdout,'(a)') '    ------- w/o intrinsic effects and w/o ZPE --------'
  write(InVar%stdout,'(a,1x,f15.5)') '                 P_th=Gamma*E_th/V          =',(P_th_tensor(1,1)+P_th_tensor(2,2)+P_th_tensor(3,3))/3.d0
  write(InVar%stdout,'(a,3(1x,f8.3),a)') '           The thermal pressure matrix is|',P_th_tensor(1,1),P_th_tensor(1,2),P_th_tensor(1,3),'|'
  write(InVar%stdout,'(a,3(1x,f8.3),a)') '                                         |',P_th_tensor(2,1),P_th_tensor(2,2),P_th_tensor(2,3),'|'
  write(InVar%stdout,'(a,3(1x,f8.3),a)') '                                         |',P_th_tensor(3,1),P_th_tensor(3,2),P_th_tensor(3,3),'|'
  write(InVar%stdout,'(a,1x,f15.5)') ' Volume (bohr^3 per unit cell):            V=',Lattice%ucvol
  write(InVar%stdout,'(a,1x,f15.5)') ' Thermal energy (eV) :                  E_th=',E_th
  write(InVar%stdout,'(a,1x,f15.5)') ' Heat capacity at constant V (k_B/f.u.): C_v=',C_v
  write(InVar%stdout,'(a,1x,f15.5)') ' Heat capacity at constant P (k_B/f.u.): C_p=',Lattice%HeatCapa_P
  write(InVar%stdout,'(a,1x,f15.5)') ' Isothermal Bulk Modulus (GPa):          B_T=',Lattice%BulkModulus_T
  write(InVar%stdout,'(a,1x,f15.5)') ' Isentropic Bulk Modulus (GPa):          B_S=',Lattice%BulkModulus_S
  write(InVar%stdout,'(a,1x,f15.5)') ' Longitudinal sound velocity (m.s-1):     Vp=',Vp
  write(InVar%stdout,'(a,1x,f15.5)') ' Transverse sound velocity (m.s-1):       Vs=',Vs
  ABI_FREE(heatcapa)
  ABI_FREE(grun_thermo)
  ABI_FREE(p_thermo1)
  ABI_FREE(p_thermo2)
  ABI_FREE(u_vib)

!FB  0.473294364993209*(2.38255605878933*1.38065e-23)/1.20512e11/3.0/((30.6135754000*0.529177e-10)**3/216)

end subroutine tdep_calc_alpha_gamma

!=====================================================================================================
subroutine tdep_write_gruneisen(Crystal,distance,Eigen2nd,Ifc,InVar,Lattice,Phi3_ref,Qpt,Rlatt_cart,Shell3at,Sym)

  implicit none

  type(crystal_t),intent(in) :: Crystal
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(in) :: Shell3at
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Eigen_Variables_type),intent(in) :: Eigen2nd
  type(Qpoints_type),intent(in) :: Qpt
  type(ifc_type),intent(in) :: Ifc
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision,intent(in) :: Phi3_ref(3,3,3,Shell3at%nshell)
  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)

  integer :: iqpt,imode,nmode,ii,jj,kk,iatcell,jatcell
  double precision :: qpt_cart(3)
  double complex, allocatable :: Gruneisen(:,:,:)
  double complex, allocatable :: Grun_mean(:)
  character(len=500) :: message
  integer :: ierr

  ierr = 0

  nmode=3*InVar%natom_unitcell
  ABI_MALLOC(Gruneisen,(3*InVar%natom_unitcell,3,3)); Gruneisen(:,:,:)=czero
  ABI_MALLOC(Grun_mean,(3*InVar%natom_unitcell))    ; Grun_mean(:)    =czero
  open(unit=53,file=trim(InVar%output_prefix)//'gruneisen.dat')
  open(unit=54,file=trim(InVar%output_prefix)//'gruneisen-ij.dat')

  do iqpt=1,Qpt%nqpt

!   Compute the Gruneisen
!   =====================
    Grun_mean(:)=czero
    qpt_cart(:)=Qpt%qpt_cart(:,iqpt)
    if ((sum(abs(Qpt%qpt_red(:,iqpt)))).lt.tol8) cycle  ! G point
    if (abs(sum(Qpt%qpt_red(:,iqpt)**2)-1.d0).lt.tol8) cycle ! Gp point
    call tdep_calc_gruneisen(distance,Eigen2nd,Gruneisen,iqpt,InVar,Lattice,Phi3_ref,qpt_cart,Rlatt_cart,Shell3at,Sym)
    do ii=1,3*InVar%natom_unitcell
      do jj=1,3
        do kk=1,3
          write(54,'(4(i5,1x),500(e15.6,1x))') iqpt,ii,jj,kk,Gruneisen(ii,jj,kk)
          if (abs(aimag(Gruneisen(ii,jj,kk))).gt.tol8) then
            write(message,'(4(i5,1x),100(e15.6,1x))') iqpt,ii,jj,kk,real(Gruneisen(ii,jj,kk))
            MSG_ERROR_NOSTOP(message,ierr)
            write(message,'(4(i5,1x),100(e15.6,1x))') iqpt,ii,jj,kk,aimag(Gruneisen(ii,jj,kk))
            MSG_ERROR_NOSTOP(message,ierr)
            MSG_ERROR('tdep_write_gruneisen : The imaginary part of the Gruneisen is not equal to zero')
         end if  
        end do  
      end do  
    end do
    do ii=1,3*InVar%natom_unitcell
      do jj=1,3
        do kk=1,3
          Grun_mean(ii)=Grun_mean(ii)+Gruneisen(ii,jj,kk)
        end do  
      end do  
    end do

!   Write the Gruneisen
!   ===================
    if (sum(abs(aimag(Grun_mean(:)))).gt.3*InVar%natom_unitcell*tol8) then
      write(message,'(i5,1x,100(e15.6,1x))') iqpt,(real(Grun_mean(ii)),ii=1,nmode)
      MSG_ERROR_NOSTOP(message,ierr)
      write(message,'(i5,1x,100(e15.6,1x))') iqpt,(aimag(Grun_mean(ii)),ii=1,nmode)
      MSG_ERROR_NOSTOP(message,ierr)
      MSG_ERROR('tdep_write_gruneisen : The imaginary part of Grun_mean is not equal to zero')
    else 
!FB      write(53,'(i5,1x,500(e15.6,1x))') iqpt,(real(Grun_mean(ii)),ii=1,nmode),((real(Grun_shell(ii,jj)),ii=1,nmode),jj=1,Shell3at%nshell)
      write(53,'(i5,1x,500(e15.6,1x))') iqpt,(real(Grun_mean(ii)),ii=1,nmode)
    end if  
  end do  
  close(53)
  close(54)
  ABI_FREE(Grun_mean)
  ABI_FREE(Gruneisen)

end subroutine tdep_write_gruneisen

!=====================================================================================================
subroutine tdep_calc_lifetime1(Crystal,distance,Eigen2nd,Ifc,InVar,Lattice,Phi3_ref,Qbz,Rlatt_cart,Shell3at,Sym)

  implicit none

  type(crystal_t),intent(in) :: Crystal
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(in) :: Shell3at
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Eigen_Variables_type),intent(in) :: Eigen2nd
  type(Qbz_type),intent(in) :: Qbz
  type(ifc_type),intent(in) :: Ifc
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision,intent(in) :: Phi3_ref(3,3,3,Shell3at%nshell)
  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)

  integer :: iq_bz,jq_bz,kq_bz,imode,nmode,ii,jj,kk,iatcell,jatcell,tot_count
  integer :: okp_count,okm_count
  double precision :: iqbz(3),jqbz(3),kqbz(3),ipjpk(3),imjpk(3)
!  double complex, allocatable :: Gruneisen(:,:,:)
!  double complex, allocatable :: Grun_mean(:)
  character(len=500) :: message
  integer :: ierr

  ierr = 0

  nmode=3*InVar%natom_unitcell
  tot_count=zero
  okp_count=zero
  okm_count=zero
  do ii=1,3
    write(InVar%stdout,'(a,3(e15.6,1x),1x)') 'Gmet=',Lattice%gmet(ii,:)
    write(InVar%stdout,'(a,3(e15.6,1x),1x)') 'Gprim=',Lattice%gprim(ii,:)
  end do
  do iq_bz=1,Qbz%nqbz
    iqbz(:)=Qbz%qbz(:,iq_bz)
!    iqbz(:)=0.2
    write(InVar%stdout,'(a,3(e15.6,1x),1x)') 'Qbz1=',iqbz(:)
    do jq_bz=1,Qbz%nqbz
      jqbz(:)=Qbz%qbz(:,jq_bz)
!      jqbz(:)=0.3
      do kq_bz=1,Qbz%nqbz
        tot_count=tot_count+1
        kqbz(:)=Qbz%qbz(:,kq_bz)
!        kqbz(:)=0.5
        ipjpk(:)=iqbz(:)+jqbz(:)+kqbz(:)
        imjpk(:)=iqbz(:)-jqbz(:)+kqbz(:)
        if (sum(abs(ipjpk(:)-int(ipjpk(:)))).lt.tol4) then
          write(InVar%stdout,'(a,i4,1x,i4,1x,i4)')'OK for qi+qj+qk=G',iq_bz,jq_bz,kq_bz  
          okp_count=okp_count+1
	else if (sum(abs(imjpk(:)-int(imjpk(:)))).lt.tol4) then
          write(InVar%stdout,'(a,i4,1x,i4,1x,i4)')'OK for qi-qj+qk=G',iq_bz,jq_bz,kq_bz  
          okm_count=okm_count+1
        else 
          cycle
        end if  
!       do imode=1,nmode
!         do jmode=1,nmode
!           do kmode=1,nmode
!            end do  
!          end do  
!        end do  
      end do  
    end do  
  end do  
  write(InVar%stdout,'(a,i10)') 'Total number of (q1,q2,q3) =',tot_count
  write(InVar%stdout,'(a,i10)') 'Partial number of (q1+q2+q3) =',okp_count
  write(InVar%stdout,'(a,i10)') 'Partial number of (q1-q2+q3) =',okm_count

end subroutine tdep_calc_lifetime1

!=====================================================================================================
end module m_tdep_phi3

!=====================================================================================================
!FBsubroutine tdep_calc_lifetime2(Crystal,distance,Eigen2nd,Ifc,InVar,Lattice,MPIdata,Phi3_ref,Qbz,Rlatt_cart,Shell3at,Sym)
!FB
!FB!Arguments ------------------------------------
!FB!scalars
!FB integer,intent(in) :: ncid,comm
!FB real(dp),intent(in) :: dosdeltae !,dossmear
!FB type(gruns_t),intent(in) :: gruns
!FB character(len=*),intent(in) :: prefix
!FB!arrays
!FB
!FB  type(crystal_t),intent(in) :: Crystal
!FB  type(Symetries_Variables_type),intent(in) :: Sym
!FB  type(Input_Variables_type),intent(in) :: InVar
!FB  type(Shell_Variables_type),intent(in) :: Shell3at
!FB  type(Lattice_Variables_type),intent(in) :: Lattice
!FB  type(Eigen_Variables_type),intent(in) :: Eigen2nd
!FB  type(MPI_enreg_type), intent(in) :: MPIdata
!FB  type(Qbz_type),intent(in) :: Qbz
!FB  type(ifc_type),intent(in) :: Ifc
!FB  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4)
!FB  double precision,intent(in) :: Phi3_ref(3,3,3,Shell3at%nshell)
!FB  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)
!FB
!FB!Local variables-------------------------------
!FB!scalars
!FB integer,parameter :: master=0,qptopt1=1,bcorr0=0
!FB integer :: nprocs,my_rank,iqibz,nqbz,nqibz,ierr,ii,nu,ncerr,nomega,cnt,unt,io
!FB real(dp) :: gavg,omega_min,omega_max,v2
!FB type(htetra_t) :: tetra
!FB character(len=500) :: msg
!FB!arrays
!FB integer :: qptrlatt(3,3)
!FB integer :: ngqpt(3)
!FB real(dp) :: shiftq(3,1)
!FB real(dp),allocatable :: gvals_qibz(:,:),wvols_qibz(:,:,:),dwdq_qibz(:,:,:)
!FB real(dp),allocatable :: qibz(:,:),qbz(:,:),wtq(:)
!FB real(dp),allocatable :: wdt(:,:),wdos(:,:),grdos(:,:),gr2dos(:,:),wibz(:),omega_mesh(:)
!FB real(dp),allocatable :: vdos(:,:),v2dos(:,:)
!FB real(dp),allocatable :: phdispl_cart_qibz(:,:,:,:,:)
!FB
!FB! ************************************************************************
!FB
!FB comm = MPIdata%comm_step
!FB nprocs  = xmpi_comm_size(comm) 
!FB my_rank = xmpi_comm_rank(comm)
!FB
!FB write(msg,'(a,(80a),4a)')ch10,('=',ii=1,80),ch10,ch10,' Calculation of Lifetime ',ch10
!FB call wrtout(InVar%stdout,msg)
!FB
!FB ! Generate the q-mesh by finding the IBZ and the corresponding weights.
!FB qptrlatt = 0
!FB do ii=1,3
!FB   qptrlatt(ii,ii) = InVar%ngqpt1(:)
!FB end do
!FB nshiftq=1
!FB shiftq(:,:)=0.0
!FB
!FB ! Get IBZ and BZ.
!FB call kpts_ibz_from_kptrlatt(Crystal, qptrlatt, qptopt1, nshiftq, shiftq, &
!FB   nqibz, qibz, wtq, nqbz, qbz)
!FB
!FB ! Build tetrahedra
!FB tetra = tetra_from_kptrlatt(Crystal, qptopt1, qptrlatt, nshiftq, shiftq, nqibz, qibz, comm, msg, ierr)
!FB if (ierr /= 0) MSG_ERROR(msg)
!FB
!FB ABI_CALLOC(wvols_qibz, (gruns%natom3, gruns%nvols, nqibz))
!FB ABI_CALLOC(gvals_qibz, (gruns%natom3, nqibz))
!FB ABI_CALLOC(dwdq_qibz, (3, gruns%natom3, nqibz))
!FB ABI_CALLOC(phdispl_cart_qibz, (2, gruns%natom3, gruns%natom3, gruns%nvols, nqibz))
!FB
!FB gavg = zero
!FB do iqibz=1,nqibz
!FB   if (mod(iqibz, nprocs) /= my_rank) cycle ! mpi-parallelism
!FB   call gruns_fourq(gruns, qibz(:,iqibz), wvols_qibz(:,:,iqibz), gvals_qibz(:,iqibz), &
!FB                    dwdq_qibz(:,:,iqibz), phdispl_cart_qibz(:,:,:,:,iqibz))
!FB   gavg = gavg + wtq(iqibz) * sum(gvals_qibz(:,iqibz))
!FB end do
!FB gavg = gavg / gruns%natom3
!FB
!FB call xmpi_sum(gavg, comm, ierr)
!FB call xmpi_sum(wvols_qibz, comm, ierr)
!FB call xmpi_sum(gvals_qibz, comm, ierr)
!FB call xmpi_sum(dwdq_qibz, comm, ierr)
!FB call xmpi_sum(phdispl_cart_qibz, comm, ierr)
!FB
!FB omega_min = gruns%ifc_vol(gruns%iv0)%omega_minmax(1)
!FB omega_max = gruns%ifc_vol(gruns%iv0)%omega_minmax(2)
!FB nomega = nint((omega_max - omega_min) / dosdeltae) + 1
!FB nomega = max(6, nomega) ! Ensure Simpson integration will be ok
!FB
!FB ABI_MALLOC(omega_mesh, (nomega))
!FB omega_mesh = arth(omega_min, dosdeltae, nomega)
!FB omega_max = omega_mesh(nomega)
!FB !write(std_out,*)"hello",omega_min,omega_max,dosdeltae,(omega_max-omega_min) / (nomega-1)
!FB ABI_MALLOC(wibz, (nqibz))
!FB ABI_MALLOC(wdt, (nomega, 2))
!FB ABI_CALLOC(wdos, (nomega, 2))
!FB ABI_CALLOC(grdos, (nomega, 2))
!FB ABI_CALLOC(gr2dos, (nomega, 2))
!FB ABI_CALLOC(vdos, (nomega, 2))
!FB ABI_CALLOC(v2dos, (nomega, 2))
!FB
!FB ! Compute DOSes.
!FB cnt = 0
!FB do iqibz=1,nqibz
!FB   do nu=1,gruns%natom3
!FB     cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! mpi-parallelism
!FB     wibz = wvols_qibz(nu, gruns%iv0, :)
!FB     call tetra%get_onewk(iqibz,bcorr0,nomega,nqibz,wibz,omega_min,omega_max,one,wdt)
!FB     wdt = wdt*wtq(iqibz)
!FB     wdos = wdos + wdt
!FB     grdos = grdos + wdt * gvals_qibz(nu,iqibz)
!FB     gr2dos = gr2dos + wdt * gvals_qibz(nu,iqibz) ** 2
!FB     v2 = sum(dwdq_qibz(1:3,nu,iqibz) ** 2)
!FB     vdos = vdos + wdt * sqrt(v2)
!FB     v2dos = v2dos + wdt * v2
!FB   end do
!FB end do
!FB
!FB call xmpi_sum(wdos, comm, ierr)
!FB call xmpi_sum(grdos, comm, ierr)
!FB call xmpi_sum(gr2dos, comm, ierr)
!FB call xmpi_sum(vdos, comm, ierr)
!FB call xmpi_sum(v2dos, comm, ierr)
!FB
!FB if (my_rank == master) then
!FB   write(unt,'(a)')'# Phonon density of states, Gruneisen DOS and phonon group velocity DOS'
!FB   write(unt,'(a)')"# Energy in Hartree, DOS in states/Hartree"
!FB   write(unt,'(a,i0)')'# Tetrahedron method with nqibz= ',nqibz
!FB   write(unt,"(a,f8.5)")"# Average Gruneisen parameter:", gavg
!FB   write(unt,'(5a)') &
!FB     "# omega PH_DOS Gruns_DOS Gruns**2_DOS Vel_DOS  Vel**2_DOS  PH_IDOS Gruns_IDOS Gruns**2_IDOS Vel_IDOS Vel**2_IDOS"
!FB   do io=1,nomega
!FB     write(unt, "(11es17.8)")omega_mesh(io), &
!FB       wdos(io,1), grdos(io,1), gr2dos(io,1), vdos(io,1), v2dos(io,1), &
!FB       wdos(io,2), grdos(io,2), gr2dos(io,2), vdos(io,2), v2dos(io,2)
!FB   end do
!FB end if
!FB
!FB
!FB ABI_FREE(qibz)
!FB ABI_FREE(wtq)
!FB ABI_FREE(qbz)
!FB ABI_FREE(wvols_qibz)
!FB ABI_FREE(gvals_qibz)
!FB ABI_FREE(dwdq_qibz)
!FB ABI_FREE(phdispl_cart_qibz)
!FB ABI_FREE(omega_mesh)
!FB ABI_FREE(wibz)
!FB ABI_FREE(wdt)
!FB ABI_FREE(wdos)
!FB ABI_FREE(grdos)
!FB ABI_FREE(gr2dos)
!FB ABI_FREE(v2dos)
!FB ABI_FREE(vdos)
!FB
!FB call tetra%free()
!FB
!FBend subroutine tdep_calc_lifetime2
!!***
