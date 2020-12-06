
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_phi2

  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_io_tools
  use m_tdep_readwrite,   only : Input_Variables_type, MPI_enreg_type
  use m_tdep_latt,        only : Lattice_Variables_type
  use m_tdep_shell,       only : Shell_Variables_type
  use m_tdep_sym,         only : Symetries_Variables_type
  use m_tdep_qpt,         only : Qpoints_type
  use m_tdep_utils,       only : Coeff_Moore_type

  implicit none

  type Eigen_Variables_type

    double precision, allocatable :: eigenval(:,:)
    double precision, allocatable :: eigenvec(:,:,:,:,:,:)
    double precision, allocatable :: dynmat(:,:,:,:,:,:)

  end type Eigen_Variables_type

  public :: tdep_calc_ftot2
  public :: tdep_calc_phi1fcoeff
  public :: tdep_calc_phi1
  public :: tdep_write_phi1
  public :: tdep_calc_phi2fcoeff
  public :: tdep_calc_phi2
  public :: tdep_write_phi2
  public :: tdep_build_phi2_33
  public :: tdep_calc_dij
  public :: tdep_write_dij
  public :: tdep_init_eigen2nd
  public :: tdep_destroy_eigen2nd
  public :: tdep_write_yaml

contains

!====================================================================================================
 subroutine tdep_calc_ftot2(Forces_TDEP,InVar,Phi1,Phi1Ui,Phi2,Phi2UiUj,ucart) 

  implicit none 

  type(Input_Variables_type),intent(in) :: InVar
  double precision, intent(in)  :: Phi2(3*InVar%natom,3*InVar%natom)
  double precision, intent(in)  :: ucart(3,InVar%natom,InVar%my_nstep)
  double precision, intent(in)  :: Phi1(3*InVar%natom)
  double precision, intent(out) :: Phi1Ui(InVar%my_nstep)
  double precision, intent(out) :: Phi2UiUj(InVar%my_nstep)
  double precision, intent(inout) :: Forces_TDEP(3*InVar%natom*InVar%my_nstep)
  
  integer :: jj,istep,jatom,katom
  double precision, allocatable :: ucart_blas(:)
  double precision, allocatable :: ftot2(:)

! Compute Forces of the model TDEP
  ABI_MALLOC(ucart_blas,(3*InVar%natom)); ucart_blas(:)=0.d0
  ABI_MALLOC(ftot2     ,(3*InVar%natom)); ftot2     (:)=0.d0
  do istep=1,InVar%my_nstep
    ucart_blas(:)=0.d0 
    ftot2     (:)=0.d0 
    do jatom=1,InVar%natom
      do jj=1,3
        ucart_blas(3*(jatom-1)+jj)=ucart(jj,jatom,istep)
      end do
    end do
    Phi1Ui(istep)=sum(Phi1(:)*ucart_blas(:))
    call DGEMM('N','N',3*InVar%natom,1,3*InVar%natom,1.d0,Phi2,3*InVar%natom,ucart_blas,3*InVar%natom,0.d0,ftot2(:),3*InVar%natom)
    call DGEMM('T','N',1,1,3*InVar%natom,1./2.d0,ftot2,3*InVar%natom,ucart_blas,3*InVar%natom,0.d0,Phi2UiUj(istep),3*InVar%natom)
    Forces_TDEP(3*InVar%natom*(istep-1)+1:3*InVar%natom*istep)=-Phi1(:)-ftot2(:)
  end do !istep  
  ABI_FREE(ucart_blas)
  ABI_FREE(ftot2)

 end subroutine tdep_calc_ftot2

!====================================================================================================
subroutine tdep_calc_phi1fcoeff(CoeffMoore,InVar,proj,Shell1at,Sym)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell1at
  type(Coeff_Moore_type), intent(inout) :: CoeffMoore
  double precision, intent(in) :: proj(3,3,Shell1at%nshell)

  integer :: ishell,ncoeff,ncoeff_prev,istep,iatom,iatshell,iat_mod
  integer :: icoeff,isym,mu,iatref
  double precision :: terme

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '############## Fill the matrices used in the pseudo-inverse #################'
  write(InVar%stdout,*) '#############################################################################'

  write(InVar%stdout,*) ' Compute the coefficients (at the 1st order) used in the Moore-Penrose...'
  do ishell=1,Shell1at%nshell
    if (Shell1at%neighbours(1,ishell)%n_interactions.eq.0) cycle
    do iatshell=1,Shell1at%neighbours(1,ishell)%n_interactions
      iatom=Shell1at%neighbours(1,ishell)%atomj_in_shell(iatshell) 
      iat_mod=mod(iatom+InVar%natom_unitcell-1,InVar%natom_unitcell)+1
      if (iat_mod==1) cycle
      iatref=Shell1at%iatref(ishell)
      isym=Shell1at%neighbours(1,ishell)%sym_in_shell(iatshell)
      ncoeff     =Shell1at%ncoeff(ishell)
      ncoeff_prev=Shell1at%ncoeff_prev(ishell)
      do mu=1,3
        do icoeff=1,ncoeff
          terme=sum(Sym%S_ref(mu,:,isym,1)*proj(:,icoeff,ishell))
          do istep=1,InVar%my_nstep
            CoeffMoore%fcoeff(mu+3*(iatom-1)+3*InVar%natom*(istep-1),icoeff+ncoeff_prev)= &
&           CoeffMoore%fcoeff(mu+3*(iatom-1)+3*InVar%natom*(istep-1),icoeff+ncoeff_prev)+terme
!           Add all the other contributions, when iat_mod==1 (due to ASR)
            CoeffMoore%fcoeff(mu+3*(iatom-iat_mod+1)+3*InVar%natom*(istep-1),icoeff+ncoeff_prev)= &
&           CoeffMoore%fcoeff(mu+3*(iatom-iat_mod+1)+3*InVar%natom*(istep-1),icoeff+ncoeff_prev)-terme
          end do !istep
        end do    
      end do  
    end do !iatshell
  end do !ishell
  write(InVar%stdout,*) ' ------- achieved'

end subroutine tdep_calc_phi1fcoeff

!====================================================================================================
subroutine tdep_calc_phi2fcoeff(CoeffMoore,InVar,proj,Shell2at,Sym,ucart)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell2at
  type(Coeff_Moore_type), intent(inout) :: CoeffMoore
  double precision, intent(in) :: ucart(3,InVar%natom,InVar%my_nstep)
  double precision, intent(in) :: proj(9,9,Shell2at%nshell)

  integer :: ishell,ncoeff,ncoeff_prev,istep,iatom,jatom,iatshell
  integer :: icoeff,isym
  integer :: mu,nu,alpha,beta,itrans
  double precision :: terme,temp
  double precision :: udiff(3),SSu(3,9)
  double precision, allocatable :: SS_ref(:,:,:,:,:)

! For each couple of atoms, transform the Phi2 (3x3) ifc matrix using the symetry operation (S)
! Note: iatom=1 is excluded in order to take into account the atomic sum rule (see below)
  ABI_MALLOC(SS_ref,(3,9,3,Sym%nsym,2)); SS_ref(:,:,:,:,:)=zero
  do isym=1,Sym%nsym
    do mu=1,3
      do alpha=1,3
        do nu=1,3
          do beta=1,3
            temp=Sym%S_ref(mu,alpha,isym,1)*Sym%S_ref(nu,beta,isym,1)
            SS_ref(mu,beta+(alpha-1)*3,nu,isym,1)=temp
            SS_ref(mu,alpha+(beta-1)*3,nu,isym,2)=temp
          end do
        end do
      end do
    end do
  end do  
        
  write(InVar%stdout,*) ' Compute the coefficients (at the 2nd order) used in the Moore-Penrose...'
  do ishell=1,Shell2at%nshell
    do iatom=1,InVar%natom
      if (Shell2at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
      do iatshell=1,Shell2at%neighbours(iatom,ishell)%n_interactions
        jatom=Shell2at%neighbours(iatom,ishell)%atomj_in_shell(iatshell) 
        if (iatom==jatom) cycle
        isym=Shell2at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
        itrans=Shell2at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
        ncoeff     =Shell2at%ncoeff(ishell)
        ncoeff_prev=Shell2at%ncoeff_prev(ishell)+CoeffMoore%ncoeff1st

        do istep=1,InVar%my_nstep
!         In order to impose the acoustic sum rule we use (u(j)-u(i))==u_j^\nu
          udiff(1)=ucart(1,jatom,istep)-ucart(1,iatom,istep)
          udiff(2)=ucart(2,jatom,istep)-ucart(2,iatom,istep)
          udiff(3)=ucart(3,jatom,istep)-ucart(3,iatom,istep)
            
!         F_i^\mu(t)=\sum_{\alpha\beta,j,\nu}S^{\mu\alpha}.S^{\nu\beta}.\Phi_{ij}^{\alpha\beta}.u_j^\nu(t)
          SSu(:,:)=zero
          do nu=1,3
            SSu(:,:)=SSu(:,:)+SS_ref(:,:,nu,isym,itrans)*udiff(nu)
          end do  
          do mu=1,3
            do icoeff=1,ncoeff
              terme=sum(SSu(mu,:)*proj(:,icoeff,ishell))
!FB              write(6,*) 'indices=', mu+3*(iatom-1)+3*InVar%natom*(istep-1),icoeff+ncoeff_prev
              CoeffMoore%fcoeff(mu+3*(iatom-1)+3*InVar%natom*(istep-1),icoeff+ncoeff_prev)= &
&             CoeffMoore%fcoeff(mu+3*(iatom-1)+3*InVar%natom*(istep-1),icoeff+ncoeff_prev)+terme
            end do    
          end do  
        
        end do !istep
      end do !iatshell
    end do !iatom
  end do !ishell
  write(InVar%stdout,*) ' ------- achieved'
  ABI_FREE(SS_ref)

end subroutine tdep_calc_phi2fcoeff

!=====================================================================================================
subroutine tdep_calc_phi1(InVar,ntotcoeff,proj,Phi1_coeff,Phi1,Shell1at,Sym)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell1at
  integer,intent(in) :: ntotcoeff
  double precision,intent(in) :: proj(3,3,Shell1at%nshell)
  double precision,intent(in) :: Phi1_coeff(ntotcoeff,1)
  double precision,intent(out) :: Phi1(3*InVar%natom)

  integer :: ishell,isym,iatom,ncoeff,ncoeff_prev
  integer :: nshell,ii,iatshell,iat_mod
  double precision,allocatable :: Phi1_3(:),Phi1_ref(:,:)

  nshell=Shell1at%nshell
  ABI_MALLOC(Phi1_ref,(3,nshell)); Phi1_ref(:,:)=zero
  ABI_MALLOC(Phi1_3,(3)) ; Phi1_3(:)=0.d0
  do ishell=1,nshell
!   Build the 3x3 IFC per shell    
    ncoeff     =Shell1at%ncoeff(ishell)
    ncoeff_prev=Shell1at%ncoeff_prev(ishell)
    do ii=1,3
      Phi1_ref(ii,ishell)=sum(proj(ii,1:ncoeff,ishell)*Phi1_coeff(ncoeff_prev+1:ncoeff_prev+ncoeff,1))
    end do  
!   Build the vector-IFC of an atom in this shell    
    if (Shell1at%neighbours(1,ishell)%n_interactions.eq.0) cycle
    do iatshell=1,Shell1at%neighbours(1,ishell)%n_interactions
      iatom=Shell1at%neighbours(1,ishell)%atomj_in_shell(iatshell)
      isym =Shell1at%neighbours(1,ishell)%sym_in_shell(iatshell)
      do ii=1,3
        Phi1_3(ii)=sum(Sym%S_ref(ii,:,isym,1)*Phi1_ref(:,ishell))
      end do  
      Phi1((iatom-1)*3+1:(iatom-1)*3+3)=Phi1_3(:)
    end do !iatshell
  end do !ishell
! Acoustic sum rule
  do ii=1,3
    do iatom=1,InVar%natom
      iat_mod=mod(iatom+InVar%natom_unitcell-1,InVar%natom_unitcell)+1
      if (iat_mod==1) cycle
      Phi1((iatom-iat_mod+1)*3+ii)=Phi1((iatom-iat_mod+1)*3+ii)-Phi1((iatom-1)*3+ii)
    end do
  end do  
  ABI_FREE(Phi1_3)
  ABI_FREE(Phi1_ref)

! Remove the rounding errors before writing (for non regression testing purposes)
  do ii=1,3*InVar%natom
    if (abs(Phi1(ii)).lt.tol8) Phi1(ii)=zero
  end do  

end subroutine tdep_calc_phi1 

!=====================================================================================================
subroutine tdep_write_phi1(InVar,Phi1)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  double precision,intent(in) :: Phi1(3*InVar%natom)

  integer :: iatcell,ii

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '#### For each shell, list of coefficients (IFC), number of neighbours... ####'
  write(InVar%stdout,*) '#############################################################################'

! Write the IFCs in the data.out file (with others specifications: 
! number of atoms in a shell, distance, Trace...)
  do iatcell=1,InVar%natom_unitcell
    write(InVar%stdout,'(a,i4)') ' ############# List of (first order) IFC for the reference atom=',iatcell
    write(InVar%stdout,'(2x,3(f9.6,1x))') (Phi1((iatcell-1)*3+ii),ii=1,3)
    write(InVar%stdout,*) ' '
  end do !iatcell  

end subroutine tdep_write_phi1 

!=====================================================================================================
subroutine tdep_calc_phi2(InVar,ntotcoeff,proj,Phi2_coeff,Phi2,Shell2at,Sym)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell2at
  integer,intent(in) :: ntotcoeff
  double precision,intent(in) :: Phi2_coeff(ntotcoeff,1),proj(9,9,Shell2at%nshell)
  double precision,intent(out) :: Phi2(3*InVar%natom,3*InVar%natom)

  integer :: ishell,isym,eatom,fatom,ncoeff,ncoeff_prev
  integer :: nshell,ii,jj,kk,ll,kappa,iatshell,itrans
  double precision,allocatable :: Phi2_33(:,:),Phi2_ref(:,:,:)

  nshell=Shell2at%nshell
  ABI_MALLOC(Phi2_ref,(3,3,nshell)); Phi2_ref(:,:,:)=zero
  ABI_MALLOC(Phi2_33,(3,3)) ; Phi2_33(:,:)=0.d0
  do ishell=1,nshell
!   Build the 3x3 IFC per shell    
    ncoeff     =Shell2at%ncoeff(ishell)
    ncoeff_prev=Shell2at%ncoeff_prev(ishell)
    kappa=0
    do ii=1,3
      do jj=1,3
        kappa=kappa+1
        Phi2_ref(ii,jj,ishell)=sum(proj(kappa,1:ncoeff,ishell)*Phi2_coeff(ncoeff_prev+1:ncoeff_prev+ncoeff,1))
      end do  
    end do  
    do eatom=1,Invar%natom
!     Build the 3x3 IFC of an atom in this shell    
      if (Shell2at%neighbours(eatom,ishell)%n_interactions.eq.0) cycle
      do iatshell=1,Shell2at%neighbours(eatom,ishell)%n_interactions
        fatom=Shell2at%neighbours(eatom,ishell)%atomj_in_shell(iatshell)
        isym =Shell2at%neighbours(eatom,ishell)%sym_in_shell(iatshell)
        itrans=Shell2at%neighbours(eatom,ishell)%transpose_in_shell(iatshell)
        if (fatom.lt.eatom) cycle
        call tdep_build_phi2_33(isym,Phi2_ref(:,:,ishell),Phi2_33,Sym,itrans) 
!       Symetrization of the Phi2 matrix
        Phi2((eatom-1)*3+1:(eatom-1)*3+3,3*(fatom-1)+1:3*(fatom-1)+3)=Phi2_33(:,:)
        do ii=1,3
          do jj=1,3
            Phi2((fatom  -1)*3+ii,3*(eatom-1)+jj)=Phi2_33(jj,ii)
          end do        
        end do  
      end do !iatshell
    end do !eatom
  end do !ishell
! Acoustic sum rule
  do eatom=1,InVar%natom
    do jj=1,3
      do kk=1,3
        do fatom=1,InVar%natom
          if (fatom==eatom) cycle
          Phi2((eatom-1)*3+jj,(eatom-1)*3+kk)=Phi2((eatom-1)*3+jj,3*(eatom-1)+kk)&
&                                               -Phi2((eatom-1)*3+jj,3*(fatom-1)+kk)
        enddo
      enddo
    enddo
  enddo
  ABI_FREE(Phi2_33)
  ABI_FREE(Phi2_ref)

! Remove the rounding errors before writing (for non regression testing purposes)
  do ii=1,3*InVar%natom
    do jj=1,3*InVar%natom
      if (abs(Phi2(ii,jj)).lt.tol8) Phi2(ii,jj)=zero
    end do
  end do  

end subroutine tdep_calc_phi2 

!=====================================================================================================
subroutine tdep_write_phi2(distance,InVar,MPIdata,Phi2,Shell2at)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(in) :: Shell2at
  type(MPI_enreg_type), intent(in) :: MPIdata
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision,intent(in) :: Phi2(3*InVar%natom,3*InVar%natom)

  integer :: iatcell,ishell,jshell,jatom
  integer :: nshell,ii,this_shell,iatshell
  double precision :: max_bound,min_bound,dist_shell
  integer,allocatable :: tab_shell(:)

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '#### For each shell, list of coefficients (IFC), number of neighbours... ####'
  write(InVar%stdout,*) '#############################################################################'

  nshell=Shell2at%nshell
! Write the IFCs in the data.out file (with others specifications: 
! number of atoms in a shell, distance, Trace...)
  ABI_MALLOC(tab_shell,(nshell)); tab_shell(:)=0
  do iatcell=1,InVar%natom_unitcell
    tab_shell(:)=0
    write(InVar%stdout,'(a,i4)') ' ############# List of (second order) IFC for the reference atom=',iatcell
!   Sort the IFC with distance in increasing order
    min_bound=-1.d0
    do ishell=1,nshell
      do jshell=1,nshell
        if ((distance(Shell2at%iatref(jshell),Shell2at%jatref(jshell),1).ge.min_bound).and.(tab_shell(jshell).eq.0)) then
          max_bound=distance(Shell2at%iatref(jshell),Shell2at%jatref(jshell),1)
          this_shell=jshell
        end if  
      end do

      do jshell=1,nshell
        if ((distance(Shell2at%iatref(jshell),Shell2at%jatref(jshell),1).lt.max_bound).and.&
&           (distance(Shell2at%iatref(jshell),Shell2at%jatref(jshell),1).ge.min_bound).and.&
&            (tab_shell(jshell).eq.0)) then
          max_bound=distance(Shell2at%iatref(jshell),Shell2at%jatref(jshell),1)
          this_shell=jshell
        end if
      end do
      tab_shell(this_shell)=1
      min_bound=max_bound
      dist_shell=distance(Shell2at%iatref(this_shell),Shell2at%jatref(this_shell),1)
      
!     Write the IFC properly  
      if (Shell2at%neighbours(iatcell,this_shell)%n_interactions.ne.0) then
        write(InVar%stdout,'(a,i4,a,i4,a,f9.6)') ' ======== NEW SHELL (ishell=',this_shell,&
&            '): There are',Shell2at%neighbours(iatcell,this_shell)%n_interactions,' atoms on this shell at distance=',dist_shell
        do iatshell=1,Shell2at%neighbours(iatcell,this_shell)%n_interactions
          jatom=Shell2at%neighbours(iatcell,this_shell)%atomj_in_shell(iatshell)
          write(InVar%stdout,'(a,i4,a,i4)') '  For jatom=',jatom,' ,with type=',mod(jatom-1,InVar%natom_unitcell)+1
          do ii=1,3
            write(InVar%stdout,'(2x,3(f9.6,1x))') Phi2((iatcell-1)*3+ii,(jatom-1)*3+1),Phi2((iatcell-1)*3+ii,(jatom-1)*3+2),&
&             Phi2((iatcell-1)*3+ii,(jatom-1)*3+3)
          end do
          write(InVar%stdout,'(a,3(1x,f11.6))') '  The components of the vector are:', distance(iatcell,jatom,2:4)
          write(InVar%stdout,'(a,(1x,f9.6))') '  Trace=',Phi2((iatcell-1)*3+1,(jatom-1)*3+1)+Phi2((iatcell-1)*3+2,&
&           (jatom-1)*3+2)+Phi2((iatcell-1)*3+3,(jatom-1)*3+3)
          write(InVar%stdout,*) ' '
        end do
      end if
    end do !ishell 
  end do !iatcell  
  ABI_FREE(tab_shell)

! Write the Phi2_unitcell.dat and Phi2.dat files
  if (InVar%debug.and.MPIdata%iam_master) then
    write(InVar%stdout,'(a)') ' See the Phi2*.dat file'
    open(unit=52,file=trim(InVar%output_prefix)//'Phi2_unitcell.dat')
    open(unit=55,file=trim(InVar%output_prefix)//'Phi2.dat')
    do jatom=1,3*InVar%natom
      if (jatom.le.3*InVar%natom_unitcell) then
        write(52,'(10000(f10.6,1x))') Phi2(jatom,:)
      end if  
      write(55,'(10000(f10.6,1x))') Phi2(jatom,:)
    end do  
    close(52)
    close(55)
  end if  

end subroutine tdep_write_phi2 

!=====================================================================================================
subroutine tdep_calc_dij(dij,eigenV,iqpt,InVar,Lattice,omega,Phi2,qpt_cart,Rlatt_cart)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  integer,intent(in) :: iqpt
  double precision,intent(in) :: Phi2(3*InVar%natom,3*InVar%natom)
  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)
  double precision,intent(in) :: qpt_cart(3)
  double precision,intent(out) :: omega (3*InVar%natom_unitcell)
  double complex  ,intent(out) :: dij   (3*InVar%natom_unitcell,3*InVar%natom_unitcell)
  double complex  ,intent(out) :: eigenV(3*InVar%natom_unitcell,3*InVar%natom_unitcell)

  integer :: LWORK,ii,jj,kk,iatom,jatom,iatcell,jatcell,itypat,jtypat,iat_mod,INFO,itemp,imode,nmode
  double precision :: phase
  double complex :: ctemp,norm
  double precision, allocatable :: RWORK(:)
  double complex, allocatable :: WORKC(:)
! double complex, allocatable :: mass_mat(:,:)

! Calculation of the dynamical matrix (Dij)
  do iatcell=1,InVar%natom_unitcell
    do jatom=1,InVar%natom
      iat_mod=mod(jatom+InVar%natom_unitcell-1,InVar%natom_unitcell)+1
      phase=0.d0
      do kk=1,3
        phase=phase+2*pi*Rlatt_cart(kk,iatcell,jatom)*qpt_cart(kk)
      end do
      do ii=1+(iatcell-1)*3,3+(iatcell-1)*3
        do jj=1,3
          dij(ii,3*(iat_mod-1)+jj)=dij(ii,3*(iat_mod-1)+jj)+dcmplx(Phi2(ii,3*(jatom-1)+jj),0.d0)*exp(dcmplx(0.d0,phase))
        end do !jj
      end do !ii 
    end do !jatom 
  end do !iatcell

! The Dij has to be an hermitian matrix
  itemp=0
  do ii=1,3*InVar%natom_unitcell
    do jj=ii,3*InVar%natom_unitcell
      if ((abs(real(dij(ii,jj))-real(dij(jj,ii))).gt.tol10).or.(abs(aimag(dij(ii,jj))+aimag(dij(jj,ii))).gt.tol10)) then
        if (InVar%debug) then
          write (InVar%stdout,'(a,1x,2(i4,1x))') 'for ii,jj=',ii,jj
          write (InVar%stdout,'(a,1x,1(f12.8,1x))') 'abs(realij-realji)=',abs(real(dij(ii,jj))-real(dij(jj,ii)))
          write (InVar%stdout,'(a,1x,1(f12.8,1x))') 'abs(imagij+imagji)=',abs(aimag(dij(ii,jj))+aimag(dij(jj,ii)))
        end if  
        itemp=itemp+1
      end if  
    end do
  end do
  if (itemp.ne.0.and.iqpt.eq.1) then
    write(InVar%stdout,*) 'WARNING: The Dij matrix is not hermitian'
    write(InVar%stdout,*) '  Probably: one shell may not have the whole number of atoms'
    write(InVar%stdout,*) '  The Dij matrix is symetrized'
  end if  

! Diagonalization of dynamical matrix Dij/sqrt(Mi*Mj)
  LWORK=2*3*InVar%natom_unitcell-1
  ABI_MALLOC(WORKC,(LWORK)); WORKC(:)=czero
  ABI_MALLOC(RWORK,(3*3*InVar%natom_unitcell-2)); RWORK(:)=zero
  do iatcell=1,InVar%natom_unitcell
    itypat=InVar%typat_unitcell(iatcell)
    do jatcell=1,InVar%natom_unitcell
      jtypat=InVar%typat_unitcell(jatcell)
      do ii=1,3
        do jj=1,3
          eigenV(ii+(iatcell-1)*3,jj+(jatcell-1)*3)=dij(ii+(iatcell-1)*3,jj+(jatcell-1)*3)/&
&                              dcmplx(dsqrt(InVar%amu(itypat)*InVar%amu(jtypat))*amu_emass,0.d0)          
        end do !jj
      end do !ii 
    end do !jatcell 
  end do !iatcell
  call ZHEEV('V','U',3*InVar%natom_unitcell,eigenV(:,:),3*InVar%natom_unitcell,omega(:),WORKC,LWORK,RWORK,INFO)

! Normalization of the eigenvectors
  nmode=3*InVar%natom_unitcell
  do imode=1,nmode
    norm=zero
    do iatom=1,InVar%natom_unitcell
      do ii=1,3
        norm=norm+eigenV(3*(iatom-1)+ii,imode)*conjg(eigenV(3*(iatom-1)+ii,imode))
      end do
    end do
    eigenV(:,imode)=eigenV(:,imode)/dsqrt(real(norm))
  end do

! Remove the squared-negative frequencies  
  do ii=1,InVar%natom_unitcell
    do jj=1,3
      if (omega((ii-1)*3+jj).lt.0.d0) then
        omega((ii-1)*3+jj)=-dsqrt(-omega((ii-1)*3+jj))
      else
        omega((ii-1)*3+jj)=dsqrt(omega((ii-1)*3+jj))
      end if
    end do  
  end do
  ABI_FREE(WORKC)
  ABI_FREE(RWORK)

end subroutine tdep_calc_dij

!=====================================================================================================
!FB subroutine tdep_write_dij(Eigen2nd,iqpt,InVar,Lattice,qpt_cart)
subroutine tdep_write_dij(Eigen2nd,iqpt,InVar,Lattice,qpt)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  integer,intent(in) :: iqpt
!FB  double precision,intent(in) :: qpt_cart(3)
  double precision,intent(in) :: qpt(3)
  type(Eigen_Variables_type),intent(in) :: Eigen2nd

  double precision, allocatable :: omega (:)
  double complex, allocatable   :: dij   (:,:)
  double complex, allocatable   :: eigenV(:,:)
  integer :: ii,jj,iatcell,jatcell

  ABI_MALLOC(omega ,(3*InVar%natom_unitcell))                       ; omega(:)   = zero
  ABI_MALLOC(dij   ,(3*InVar%natom_unitcell,3*InVar%natom_unitcell)); dij(:,:)   =czero
  ABI_MALLOC(eigenV,(3*InVar%natom_unitcell,3*InVar%natom_unitcell)); eigenV(:,:)=czero
  omega(:)=Eigen2nd%eigenval(:,iqpt)
  do iatcell=1,InVar%natom_unitcell
    do jatcell=1,InVar%natom_unitcell
      do ii=1,3
        do jj=1,3
          dij   ((iatcell-1)*3+ii,(jatcell-1)*3+jj)=dcmplx(Eigen2nd%dynmat  (1,ii,iatcell,jj,jatcell,iqpt),&
&                                                          Eigen2nd%dynmat  (2,ii,iatcell,jj,jatcell,iqpt))
          eigenV((iatcell-1)*3+ii,(jatcell-1)*3+jj)=dcmplx(Eigen2nd%eigenvec(1,ii,iatcell,jj,jatcell,iqpt),&
&                                                          Eigen2nd%eigenvec(2,ii,iatcell,jj,jatcell,iqpt))
        end do !ii 
      end do !jj 
    end do !jatcell 
  end do !iatcell

! Print the dynamical matrix (Dij)
!FB  write(52,'(a,1x,3(f10.6,1x))') 'For qpt=',qpt_cart(:)*Lattice%acell_unitcell(:)
  write(52,'(a,1x,3(f10.6,1x))') 'For qpt=',qpt(:)
  write(52,'(a,i4,a)') '  Dij(',iqpt,'real)='
  do iatcell=1,InVar%natom_unitcell
    write(52,'(100(f10.6,1x))') real(dij(1+(iatcell-1)*3,:))
    write(52,'(100(f10.6,1x))') real(dij(2+(iatcell-1)*3,:))
    write(52,'(100(f10.6,1x))') real(dij(3+(iatcell-1)*3,:))
  end do  
  write(52,'(a,i4,a)') '  Dij(',iqpt,'imag)='
  do iatcell=1,InVar%natom_unitcell
    write(52,'(100(f10.6,1x))') aimag(dij(1+(iatcell-1)*3,:))
    write(52,'(100(f10.6,1x))') aimag(dij(2+(iatcell-1)*3,:))
    write(52,'(100(f10.6,1x))') aimag(dij(3+(iatcell-1)*3,:))
  end do  
  write(52,*)' '

! Print the frequencies (omega)
  if (InVar%Enunit.eq.0) write(53,'(i5,1x,100(f15.6,1x))') iqpt,(omega(ii)*Ha_eV*1000,ii=1,3*InVar%natom_unitcell)
  if (InVar%Enunit.eq.1) write(53,'(i5,1x,100(f15.6,1x))') iqpt,(omega(ii)*Ha_cmm1   ,ii=1,3*InVar%natom_unitcell)
  if (InVar%Enunit.eq.2) write(53,'(i5,1x,100(f15.6,1x))') iqpt,(omega(ii)           ,ii=1,3*InVar%natom_unitcell)

! Print the eigenvectors (eigenV) 
  write(51,*) 'For iqpt=',iqpt
  do ii=1,3*InVar%natom_unitcell
    write(51,*) 'Mode number',ii,' energy',omega(ii)
    write(51,*) '  Real:'
    write(51,*) real(eigenV(:,ii))
    write(51,*) '  Imag:'
    write(51,*) aimag(eigenV(:,ii))
  end do  
  write(51,*) ' '
  ABI_FREE(omega)
  ABI_FREE(dij)
  ABI_FREE(eigenV)

end subroutine tdep_write_dij

!=====================================================================================================
subroutine tdep_build_phi2_33(isym,Phi2_ref,Phi2_33,Sym,itrans) 

  implicit none

  type(Symetries_Variables_type),intent(in) :: Sym
! type(Input_Variables_type),intent(in) :: InVar
  double precision, intent(in) :: Phi2_ref(3,3)
  double precision, intent(out) :: Phi2_33(3,3)
  integer,intent(in) :: isym,itrans

  double precision :: Phi2_tmp(3,3),tmp1(3,3)

! Transform in the new basis wrt S_ref
  call DGEMM('N','N',3,3,3,1.d0,Sym%S_ref(:,:,isym,1),3,Phi2_ref,3,0.d0,Phi2_tmp,3)
  call DGEMM('N','N',3,3,3,1.d0,Phi2_tmp,3,Sym%S_inv(:,:,isym,1),3,0.d0,Phi2_33,3)

  if ((itrans.lt.1).or.(itrans.gt.2)) then
    MSG_BUG('This value of the symmetry index is not permitted')
  end if
! Transpose the 3x3 matrix if required
  if (itrans.eq.2) then
    tmp1(:,:)=Phi2_33(:,:)
    Phi2_33(1,2)=tmp1(2,1)
    Phi2_33(1,3)=tmp1(3,1)
    Phi2_33(2,3)=tmp1(3,2)
    Phi2_33(2,1)=tmp1(1,2)
    Phi2_33(3,1)=tmp1(1,3)
    Phi2_33(3,2)=tmp1(2,3)
  end if  

end subroutine tdep_build_phi2_33

!=====================================================================================================
subroutine tdep_init_eigen2nd(Eigen2nd,natom_unitcell,nqpt)

  implicit none

  integer, intent(in) :: natom_unitcell,nqpt
  type(Eigen_Variables_type),intent(out) :: Eigen2nd

  ABI_MALLOC(Eigen2nd%eigenval,(3*natom_unitcell,nqpt));                    Eigen2nd%eigenval(:,:)        =zero
  ABI_MALLOC(Eigen2nd%eigenvec,(2,3,natom_unitcell,3,natom_unitcell,nqpt)); Eigen2nd%eigenvec(:,:,:,:,:,:)=zero
  ABI_MALLOC(Eigen2nd%dynmat  ,(2,3,natom_unitcell,3,natom_unitcell,nqpt)); Eigen2nd%dynmat  (:,:,:,:,:,:)=zero

end subroutine tdep_init_eigen2nd

!=====================================================================================================
subroutine tdep_destroy_eigen2nd(Eigen2nd)

  implicit none

  type(Eigen_Variables_type),intent(inout) :: Eigen2nd

  ABI_FREE(Eigen2nd%eigenval)
  ABI_FREE(Eigen2nd%eigenvec)
  ABI_FREE(Eigen2nd%dynmat)

end subroutine tdep_destroy_eigen2nd

!=====================================================================================================
subroutine tdep_write_yaml(Eigen2nd,Qpt,Prefix)

  implicit none

  type(Eigen_Variables_type),intent(in) :: Eigen2nd
  type(Qpoints_type),intent(in) :: Qpt
  character(len=*) :: Prefix
! type(Lattice_Variables_type),intent(in) :: Lattice

  integer :: ii,jj,iatcell,jatcell,iqpt,imode,nmode
  double precision :: distance
  double complex, allocatable   :: eigenV(:,:)
  
  nmode=size(Eigen2nd%eigenval,dim=1)
  open(unit=52,file=trim(Prefix)//'phonon-bands.yaml')
  write(52,'(a,i4)') 'nqpoint:',Qpt%nqpt 
  write(52,'(a,i4)') 'npath:',Qpt%qpt_tot-1
  write(52,'(a)')    'segment_nqpoint:'
  do ii=1,Qpt%qpt_tot-1
    write(52,'(a,i4)') '- ',Qpt%lgth_segments(ii)
  end do
  write(52,'(a,i4)')    'natom:',nmode/3
  write(52,'(a)')    'phonon:'
  distance=0.d0
  do iqpt=1,Qpt%nqpt
    write(52,'(a,3(f15.6,1x,a))') '- q-position: [',Qpt%qpt_red(1,iqpt),',',Qpt%qpt_red(2,iqpt),',',Qpt%qpt_red(3,iqpt),']'
    if (iqpt.gt.1) distance=distance+sqrt(sum((Qpt%qpt_cart(:,iqpt)-Qpt%qpt_cart(:,iqpt-1))**2))
    write(52,'(a,f15.6)') '  distance:',distance
    do ii=1,Qpt%qpt_tot
      if (sum(abs(Qpt%qpt_red(:,iqpt)-Qpt%special_red(ii,:))).lt.tol8) then
        write(52,'(3a)') "  label: '",trim(Qpt%special_qpt(ii)),"'"
        exit
      end if
    end do !ii
    ABI_MALLOC(eigenV,(nmode,nmode)); eigenV(:,:)=czero
    do iatcell=1,nmode/3
      do jatcell=1,nmode/3
        do ii=1,3
          do jj=1,3
            eigenV((iatcell-1)*3+ii,(jatcell-1)*3+jj)=dcmplx(Eigen2nd%eigenvec(1,ii,iatcell,jj,jatcell,iqpt),&
&                                                            Eigen2nd%eigenvec(2,ii,iatcell,jj,jatcell,iqpt))
          end do !ii 
        end do !jj 
      end do !jatcell 
    end do !iatcell
    write(52,'(a)')    '  band:'
    do imode=1,nmode
      write(52,'(a,i4)')    '  - #',imode
      write(52,'(a,f15.6)') '    frequency:',Eigen2nd%eigenval(imode,iqpt)*Ha_THz
      write(52,'(a)') '    eigenvector:'
      do iatcell=1,nmode/3
        write(52,'(a,i4)') "    - # atom ", iatcell
        do ii=1,3
          write(52,'(a,f18.9,a,f18.9,a)') "      - [",real(eigenV((iatcell-1)*3+ii,imode)),','&
&           ,aimag(eigenV((iatcell-1)*3+ii,imode)),']'
        end do
      end do
    end do !nmode  
    ABI_FREE(eigenV)
    write(52,'(a)') ''
  end do

  close(52)


end subroutine tdep_write_yaml

!=====================================================================================================
end module m_tdep_phi2
