#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_psij

  use defs_basis
  use m_crystal,          only : crystal_t
  use m_ddb,              only : ddb_type
  use m_ifc,              only : ifc_type,ifc_fourq
  use m_tdep_qpt,         only : Qpoints_type
  use m_tdep_readwrite,   only : Input_Variables_type
  use m_tdep_latt,        only : Lattice_Variables_type
  use m_tdep_shell,       only : Shell_Variables_type
  use m_tdep_sym,         only : Symetries_Variables_type
  use m_tdep_phij,        only : Eigen_Variables_type,Eigen_Variables_type,tdep_init_eigen2nd,tdep_destroy_eigen2nd

  implicit none

  public :: tdep_calc_psijfcoeff
  public :: tdep_build_psijNNN
  public :: tdep_build_psij333
  public :: tdep_calc_gruneisen
  public :: tdep_calc_alpha_gamma

contains

!====================================================================================================
subroutine tdep_calc_psijfcoeff(InVar,ntotcoeff,proj,Shell3at,Sym,ucart,fcoeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'tdep_calc_psijfcoeff'
!End of the abilint section

  implicit none

  integer, intent(in) :: ntotcoeff
  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell3at
  double precision, intent(in) :: ucart(3,InVar%natom,InVar%nstep)
  double precision, intent(in) :: proj(27,27,Shell3at%nshell)
  double precision, intent(out) :: fcoeff(3*InVar%natom*InVar%nstep,ntotcoeff)

  integer :: ishell,ncoeff,ncoeff_prev,istep,iatom,jatom,katom
  integer :: icoeff,isym,iatshell,trans
  integer :: mu,nu,xi,alpha,beta,gama
  double precision :: terme,temp
  double precision :: udiff_ki(3),udiff_ji(3),SSSu(3,27)
  double precision, allocatable :: SSS_ref(:,:,:,:,:,:)

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '########################## Compute the pseudo-inverse #######################'
  write(InVar%stdout,*) '#############################################################################'

! For each couple of atoms, transform the Psij (3x3) ifc matrix using the symetry operation (S)
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
        
  write(InVar%stdout,*) ' Compute the coefficients used in the Moore-Penrose...'
  do ishell=1,Shell3at%nshell
    do iatom=1,InVar%natom
      if (Shell3at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
      do iatshell=1,Shell3at%neighbours(iatom,ishell)%n_interactions
        jatom=Shell3at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
        katom=Shell3at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
!FB        if (iatom==katom.and.iatom==jatom) cycle
        if (iatom==katom) cycle
        isym =Shell3at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
        trans=Shell3at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
        ncoeff     =Shell3at%ncoeff(ishell)
        ncoeff_prev=Shell3at%ncoeff_prev(ishell)
  
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
!FB              SSSu(:,:)=SSSu(:,:)+SSS_ref(:,:,nu,xi,isym,trans)*udiff_ki(xi)*udiff_ji(nu)
!FB              SSSu(:,:)=SSSu(:,:)+SSS_ref(:,:,nu,xi,isym,trans)*(ucart(nu,jatom,istep)*ucart(xi,katom,istep)+ucart(nu,iatom,istep)*ucart(xi,iatom,istep))
              SSSu(:,:)=SSSu(:,:)+SSS_ref(:,:,nu,xi,isym,trans)*udiff_ki(xi)*ucart(nu,jatom,istep)
            end do  
          end do  
          do mu=1,3
            do icoeff=1,ncoeff
              terme=sum(SSSu(mu,:)*proj(:,icoeff,ishell))
              fcoeff(mu+3*(iatom-1)+3*InVar%natom*(istep-1),icoeff+ncoeff_prev)=fcoeff(mu+3*(iatom-1)+3*InVar%natom*(istep-1),icoeff+ncoeff_prev)+terme
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
subroutine tdep_build_psijNNN(distance,InVar,ntotcoeff,proj,Psij_coeff,Psij_NN,Shell3at,Sym)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'tdep_build_psijNNN'
!End of the abilint section

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell3at
  integer,intent(in) :: ntotcoeff
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4),proj(27,27,Shell3at%nshell)
  double precision,intent(in) :: Psij_coeff(ntotcoeff,1)
  double precision,intent(out) :: Psij_NN(3*InVar%natom,3*InVar%natom,3*InVar%natom)

  integer :: iatcell,ishell,isym,iatom,jatom,katom,eatom,fatom,gatom,ncoeff,ncoeff_prev
  integer :: iatref,jatref,katref,iatshell,nshell,trans
  integer :: ii,jj,kk,ll,kappa
  double precision :: norm
  double precision, allocatable :: Psij_333(:,:,:),Psij_ref(:,:,:,:)

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '#### For each shell, list of coefficients (IFC), number of neighbours... ####'
  write(InVar%stdout,*) '#############################################################################'

  nshell=Shell3at%nshell
!==========================================================================================
!======== 1/ Build the Psij_NN ============================================================
!==========================================================================================
  ABI_MALLOC(Psij_ref,(3,3,3,nshell)) ; Psij_ref(:,:,:,:)=0.d0
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
    do eatom=1,Invar%natom
!     Build the 3x3x3 IFC of an atom in this shell    
      if (Shell3at%neighbours(eatom,ishell)%n_interactions.eq.0) cycle
      do iatshell=1,Shell3at%neighbours(eatom,ishell)%n_interactions
        fatom=Shell3at%neighbours(eatom,ishell)%atomj_in_shell(iatshell)
        gatom=Shell3at%neighbours(eatom,ishell)%atomk_in_shell(iatshell)
        isym =Shell3at%neighbours(eatom,ishell)%sym_in_shell(iatshell)
        trans=Shell3at%neighbours(eatom,ishell)%transpose_in_shell(iatshell)
        if ((eatom.eq.fatom).and.(fatom.eq.gatom)) cycle
        call tdep_build_psij333(eatom,fatom,gatom,isym,InVar,Psij_ref(:,:,:,ishell),Psij_333,Sym,trans) 
!       Symetrization of the Psij_NN matrix
        do ii=1,3
          do jj=1,3
            do kk=1,3
              Psij_NN(3*(eatom-1)+ii,3*(fatom-1)+jj,3*(gatom-1)+kk)=Psij_333(ii,jj,kk) !\Psi_efg
              Psij_NN(3*(eatom-1)+ii,3*(gatom-1)+kk,3*(fatom-1)+jj)=Psij_333(ii,jj,kk) !\Psi_egf
              Psij_NN(3*(fatom-1)+jj,3*(eatom-1)+ii,3*(gatom-1)+kk)=Psij_333(ii,jj,kk) !\Psi_feg
              Psij_NN(3*(fatom-1)+jj,3*(gatom-1)+kk,3*(eatom-1)+ii)=Psij_333(ii,jj,kk) !\Psi_fge
              Psij_NN(3*(gatom-1)+kk,3*(eatom-1)+ii,3*(fatom-1)+jj)=Psij_333(ii,jj,kk) !\Psi_gef
              Psij_NN(3*(gatom-1)+kk,3*(fatom-1)+jj,3*(eatom-1)+ii)=Psij_333(ii,jj,kk) !\Psi_gfe
            end do  
          end do        
        end do  
      end do !iatshell
    end do !eatom
  end do !ishell

! Set the (efe) & (eef) terms to zero. These terms will be initialized with the accoustic sum rule
  do eatom=1,InVar%natom
    do fatom=1,InVar%natom
      do ii=1,3
        do jj=1,3
          do kk=1,3
            Psij_NN(3*(eatom-1)+ii,3*(fatom-1)+jj,3*(eatom-1)+kk)=zero
            Psij_NN(3*(eatom-1)+ii,3*(eatom-1)+jj,3*(fatom-1)+kk)=zero
          end do
        end do
      end do
    end do
  end do
  do eatom=1,InVar%natom
    do fatom=1,InVar%natom
      do gatom=1,InVar%natom
        do ii=1,3
          do jj=1,3
!FB            do kk=1,3
!FB              if (fatom.ne.eatom.and.gatom.ne.eatom) then
!FB                Psij_NN(3*(eatom-1)+ii,3*(fatom-1)+jj,3*(eatom-1)+kk)=&
!FB&               Psij_NN(3*(eatom-1)+ii,3*(fatom-1)+jj,3*(eatom-1)+kk)-&
!FB&               Psij_NN(3*(eatom-1)+ii,3*(fatom-1)+jj,3*(gatom-1)+kk)
!FB
!FB                Psij_NN(3*(eatom-1)+ii,3*(eatom-1)+jj,3*(gatom-1)+kk)=&
!FB&               Psij_NN(3*(eatom-1)+ii,3*(eatom-1)+jj,3*(gatom-1)+kk)-&
!FB&               Psij_NN(3*(eatom-1)+ii,3*(fatom-1)+jj,3*(gatom-1)+kk)
!FB              end if
            do kk=1,3
              if (gatom.eq.eatom) cycle
              Psij_NN(3*(eatom-1)+ii,3*(fatom-1)+jj,3*(eatom-1)+kk)=&
&             Psij_NN(3*(eatom-1)+ii,3*(fatom-1)+jj,3*(eatom-1)+kk)-&
&             Psij_NN(3*(eatom-1)+ii,3*(fatom-1)+jj,3*(gatom-1)+kk)
            end do 
            do kk=1,3
              if (fatom.eq.eatom) cycle
              Psij_NN(3*(eatom-1)+ii,3*(eatom-1)+jj,3*(gatom-1)+kk)=&
&             Psij_NN(3*(eatom-1)+ii,3*(eatom-1)+jj,3*(gatom-1)+kk)-&
&             Psij_NN(3*(eatom-1)+ii,3*(fatom-1)+jj,3*(gatom-1)+kk)
            end do  
          end do
        end do
      end do
    end do
  end do

!FB! Accoustic sum ruleS (whatever eatom and fatom)
!FB  do ishell=1,Shell3at%nshell
!FB    do eatom=1,InVar%natom
!FB      if (Shell3at%neighbours(eatom,ishell)%n_interactions.eq.0) cycle
!FB      do iatshell=1,Shell3at%neighbours(eatom,ishell)%n_interactions
!FB        fatom=Shell3at%neighbours(eatom,ishell)%atomj_in_shell(iatshell)
!FB        gatom=Shell3at%neighbours(eatom,ishell)%atomk_in_shell(iatshell)
!FB        if (gatom==eatom) cycle
!FB        if (fatom==eatom) cycle
!FB        do ii=1,3
!FB          do jj=1,3
!FB            do kk=1,3
!FB              Psij_NN(3*(eatom-1)+ii,3*(fatom-1)+jj,3*(eatom-1)+kk)=&
!FB&             Psij_NN(3*(eatom-1)+ii,3*(fatom-1)+jj,3*(eatom-1)+kk)-&
!FB&             Psij_NN(3*(eatom-1)+ii,3*(fatom-1)+jj,3*(gatom-1)+kk)
!FB
!FB              Psij_NN(3*(eatom-1)+ii,3*(eatom-1)+jj,3*(gatom-1)+kk)=&
!FB&             Psij_NN(3*(eatom-1)+ii,3*(eatom-1)+jj,3*(gatom-1)+kk)-&
!FB&             Psij_NN(3*(eatom-1)+ii,3*(fatom-1)+jj,3*(gatom-1)+kk)
!FB            enddo !kk
!FB          enddo !jj
!FB        enddo !ii
!FB      enddo !iatshell
!FB    end do !eatom
!FB  end do !ishell

! Set the (eee) terms to zero. These terms will be initialized with the accoustic sum rule
  do eatom=1,InVar%natom
    do ii=1,3
      do jj=1,3
        do kk=1,3
          Psij_NN(3*(eatom-1)+ii,3*(eatom-1)+jj,3*(eatom-1)+kk)=zero
        end do
      end do
    end do
  end do
!FB  do ishell=1,Shell3at%nshell
!FB    do eatom=1,InVar%natom
!FB      if (Shell3at%neighbours(eatom,ishell)%n_interactions.eq.0) cycle
!FB      do iatshell=1,Shell3at%neighbours(eatom,ishell)%n_interactions
!FB        fatom=Shell3at%neighbours(eatom,ishell)%atomj_in_shell(iatshell)
!FB        gatom=Shell3at%neighbours(eatom,ishell)%atomk_in_shell(iatshell)
  do eatom=1,InVar%natom
    do gatom=1,InVar%natom
      if (gatom.eq.eatom) cycle
      do ii=1,3
        do jj=1,3
          do kk=1,3
            Psij_NN(3*(eatom-1)+ii,3*(eatom-1)+jj,3*(eatom-1)+kk)=&
&           Psij_NN(3*(eatom-1)+ii,3*(eatom-1)+jj,3*(eatom-1)+kk)-&
&           Psij_NN(3*(eatom-1)+ii,3*(eatom-1)+jj,3*(gatom-1)+kk)
          enddo !kk
        enddo !jj
      enddo !ii
    end do !gatom
  end do !eatom

  ABI_FREE(Psij_333)
  ABI_FREE(Psij_ref)

! Check the acoustic sum rules
  do eatom=1,InVar%natom
    do fatom=1,InVar%natom
      do ii=1,3
        do jj=1,3
          do kk=1,3
            norm=zero
            do gatom=1,InVar%natom
              norm=norm+Psij_NN(3*(eatom-1)+ii,3*(fatom-1)+jj,3*(gatom-1)+kk)
            end do !gatom
            if (abs(norm).gt.tol8) then
              write(std_err,*) ' BUG : the acoustic sum rule is not fulfilled'
              write(std_err,'(5(i3,x),1(e17.10,x))') ii,jj,kk,eatom,fatom,norm
              stop -1
            end if
          end do !kk
        end do !jj
      end do !ii  
    end do !fatom  
  end do !eatom

!==========================================================================================
!======== 2/ Write the Psij_NN in output ==================================================
!==========================================================================================
! Remove the rounding errors before writing (for non regression testing purposes)
  do ii=1,3*InVar%natom
    do jj=1,3*InVar%natom
      do kk=1,3*InVar%natom
        if (abs(Psij_NN(ii,jj,kk)).lt.tol8) Psij_NN(ii,jj,kk)=zero
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
        write(InVar%stdout,'(a,i4,a,i4)') '  For iatcell=',iatref,' ,with type=',mod(iatref-1,InVar%natom_unitcell)+1
        write(InVar%stdout,'(a,i4,a,i4)') '  For jatom  =',jatom ,' ,with type=',mod(jatom -1,InVar%natom_unitcell)+1
        write(InVar%stdout,'(a,i4,a,i4)') '  For katom  =',katom ,' ,with type=',mod(katom -1,InVar%natom_unitcell)+1
        write(InVar%stdout,'(a)') '  \Psi^{\alpha\beta x}='
        do ii=1,3
          write(InVar%stdout,'(2x,3(f9.6,x))') (Psij_NN((iatref-1)*3+ii,(jatom-1)*3+jj,(katom-1)*3+1),jj=1,3)
        end do
        write(InVar%stdout,'(a)') '  \Psi^{\alpha\beta y}='
        do ii=1,3
          write(InVar%stdout,'(2x,3(f9.6,x))') (Psij_NN((iatref-1)*3+ii,(jatom-1)*3+jj,(katom-1)*3+2),jj=1,3)
        end do
        write(InVar%stdout,'(a)') '  \Psi^{\alpha\beta z}='
        do ii=1,3
          write(InVar%stdout,'(2x,3(f9.6,x))') (Psij_NN((iatref-1)*3+ii,(jatom-1)*3+jj,(katom-1)*3+3),jj=1,3)
        end do
        write(InVar%stdout,'(a,3(x,f9.6))') '  (i,j) vector components:', distance(iatref,jatom,2:4)
        write(InVar%stdout,'(a,3(x,f9.6))') '  (j,k) vector components:', distance(jatom ,katom,2:4)
        write(InVar%stdout,'(a,3(x,f9.6))') '  (k,i) vector components:', distance(katom,iatref,2:4)
        write(InVar%stdout,*) ' '
      end do !iatshell
    end if !n_interactions 
  end do !ishell 

! Write the Psij_unitcell.dat and Psij_NN.dat files
  if (InVar%debug) then
    write(InVar%stdout,'(a)') ' See the Psij*.dat file'
    open(unit=52,file='Psij_unitcell.dat')
    open(unit=55,file='Psij_NN.dat')
    do iatcell=1,3*InVar%natom
      if (iatcell.le.3*InVar%natom_unitcell) then
        write(52,'(10000(f10.6,x))') Psij_NN(iatcell,:,1),Psij_NN(iatcell,:,2),Psij_NN(iatcell,:,3)
      end if  
      write(55,'(10000(f10.6,x))') Psij_NN(iatcell,:,1),Psij_NN(iatcell,:,2),Psij_NN(iatcell,:,3)
    end do  
    close(52)
    close(55)
  end if  

end subroutine tdep_build_psijNNN

!=====================================================================================================
subroutine tdep_calc_gruneisen(distance,Eigen2nd,Gruneisen,iqpt,InVar,Lattice,Psij_NN,qpt_cart,Rlatt_cart,Shell3at)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'tdep_calc_gruneisen'
!End of the abilint section

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(in) :: Shell3at
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Eigen_Variables_type),intent(in) :: Eigen2nd
  integer,intent(in) :: iqpt
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision,intent(in) :: Psij_NN(3*InVar%natom,3*InVar%natom,3*InVar%natom)
  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)
  double precision,intent(in) :: qpt_cart(3)
  double complex  ,intent(out) :: Gruneisen(3*InVar%natom_unitcell)

  integer :: ii,jj,kk,iatcell,jatom,katom,itypat,jtypat,natom_unitcell
  integer :: imode,nmode,ncomp,jat_mod,jatcell,iatshell,ishell
  double precision :: phase
  double complex, allocatable :: mass_mat(:,:),eigen_prod(:,:,:,:,:),Grun_shell(:,:)

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
  Gruneisen(:)=zero
  do ishell=1,Shell3at%nshell
    do iatcell=1,InVar%natom_unitcell
      itypat=InVar%typat_unitcell(iatcell)
      do iatshell=1,Shell3at%neighbours(iatcell,ishell)%n_interactions
        jatom=Shell3at%neighbours(iatcell,ishell)%atomj_in_shell(iatshell)
        katom=Shell3at%neighbours(iatcell,ishell)%atomk_in_shell(iatshell)
        jat_mod=mod(jatom+InVar%natom_unitcell-1,InVar%natom_unitcell)+1
        jtypat=InVar%typat_unitcell(jat_mod)
        phase=0.d0
        do jj=1,3
          phase=phase+2*pi*Rlatt_cart(jj,iatcell,jatom)*qpt_cart(jj)
        end do
        do ii=1,3
          do jj=1,3
            do kk=1,3
              do imode=1,nmode
                Grun_shell(imode,ishell)=Grun_shell(imode,ishell)-dcmplx(Psij_NN(3*(iatcell-1)+ii,3*(jatom-1)+jj,3*(katom-1)+kk),0.d0)*&
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
  ABI_FREE(eigen_prod)
  do ishell=1,Shell3at%nshell
    Gruneisen(:)=Gruneisen(:)+Grun_shell(:,ishell)
  end do  
  ABI_FREE(Grun_shell)

end subroutine tdep_calc_gruneisen

!=====================================================================================================
subroutine tdep_build_psij333(eatom,fatom,gatom,isym,InVar,Psij_ref,Psij_333,Sym,trans) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'tdep_build_psij333'
!End of the abilint section

  implicit none

  type(Symetries_Variables_type),intent(in) :: Sym
  type(Input_Variables_type),intent(in) :: InVar
  double precision, intent(in) :: Psij_ref(3,3,3)
  double precision, intent(out) :: Psij_333(3,3,3)
  integer,intent(in) :: eatom,fatom,gatom,isym,trans

  integer :: alpha,beta,gama
  integer :: ii,jj,kk,ee,ff,gg,mu,nu,xi
  double precision :: Psij_tmp(3,3,3)


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
    write(InVar%stdout,'(a)') '  BUG: this value of the symmetry index is not permitted'
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
subroutine tdep_calc_alpha_gamma(Crystal,distance,DDB,Ifc,InVar,Lattice,Psij_NN,Rlatt_cart,Shell3at,Sym)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'tdep_calc_alpha_gamma'
 use interfaces_41_geometry
!End of the abilint section

  implicit none

  type(crystal_t),intent(in) :: Crystal
  type(ddb_type),intent(in) :: DDB
  type(ifc_type),intent(in) :: Ifc
  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Shell_Variables_type),intent(in) :: Shell3at
  type(Symetries_Variables_type),intent(in) :: Sym
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision,intent(in) :: Psij_NN(3*InVar%natom,3*InVar%natom,3*InVar%natom)
  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)

  integer :: iq_ibz,nqbz,nqibz,ii,iatcell,jatcell,jj,nmode
  integer, allocatable :: ibz2bz(:)
  double precision :: k_B,wovert,xx,C_v,Gama,alpha_v
  double precision, allocatable :: qbz(:,:),wtq(:),wtq_folded(:),wtqibz(:),qibz(:,:),qibz_cart(:,:)
  double precision, allocatable :: omega(:,:),displ(:,:),out_eigvec(:,:,:,:,:,:),heatcapa(:),grun_thermo(:)
  double complex :: Gruneisen(3*InVar%natom_unitcell)
  type(Eigen_Variables_type) :: Eigen2nd_tmp

! Initialize useful quantities stored in the DDB file  
  nqbz=DDB%nblok
  ABI_MALLOC(qbz,(3,nqbz)); qbz(:,:)=zero
  qbz(:,:)=DDB%qpt(1:3,:)

! Reduce the number of such points by symmetrization.
  ABI_ALLOCATE(ibz2bz,(nqbz))
  ABI_ALLOCATE(wtq,(nqbz))
  ABI_ALLOCATE(wtq_folded,(nqbz))
  wtq(:)=one/nqbz         ! Weights sum up to one

!FB  write(InVar%stdlog,*) 'nqbz = ', nqbz
  call symkpt(0,Crystal%gmet,ibz2bz,InVar%stdlog,qbz,nqbz,nqibz,Sym%nsym,Crystal%symrec,1,wtq,wtq_folded)
!FB  write(InVar%stdlog,*) 'nqibz = ', nqibz
 
  ABI_ALLOCATE(wtqibz   ,(nqibz))
  ABI_ALLOCATE(qibz     ,(3,nqibz))
  ABI_ALLOCATE(qibz_cart,(3,nqibz))
  do iq_ibz=1,nqibz
    wtqibz(iq_ibz)=wtq_folded(ibz2bz(iq_ibz))
    qibz(:,iq_ibz)=qbz(:,ibz2bz(iq_ibz))
  end do
  ABI_DEALLOCATE(wtq_folded)

! Loop over irreducible q-points
! =======================
  nmode=3*InVar%natom_unitcell
  k_B=kb_HaK*Ha_eV
  wovert=1.d0/(2*InVar%temperature*k_B)
  ABI_MALLOC(heatcapa,(nmode))   ; heatcapa   (:)=0.d0 
  ABI_MALLOC(grun_thermo,(nmode)); grun_thermo(:)=0.d0 
  ABI_MALLOC(displ,(2*3*InVar%natom_unitcell*3*InVar%natom_unitcell,nqibz)); displ(:,:)=zero
  ABI_MALLOC(out_eigvec,(2,3,InVar%natom_unitcell,3,InVar%natom_unitcell,nqibz)); out_eigvec(:,:,:,:,:,:)=zero
  ABI_MALLOC(omega,(3*InVar%natom_unitcell,nqibz)); omega(:,:)=zero
  call tdep_init_eigen2nd(Eigen2nd_tmp,InVar%natom_unitcell,nqibz)
  do iq_ibz=1,nqibz
!   Compute the frequencies
!   =======================
    call ifc_fourq(Ifc,Crystal,qibz(:,iq_ibz),omega(:,iq_ibz),displ(:,iq_ibz),out_eigvec=out_eigvec(:,:,:,:,:,iq_ibz))
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
!   Print the frequencies (omega)
!FB    write(InVar%stdlog,'(i5,x,100(f15.6,x))') iq_ibz,(omega(ii,iq_ibz)*Ha_eV*1000,ii=1,nmode)
!FB    write(InVar%stdlog,*) 'wtq=',wtqibz(iq_ibz)
!   Print the eigenvectors (eigenV) 
!FB    do ii=1,nmode
!FB      write(InVar%stdlog,*) 'Mode number',ii,' energy',omega(ii,iq_ibz)
!FB      write(InVar%stdlog,*) '  Real:'
!FB      write(InVar%stdlog,*) real(Eigen2nd_tmp%eigenvec(:,ii,iq_ibz))
!FB      write(InVar%stdlog,*) '  Imag:'
!FB      write(InVar%stdlog,*) imag(Eigen2nd_tmp%eigenvec(:,ii,iq_ibz))
!FB    end do  
    call tdep_calc_gruneisen(distance,Eigen2nd_tmp,Gruneisen,iq_ibz,InVar,Lattice,Psij_NN,qibz_cart(:,iq_ibz),Rlatt_cart,Shell3at)
!FB    write(InVar%stdlog,*) 'Gruneisen='
!FB    write(InVar%stdlog,'(i5,x,100(e15.6,x))') iq_ibz,(real(Gruneisen(ii)),ii=1,nmode)

!   Compute the heat capacity and thermodynamical gruneisen parameter
!   =================================================================
    do ii=1,nmode
      xx=omega(ii,iq_ibz)*Ha_eV
      if (xx.le.0) cycle
!FB      write(6,*)'wovert,xx=',wovert,xx
      heatcapa(ii)   =heatcapa(ii)   +(wovert*xx/sinh(wovert*xx))**2*wtqibz(iq_ibz)
      grun_thermo(ii)=grun_thermo(ii)+(wovert*xx/sinh(wovert*xx))**2*Gruneisen(ii)*wtqibz(iq_ibz)
!FB      write(InVar%stdlog,'(a,i5,x,100(f15.6,x))')'mode,heatcapa(mode)='   ,ii,heatcapa(ii)
!FB      write(InVar%stdlog,'(a,i5,x,100(f15.6,x))')'mode,grun_thermo(mode)=',ii,grun_thermo(ii)
    end do  
  end do
  call tdep_destroy_eigen2nd(Eigen2nd_tmp)
  ABI_DEALLOCATE(displ)
  ABI_DEALLOCATE(out_eigvec)
  ABI_DEALLOCATE(omega)
  ABI_DEALLOCATE(ibz2bz)
  ABI_DEALLOCATE(wtq)
  ABI_DEALLOCATE(qbz)
  ABI_DEALLOCATE(wtqibz)
  ABI_DEALLOCATE(qibz)
  ABI_DEALLOCATE(qibz_cart)

! Summary
! =======
  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '################## Gruneisen parameter, Thermal expansion... ################'
  write(InVar%stdout,*) '#############################################################################'
!FB  heatcapa   (:)=heatcapa   (:)/InVar%natom_unitcell
!FB  write(InVar%stdlog,*) 'Summary:'
!FB  write(InVar%stdlog,'(a,x,100(f15.6,x))') 'C_v(per_mode \& per unit cell)=',(heatcapa(ii)   ,ii=1,nmode)
!FB  write(InVar%stdlog,'(a,x,100(f15.6,x))') 'Gamma(per mode)=',(grun_thermo(ii),ii=1,nmode)
  C_v =sum(heatcapa(:))
  Gama=sum(grun_thermo(:))/C_v
  alpha_v=Gama*C_v*kb_HaK*Ha_J/(Lattice%BulkModulus*1.d9)/(Lattice%ucvol*Bohr_Ang**3*1.d-30)
  write(InVar%stdout,*) 'Specific heat : C_v  (in k_B per unit cell)=',C_v
  write(InVar%stdout,*) 'Gruneisen parameter : Gamma=',Gama
  write(InVar%stdout,*) 'Bulk Modulus (in GPa): B=',Lattice%BulkModulus
  write(InVar%stdout,*) 'Volume (in bohr^3 per unit cell): V=',Lattice%ucvol
  write(InVar%stdout,*) 'Thermal expansion : alpha_v=',alpha_v
  ABI_FREE(heatcapa)

!FB  0.473294364993209*(2.38255605878933*1.38065e-23)/1.20512e11/3.0/((30.6135754000*0.529177e-10)**3/216)

end subroutine tdep_calc_alpha_gamma

!=====================================================================================================
subroutine tdep_write_gruneisen(distance,Eigen2nd,InVar,Lattice,Psij_NN,Qpt,Rlatt_cart,Shell3at)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'tdep_write_gruneisen'
!End of the abilint section

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Shell_Variables_type),intent(in) :: Shell3at
  type(Lattice_Variables_type),intent(in) :: Lattice
  type(Eigen_Variables_type),intent(in) :: Eigen2nd
  type(Qpoints_type),intent(in) :: Qpt
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision,intent(in) :: Psij_NN(3*InVar%natom,3*InVar%natom,3*InVar%natom)
  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)

  integer :: iqpt,nmode,ii,jj
  double precision :: qpt_cart(3)
  double complex, allocatable :: Gruneisen(:)
  
  nmode=3*InVar%natom_unitcell
  ABI_MALLOC(Gruneisen,(3*InVar%natom_unitcell)); Gruneisen(:)=zero
  open(unit=53,file='gruneisen.dat')
  do iqpt=1,Qpt%nqpt
    qpt_cart(:)=Qpt%qpt_cart(:,iqpt)
    call tdep_calc_gruneisen(distance,Eigen2nd,Gruneisen,iqpt,InVar,Lattice,Psij_NN,qpt_cart,Rlatt_cart,Shell3at)
!   Write the Gruneisen
    if (sum(abs(qpt_cart(:))).gt.tol8) then 
      if (sum(abs(imag(Gruneisen(:)))).gt.tol8) then
        write(6,*) 'BUG : the imaginary part of the Gruneisen is not equal to zero'
        write(53,'(i5,x,100(e15.6,x))') iqpt,(real(Gruneisen(ii)),ii=1,nmode),(imag(Gruneisen(ii)),ii=1,nmode)
      else 
!FB        write(53,'(i5,x,500(e15.6,x))') iqpt,(real(Gruneisen(ii)),ii=1,nmode),((real(Grun_shell(ii,jj)),ii=1,nmode),jj=1,Shell3at%nshell)
        write(53,'(i5,x,500(e15.6,x))') iqpt,(real(Gruneisen(ii)),ii=1,nmode)
      end if  
    end if  
  end do  
  close(53)
  ABI_FREE(Gruneisen)

end subroutine tdep_write_gruneisen

!=====================================================================================================
end module m_tdep_psij
