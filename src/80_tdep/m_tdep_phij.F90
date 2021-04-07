
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_phij

  use defs_basis
  use m_errors
  use m_abicore
  use m_tdep_readwrite,   only : Input_Variables_type
  use m_tdep_latt,        only : Lattice_Variables_type
  use m_tdep_shell,       only : Shell_Variables_type
  use m_tdep_sym,         only : Symetries_Variables_type
  use m_tdep_qpt,         only : Qpoints_type
  use m_tdep_utils,       only : Coeff_Moore_type

  implicit none

  type Eigen_Variables_type

    double precision, allocatable :: eigenval(:,:)
    double complex, allocatable :: eigenvec(:,:,:)

  end type Eigen_Variables_type

  public :: tdep_calc_pijfcoeff
  public :: tdep_build_pijN
  public :: tdep_calc_phijfcoeff
  public :: tdep_build_phijNN
  public :: tdep_build_phij33
  public :: tdep_calc_dij
  public :: tdep_write_dij
  public :: tdep_init_eigen2nd
  public :: tdep_destroy_eigen2nd
  public :: tdep_write_yaml

contains

!====================================================================================================
subroutine tdep_calc_pijfcoeff(CoeffMoore,InVar,proj,Shell1at,Sym)

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
          do istep=1,InVar%nstep
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

end subroutine tdep_calc_pijfcoeff

!====================================================================================================
subroutine tdep_calc_phijfcoeff(CoeffMoore,InVar,proj,Shell2at,Sym,ucart)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell2at
  type(Coeff_Moore_type), intent(inout) :: CoeffMoore
  double precision, intent(in) :: ucart(3,InVar%natom,InVar%nstep)
  double precision, intent(in) :: proj(9,9,Shell2at%nshell)

  integer :: ishell,ncoeff,ncoeff_prev,istep,iatom,jatom,iatshell
  integer :: icoeff,isym
  integer :: mu,nu,alpha,beta,trans
  double precision :: terme,temp
  double precision :: udiff(3),SSu(3,9)
  double precision, allocatable :: SS_ref(:,:,:,:,:)

! For each couple of atoms, transform the Phij (3x3) ifc matrix using the symetry operation (S)
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
        trans=Shell2at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
        ncoeff     =Shell2at%ncoeff(ishell)
        ncoeff_prev=Shell2at%ncoeff_prev(ishell)+CoeffMoore%ncoeff1st

        do istep=1,InVar%nstep

!         In order to impose the acoustic sum rule we use (u(j)-u(i))==u_j^\nu
          udiff(1)=ucart(1,jatom,istep)-ucart(1,iatom,istep)
          udiff(2)=ucart(2,jatom,istep)-ucart(2,iatom,istep)
          udiff(3)=ucart(3,jatom,istep)-ucart(3,iatom,istep)
            
!         F_i^\mu(t)=\sum_{\alpha\beta,j,\nu}S^{\mu\alpha}.S^{\nu\beta}.\Phi_{ij}^{\alpha\beta}.u_j^\nu(t)
          SSu(:,:)=zero
          do nu=1,3
            SSu(:,:)=SSu(:,:)+SS_ref(:,:,nu,isym,trans)*udiff(nu)
          end do  
          do mu=1,3
            do icoeff=1,ncoeff
              terme=sum(SSu(mu,:)*proj(:,icoeff,ishell))
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

end subroutine tdep_calc_phijfcoeff

!=====================================================================================================
subroutine tdep_build_pijN(InVar,ntotcoeff,proj,Pij_coeff,Pij_N,Shell1at,Sym)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell1at
  integer,intent(in) :: ntotcoeff
  double precision,intent(in) :: proj(3,3,Shell1at%nshell)
  double precision,intent(in) :: Pij_coeff(ntotcoeff,1)
  double precision,intent(out) :: Pij_N(3*InVar%natom)

  integer :: iatcell,ishell,isym,iatom,ncoeff,ncoeff_prev
  integer :: nshell,ii,iatshell,iat_mod
  double precision :: sum1
  double precision,allocatable :: Pij_3(:),Pij_ref(:,:)

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '#### For each shell, list of coefficients (IFC), number of neighbours... ####'
  write(InVar%stdout,*) '#############################################################################'

  nshell=Shell1at%nshell
!==========================================================================================
!======== 1/ Build the Pij_N ============================================================
!==========================================================================================
  ABI_MALLOC(Pij_ref,(3,nshell)); Pij_ref(:,:)=zero
  ABI_MALLOC(Pij_3,(3)) ; Pij_3(:)=0.d0
  do ishell=1,nshell
!   Build the 3x3 IFC per shell    
    ncoeff     =Shell1at%ncoeff(ishell)
    ncoeff_prev=Shell1at%ncoeff_prev(ishell)
    do ii=1,3
      Pij_ref(ii,ishell)=sum(proj(ii,1:ncoeff,ishell)*Pij_coeff(ncoeff_prev+1:ncoeff_prev+ncoeff,1))
    end do  
!   Build the vector-IFC of an atom in this shell    
    if (Shell1at%neighbours(1,ishell)%n_interactions.eq.0) cycle
    do iatshell=1,Shell1at%neighbours(1,ishell)%n_interactions
      iatom=Shell1at%neighbours(1,ishell)%atomj_in_shell(iatshell)
      isym =Shell1at%neighbours(1,ishell)%sym_in_shell(iatshell)
      do ii=1,3
        Pij_3(ii)=sum(Sym%S_ref(ii,:,isym,1)*Pij_ref(:,ishell))
      end do  
      Pij_N((iatom-1)*3+1:(iatom-1)*3+3)=Pij_3(:)
    end do !iatshell
  end do !ishell
! Acoustic sum rule
  do ii=1,3
    do iatom=1,InVar%natom
      iat_mod=mod(iatom+InVar%natom_unitcell-1,InVar%natom_unitcell)+1
      if (iat_mod==1) cycle
      Pij_N((iatom-iat_mod+1)*3+ii)=Pij_N((iatom-iat_mod+1)*3+ii)-Pij_N((iatom-1)*3+ii)
    end do
  end do  
  
  do ii=1,3
    sum1=zero
    do iatom=1,InVar%natom_unitcell
      sum1=sum1+Pij_N((iatom-1)*3+ii)
    enddo
    if (sum1.gt.tol8) then
      ABI_WARNING('The acoustic sum rule is not fulfilled at the 1st order')
    end if
  enddo
  ABI_FREE(Pij_3)
  ABI_FREE(Pij_ref)

!==========================================================================================
!======== 2/ Write the Pij_N in output ==================================================
!==========================================================================================
! Remove the rounding errors before writing (for non regression testing purposes)
  do ii=1,3*InVar%natom
    if (abs(Pij_N(ii)).lt.tol8) Pij_N(ii)=zero
  end do  

! Write the IFCs in the data.out file (with others specifications: 
! number of atoms in a shell, distance, Trace...)
  do iatcell=1,InVar%natom_unitcell
    write(InVar%stdout,'(a,i4)') ' ############# List of (first order) IFC for the reference atom=',iatcell
    write(InVar%stdout,'(2x,3(f9.6,1x))') (Pij_N((iatcell-1)*3+ii),ii=1,3)
    write(InVar%stdout,*) ' '
  end do !iatcell  

end subroutine tdep_build_pijN 

!=====================================================================================================
subroutine tdep_build_phijNN(distance,InVar,ntotcoeff,proj,Phij_coeff,Phij_NN,Shell2at,Sym)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in) :: Shell2at
  integer,intent(in) :: ntotcoeff
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4),proj(9,9,Shell2at%nshell)
  double precision,intent(in) :: Phij_coeff(ntotcoeff,1)
  double precision,intent(out) :: Phij_NN(3*InVar%natom,3*InVar%natom)

  integer :: iatcell,ishell,jshell,isym,iatom,jatom,eatom,fatom,ncoeff,ncoeff_prev
  integer :: nshell,ii,jj,kk,ll,this_shell,kappa,iatshell,iatref,trans
  double precision :: max_bound,min_bound,dist_shell,delta
  integer,allocatable :: tab_shell(:),counter(:,:,:)
  double precision,allocatable :: Phij_33(:,:),Phij_shell(:,:,:),correction(:,:,:),Phij_ref(:,:,:)

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '#### For each shell, list of coefficients (IFC), number of neighbours... ####'
  write(InVar%stdout,*) '#############################################################################'

  nshell=Shell2at%nshell
!==========================================================================================
!======== 1/ Build the Phij_NN ============================================================
!==========================================================================================
  ABI_MALLOC(Phij_ref,(3,3,nshell)); Phij_ref(:,:,:)=zero
  ABI_MALLOC(Phij_33,(3,3)) ; Phij_33(:,:)=0.d0
  do ishell=1,nshell
!   Build the 3x3 IFC per shell    
    ncoeff     =Shell2at%ncoeff(ishell)
    ncoeff_prev=Shell2at%ncoeff_prev(ishell)
    kappa=0
    do ii=1,3
      do jj=1,3
        kappa=kappa+1
        Phij_ref(ii,jj,ishell)=sum(proj(kappa,1:ncoeff,ishell)*Phij_coeff(ncoeff_prev+1:ncoeff_prev+ncoeff,1))
      end do  
    end do  
    do eatom=1,Invar%natom
!     Build the 3x3 IFC of an atom in this shell    
      if (Shell2at%neighbours(eatom,ishell)%n_interactions.eq.0) cycle
      do iatshell=1,Shell2at%neighbours(eatom,ishell)%n_interactions
        fatom=Shell2at%neighbours(eatom,ishell)%atomj_in_shell(iatshell)
        isym =Shell2at%neighbours(eatom,ishell)%sym_in_shell(iatshell)
        trans=Shell2at%neighbours(eatom,ishell)%transpose_in_shell(iatshell)
        if (fatom.lt.eatom) cycle
        call tdep_build_phij33(isym,Phij_ref(:,:,ishell),Phij_33,Sym,trans) 
!       Symetrization of the Phij_NN matrix
        Phij_NN((eatom-1)*3+1:(eatom-1)*3+3,3*(fatom-1)+1:3*(fatom-1)+3)=Phij_33(:,:)
        do ii=1,3
          do jj=1,3
            Phij_NN((fatom  -1)*3+ii,3*(eatom-1)+jj)=Phij_33(jj,ii)
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
          Phij_NN((eatom-1)*3+jj,(eatom-1)*3+kk)=Phij_NN((eatom-1)*3+jj,3*(eatom-1)+kk)&
&                                               -Phij_NN((eatom-1)*3+jj,3*(fatom-1)+kk)
        enddo
      enddo
    enddo
  enddo
  ABI_FREE(Phij_33)
  ABI_FREE(Phij_ref)

!==========================================================================================
!======== 2/ Symetrize the IFC for low symetry systems (if required) ======================
!==========================================================================================
! Detect the IFC terms which are non-symetric
! Modify the correspondingly non-symetric terms if the input variable Impose_symetry=2 or 3
  if (((InVar%Impose_symetry.eq.2).or.(InVar%Impose_symetry.eq.3))) then 
    ABI_MALLOC(correction,(InVar%natom,3,3)); correction(:,:,:)=zero
    ABI_MALLOC(counter,(InVar%natom,3,3));    counter(:,:,:)   =zero
    ABI_MALLOC(Phij_shell,(3,3,nshell)); Phij_shell(:,:,:)=zero
!   Compute the non-symetric contribution for each shell     
    do ishell=1,nshell
      do kk=1,3
        iatref=Shell2at%iatref(ishell)
        if (Shell2at%neighbours(iatref,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell2at%neighbours(iatref,ishell)%n_interactions
          jatom=Shell2at%neighbours(iatref,ishell)%atomj_in_shell(iatshell)
!FB20171213          if (iatom==jatom) cycle
          if (iatref==jatom) cycle 
          do ll=1,3
            Phij_shell(kk,ll,ishell)=Phij_shell(kk,ll,ishell)+Phij_NN((iatref-1)*3+kk,(jatom-1)*3+ll)
          end do !ll
        end do !iatshell
      end do !kk
!     Write the non-symetric contribution
      do kk=1,3
        do ll=1,3
          if (abs(Phij_shell(kk,ll,ishell)-Phij_shell(ll,kk,ishell)).gt.tol8) then
            write(InVar%stdout,'(a,i4,a,i4,1x,i4)') '  WARNING: For shell',ishell ,' and directions (kk,ll)=',kk,ll
            write(InVar%stdout,*) '          the shell gives a non symetric contribution to Dij:'
            write(InVar%stdout,'(a,1x,f9.6,1x,f9.6)') '           the Dij(qpt=0) and Dji(qpt=0)=',&
&             Phij_shell(kk,ll,ishell),Phij_shell(ll,kk,ishell)
            write(InVar%stdout,*) '          This could lead to a non-hermitian Dij matrix.'
            write(InVar%stdout,*) ' '
          end if  
        end do !ll
      end do !kk
!     Compute the total non-symetric contribution for each line of the IFC
!FB      correction(:,:,:)=zero
!FB      counter(:,:,:)   =zero
      do iatom=1,InVar%natom
        do kk=1,3
          if (Shell2at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
          do iatshell=1,Shell2at%neighbours(iatom,ishell)%n_interactions
            jatom=Shell2at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
            if (iatom==jatom) cycle
            do ll=1,3
              if (abs(Phij_shell(kk,ll,ishell)-Phij_shell(ll,kk,ishell)).gt.tol8) then
                correction(iatom,kk,ll)=correction(iatom,kk,ll)+Phij_NN((iatom-1)*3+kk,(jatom-1)*3+ll)
                counter(iatom,kk,ll)   =counter(iatom,kk,ll)+1
              end if  
            end do !ll 
          end do !jatom  
        end do !kk 
!       Verify that the number of non-symetric contributions is symetric      
        do kk=1,3
          do ll=1,3
            if (counter(iatom,kk,ll).ne.counter(iatom,ll,kk)) then
              ABI_BUG('The correction cannot be applied')
            end if
          end do !ll
        end do !kk 
      end do !iatom
    end do !ishell
    do ishell=1,nshell
      do iatom=1,InVar%natom
!       Apply the non-symetric contributions to all non-symetric terms in the IFC      
        do kk=1,3
          do iatshell=1,Shell2at%neighbours(iatom,ishell)%n_interactions
            jatom=Shell2at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
            do ll=1,3
              if (iatom==jatom) cycle
              if (abs(Phij_shell(kk,ll,ishell)-Phij_shell(ll,kk,ishell)).gt.tol8) then
                delta=(correction(iatom,kk,ll)-correction(iatom,ll,kk))/2.d0/counter(iatom,kk,ll)
                if (InVar%debug.and.iatom.le.InVar%natom_unitcell) then
                  write(InVar%stdout,'(a,f9.6,a,3(i4,1x))') 'Correction=',delta,' for iatom,jatom,shell=',iatom,jatom,ishell
                end if  
                Phij_NN((iatom-1)*3+kk,(iatom-1)*3+ll)=Phij_NN((iatom-1)*3+kk,(iatom-1)*3+ll)+delta
                Phij_NN((iatom-1)*3+kk,(jatom-1)*3+ll)=Phij_NN((iatom-1)*3+kk,(jatom-1)*3+ll)-delta
              end if !Phij_shell
            end do !ll
          end do !jatom
        end do !kk
      end do !iatom 
    end do !ishell  
    ABI_FREE(correction)
    ABI_FREE(counter)
    ABI_FREE(Phij_shell)
  end if !Impose_symetry 

!==========================================================================================
!======== 3/ Write the Phij_NN in output ==================================================
!==========================================================================================
! Remove the rounding errors before writing (for non regression testing purposes)
  do ii=1,3*InVar%natom
    do jj=1,3*InVar%natom
      if (abs(Phij_NN(ii,jj)).lt.tol8) Phij_NN(ii,jj)=zero
    end do
  end do  

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
            write(InVar%stdout,'(2x,3(f9.6,1x))') Phij_NN((iatcell-1)*3+ii,(jatom-1)*3+1),Phij_NN((iatcell-1)*3+ii,(jatom-1)*3+2),&
&             Phij_NN((iatcell-1)*3+ii,(jatom-1)*3+3)
          end do
          write(InVar%stdout,'(a,3(1x,f11.6))') '  The components of the vector are:', distance(iatcell,jatom,2:4)
          write(InVar%stdout,'(a,(1x,f9.6))') '  Trace=',Phij_NN((iatcell-1)*3+1,(jatom-1)*3+1)+Phij_NN((iatcell-1)*3+2,&
&           (jatom-1)*3+2)+Phij_NN((iatcell-1)*3+3,(jatom-1)*3+3)
          write(InVar%stdout,*) ' '
        end do
      end if
    end do !ishell 
  end do !iatcell  
  ABI_FREE(tab_shell)

! Write the Phij_unitcell.dat and Phij_NN.dat files
  if (InVar%debug) then
    write(InVar%stdout,'(a)') ' See the Phij*.dat file'
    open(unit=52,file=trim(InVar%output_prefix)//'Phij_unitcell.dat')
    open(unit=55,file=trim(InVar%output_prefix)//'Phij_NN.dat')
    do jatom=1,3*InVar%natom
      if (jatom.le.3*InVar%natom_unitcell) then
        write(52,'(10000(f10.6,1x))') Phij_NN(jatom,:)
      end if  
      write(55,'(10000(f10.6,1x))') Phij_NN(jatom,:)
    end do  
    close(52)
    close(55)
  end if  

end subroutine tdep_build_phijNN 

!=====================================================================================================
subroutine tdep_calc_dij(dij,eigenV,iqpt,InVar,Lattice,omega,Phij_NN,qpt_cart,Rlatt_cart)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  integer,intent(in) :: iqpt
  double precision,intent(in) :: Phij_NN(3*InVar%natom,3*InVar%natom)
  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)
  double precision,intent(in) :: qpt_cart(3)
  double precision,intent(out) :: omega (3*InVar%natom_unitcell)
  double complex  ,intent(out) :: dij   (3*InVar%natom_unitcell,3*InVar%natom_unitcell)
  double complex  ,intent(out) :: eigenV(3*InVar%natom_unitcell,3*InVar%natom_unitcell)

  integer :: LWORK,ii,jj,kk,iatom,jatom,iatcell,itypat,jtypat,iat_mod,INFO,itemp,imode,nmode
  double precision :: phase
  double complex :: ctemp,norm
  double precision, allocatable :: omega2(:)
  double precision, allocatable :: RWORK(:)
  double complex, allocatable :: WORKC(:)
! double complex, allocatable :: mass_mat(:,:)

! Calculation of the dynamical matrix (Dij)
  do iatcell=1,InVar%natom_unitcell
    itypat=InVar%typat_unitcell(iatcell)
    do jatom=1,InVar%natom
      iat_mod=mod(jatom+InVar%natom_unitcell-1,InVar%natom_unitcell)+1
      jtypat=InVar%typat_unitcell(iat_mod)
      phase=0.d0
      do kk=1,3
        phase=phase+2*pi*Rlatt_cart(kk,iatcell,jatom)*qpt_cart(kk)
      end do
      do ii=1+(iatcell-1)*3,3+(iatcell-1)*3
        do jj=1,3
          dij(ii,3*(iat_mod-1)+jj)=dij(ii,3*(iat_mod-1)+jj)+dcmplx(Phij_NN(ii,3*(jatom-1)+jj),0.d0)*exp(dcmplx(0.d0,phase))/&
&                              dcmplx(dsqrt(InVar%amu(itypat)*InVar%amu(jtypat))*amu_emass,0.d0)          
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
        if ((InVar%Impose_symetry.eq.0).or.(InVar%Impose_symetry.eq.2)) then
          write(InVar%stdout,*) 'STOP: The Dij matrix is not hermitian'
          write(InVar%stdout,'(a,1x,3(f10.6,1x))') 'For qpt=',qpt_cart(:)*Lattice%acell_unitcell(:)
          write(InVar%stdout,'(a,i4,a)') '  Dij(',iqpt,'real)='
          do iatcell=1,InVar%natom_unitcell
            write(InVar%stdout,'(100(f12.8,1x))') real(dij(1+(iatcell-1)*3,:))
            write(InVar%stdout,'(100(f12.8,1x))') real(dij(2+(iatcell-1)*3,:))
            write(InVar%stdout,'(100(f12.8,1x))') real(dij(3+(iatcell-1)*3,:))
          end do  
          write(InVar%stdout,'(a,i4,a)') '  Dij(',iqpt,'imag)='
          do iatcell=1,InVar%natom_unitcell
            write(InVar%stdout,'(100(f12.8,1x))') aimag(dij(1+(iatcell-1)*3,:))
            write(InVar%stdout,'(100(f12.8,1x))') aimag(dij(2+(iatcell-1)*3,:))
            write(InVar%stdout,'(100(f12.8,1x))') aimag(dij(3+(iatcell-1)*3,:))
          end do  
          ABI_ERROR('The Dij matrix is not hermitian')
        else if ((InVar%Impose_symetry.eq.1).or.(InVar%Impose_symetry.eq.3)) then
          ctemp=(dij(ii,jj)+conjg(dij(jj,ii)))/2.d0
          dij(ii,jj)=ctemp
          dij(jj,ii)=conjg(ctemp)
        end if  
      end if  
    end do
  end do
  if (itemp.ne.0.and.iqpt.eq.1) then
    write(InVar%stdout,*) 'WARNING: The Dij matrix is not hermitian'
    write(InVar%stdout,*) '  Probably: one shell may not have the whole number of atoms'
    write(InVar%stdout,*) '  The Dij matrix is symetrized'
  end if  

!FB! Remove the rounding errors before writing (for non regression testing purposes)
!FB  do ii=1,3*InVar%natom_unitcell
!FB    do jj=1,3*InVar%natom_unitcell
!FB      if (abs(real(dij(ii,jj))).lt.tol8) dij(ii,jj)=dcmplx(zero,imag(dij(ii,jj)))
!FB      if (abs(imag(dij(ii,jj))).lt.tol8) dij(ii,jj)=dcmplx(real(dij(ii,jj)),zero)
!FB    end do
!FB  end do  

! Diagonalization of dynamical matrix
  LWORK=2*3*InVar%natom_unitcell-1
  ABI_MALLOC(WORKC,(LWORK)); WORKC(:)=czero
  ABI_MALLOC(RWORK,(3*3*InVar%natom_unitcell-2)); RWORK(:)=zero
!FB  ABI_MALLOC(mass_mat,(3*InVar%natom_unitcell,3*InVar%natom_unitcell)); mass_mat(:,:)=dcmplx(0.d0,0.d0)
!FB  do iatom=1,InVar%natom_unitcell
!FB    itypat=InVar%typat_unitcell(iatom)
!FB    do ii=1,3
!FB      mass_mat(3*(iatom-1)+ii,3*(iatom-1)+ii)=dcmplx(InVar%amu(itypat)*amu_emass,0.d0)
!FB    end do
!FB  end do  
!FB  eigenV(:,:)=dij(:,:)
!FB  call ZHEGV(1,'V','U',3*InVar%natom_unitcell,eigenV(:,:),3*InVar%natom_unitcell,mass_mat(:,:),3*InVar%natom_unitcell,omega(:),WORKC,LWORK,RWORK,INFO)
!FB  ABI_FREE(mass_mat)
  eigenV(:,:)=dij(:,:)
  call ZHEEV('V','U',3*InVar%natom_unitcell,eigenV(:,:),3*InVar%natom_unitcell,omega(:),WORKC,LWORK,RWORK,INFO)

! Normalization of the eigenvectors
  nmode=3*InVar%natom_unitcell
  ABI_MALLOC(omega2,(nmode)); omega2(:)=zero
  do imode=1,nmode
    norm=zero
    do iatom=1,InVar%natom_unitcell
      do ii=1,3
        norm=norm+eigenV(3*(iatom-1)+ii,imode)*conjg(eigenV(3*(iatom-1)+ii,imode))
        do jatom=1,InVar%natom_unitcell
          do jj=1,3
            omega2(imode)=omega2(imode)+&
&             conjg(eigenV(3*(iatom-1)+ii,imode))*dij(3*(iatom-1)+ii,3*(jatom-1)+jj)*eigenV(3*(jatom-1)+jj,imode)
          end do
        end do  
      end do
    end do
    eigenV(:,imode)=eigenV(:,imode)/dsqrt(real(norm))
  end do
!FB  write(6,'(i5,x,100(f15.6,x))') iqpt,(omega2(ii)*Ha_cmm1**2,ii=1,nmode)
  ABI_FREE(omega2)

! Remove the negative frequencies  
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
subroutine tdep_write_dij(dij,eigenV,iqpt,InVar,Lattice,omega,qpt_cart)

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type),intent(in) :: Lattice
  integer,intent(in) :: iqpt
  double precision,intent(in) :: qpt_cart(3)
  double precision,intent(in) :: omega (3*InVar%natom_unitcell)
  double complex  ,intent(in) :: dij   (3*InVar%natom_unitcell,3*InVar%natom_unitcell)
  double complex  ,intent(in) :: eigenV(3*InVar%natom_unitcell,3*InVar%natom_unitcell)

  integer :: ii,iatcell

! Print the dynamical matrix (Dij)
  write(52,'(a,1x,3(f10.6,1x))') 'For qpt=',qpt_cart(:)*Lattice%acell_unitcell(:)
  write(52,'(a,i4,a)') '  Dij(',iqpt,'real)*1.d6='
  do iatcell=1,InVar%natom_unitcell
    write(52,'(100(f10.6,1x))') real(dij(1+(iatcell-1)*3,:))*1.d6
    write(52,'(100(f10.6,1x))') real(dij(2+(iatcell-1)*3,:))*1.d6
    write(52,'(100(f10.6,1x))') real(dij(3+(iatcell-1)*3,:))*1.d6
  end do  
  write(52,'(a,i4,a)') '  Dij(',iqpt,'imag)*1.d6='
  do iatcell=1,InVar%natom_unitcell
    write(52,'(100(f10.6,1x))') aimag(dij(1+(iatcell-1)*3,:))*1.d6
    write(52,'(100(f10.6,1x))') aimag(dij(2+(iatcell-1)*3,:))*1.d6
    write(52,'(100(f10.6,1x))') aimag(dij(3+(iatcell-1)*3,:))*1.d6
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

end subroutine tdep_write_dij

!=====================================================================================================
subroutine tdep_build_phij33(isym,Phij_ref,Phij_33,Sym,trans) 

  implicit none

  type(Symetries_Variables_type),intent(in) :: Sym
! type(Input_Variables_type),intent(in) :: InVar
  double precision, intent(in) :: Phij_ref(3,3)
  double precision, intent(out) :: Phij_33(3,3)
  integer,intent(in) :: isym,trans

  double precision :: Phij_tmp(3,3),tmp1(3,3)


! Transform in the new basis wrt S_ref
  call DGEMM('N','N',3,3,3,1.d0,Sym%S_ref(:,:,isym,1),3,Phij_ref,3,0.d0,Phij_tmp,3)
  call DGEMM('N','N',3,3,3,1.d0,Phij_tmp,3,Sym%S_inv(:,:,isym,1),3,0.d0,Phij_33,3)

  if ((trans.lt.1).or.(trans.gt.2)) then
    ABI_BUG('This value of the symmetry index is not permitted')
  end if
! Transpose the 3x3 matrix if required
  if (trans.eq.2) then
    tmp1(:,:)=Phij_33(:,:)
    Phij_33(1,2)=tmp1(2,1)
    Phij_33(1,3)=tmp1(3,1)
    Phij_33(2,3)=tmp1(3,2)
    Phij_33(2,1)=tmp1(1,2)
    Phij_33(3,1)=tmp1(1,3)
    Phij_33(3,2)=tmp1(2,3)
  end if  

end subroutine tdep_build_phij33

!=====================================================================================================
subroutine tdep_init_eigen2nd(Eigen2nd,natom_unitcell,nqpt)

  implicit none

  integer, intent(in) :: natom_unitcell,nqpt
  type(Eigen_Variables_type),intent(out) :: Eigen2nd

  ABI_MALLOC(Eigen2nd%eigenval,(3*natom_unitcell,nqpt));                  Eigen2nd%eigenval(:,:)=  zero
  ABI_MALLOC(Eigen2nd%eigenvec,(3*natom_unitcell,3*natom_unitcell,nqpt)); Eigen2nd%eigenvec(:,:,:)=zero

end subroutine tdep_init_eigen2nd

!=====================================================================================================
subroutine tdep_destroy_eigen2nd(Eigen2nd)

  implicit none

  type(Eigen_Variables_type),intent(inout) :: Eigen2nd

  ABI_FREE(Eigen2nd%eigenval)
  ABI_FREE(Eigen2nd%eigenvec)

end subroutine tdep_destroy_eigen2nd

!=====================================================================================================
subroutine tdep_write_yaml(Eigen2nd,Qpt,Prefix)

  implicit none

  type(Eigen_Variables_type),intent(in) :: Eigen2nd
  type(Qpoints_type),intent(in) :: Qpt
  character(len=*) :: Prefix
! type(Lattice_Variables_type),intent(in) :: Lattice

  integer :: ii,iqpt,imode,nmode,idir,iatom
  double precision :: distance
  
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
    write(52,'(a)')    '  band:'
    do imode=1,nmode
      write(52,'(a,i4)')    '  - #',imode
      write(52,'(a,f15.6)') '    frequency:',Eigen2nd%eigenval(imode,iqpt)*Ha_THz
      write(52,'(a)') '    eigenvector:'
      do iatom=1,nmode/3
        write(52,'(a,i4)') "    - # atom ", iatom
        do idir=1,3
          write(52,'(a,f18.9,a,f18.9,a)') "      - [",real(Eigen2nd%eigenvec((iatom-1)*3+idir,imode,iqpt)),','&
&           ,aimag(Eigen2nd%eigenvec((iatom-1)*3+idir,imode,iqpt)),']'
        end do
      end do
    end do !nmode  
    write(52,'(a)') ''
  end do

  close(52)


end subroutine tdep_write_yaml

!=====================================================================================================
end module m_tdep_phij
