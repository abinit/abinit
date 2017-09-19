#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_phij

  use defs_basis
  use m_readwrite,   only : Input_Variables_type
  use m_latt,        only : Lattice_Variables_type
  use m_sym,         only : Symetries_Variables_type
  use m_qpt,         only : Qpoints_type

  implicit none

  public :: calc_phij_fcoeff
  public :: calc_phij_nbcoeff

contains

!====================================================================================================
subroutine calc_phij_fcoeff(InVar,bond_ref,nshell,ntotcoeff,proj,shell,Sym,ucart,fcoeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_phij_fcoeff'
!End of the abilint section

  implicit none

  integer, intent(in) :: ntotcoeff,nshell
  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  integer, intent(in) :: shell(nshell,4)
  integer, intent(in) :: bond_ref(InVar%natom,InVar%natom,3)
  double precision, intent(in) :: ucart(3,InVar%natom,InVar%nstep_remain)
  double precision, intent(in) :: proj(9,9,nshell)
  double precision, intent(out) :: fcoeff(3*InVar%natom*InVar%nstep_remain,ntotcoeff)

  integer :: ishell,ncoeff,ncoeff_prev,istep,katom,latom
  integer :: icoeff,isym
  integer :: mu,nu,alpha,beta
  double precision :: terme
  double precision :: udiff(3),SSu(3,9)

! For each couple of atoms, transform the Phij (3x3) ifc matrix using the symetry operation (matR)
! Note: katom=1 is excluded in order to take into account the atomic sum rule (see below)
  do katom=1,InVar%natom
    do latom=1,InVar%natom
      if (katom==latom) cycle
      isym=abs(Sym%matR(katom,latom))
      ishell=bond_ref(katom,latom,3)
      if (ishell==0) cycle
      ncoeff     =shell(ishell,1)
      ncoeff_prev=shell(ishell,2)

      do istep=1,InVar%nstep_remain

!       In order to impose the acoustic sum rule we use (u(l)-u(k))
        udiff(1)=ucart(1,latom,istep)-ucart(1,katom,istep)
        udiff(2)=ucart(2,latom,istep)-ucart(2,katom,istep)
        udiff(3)=ucart(3,latom,istep)-ucart(3,katom,istep)
      
        SSu(:,:)=zero
            
   !    F_j^{\mu\nu}(t)=\sum_{\alpha\beta,i}\Phi_{ij}^{\alpha\beta}.S^{\alpha\mu}.S^{\beta\nu}.u_i^\nu(t)
        do mu=1,3
          do alpha=1,3
            do nu=1,3
              do beta=1,3
                if (Sym%matR(katom,latom).gt.0) then
                  SSu(mu,beta+(alpha-1)*3)=SSu(mu,beta+(alpha-1)*3)+Sym%matR_inv(alpha,mu,isym,1)*Sym%matR_inv(beta,nu,isym,1)*udiff(nu)
                else
                  SSu(mu,alpha+(beta-1)*3)=SSu(mu,alpha+(beta-1)*3)+Sym%matR_inv(alpha,mu,isym,1)*Sym%matR_inv(beta,nu,isym,1)*udiff(nu)
                end if
              end do
            end do
          end do
        end do
        do mu=1,3
          do icoeff=1,ncoeff
            terme=sum(SSu(mu,:)*proj(:,icoeff,ishell))
            fcoeff(mu+3*(katom-1)+3*InVar%natom*(istep-1),icoeff+ncoeff_prev)=fcoeff(mu+3*(katom-1)+3*InVar%natom*(istep-1),icoeff+ncoeff_prev)+terme
          end do    
        end do  
      
      end do !istep
    end do !latom
  end do !katom

end subroutine calc_phij_fcoeff

!=====================================================================================================
subroutine calc_phij_nbcoeff(distance,iatcell,InVar,ishell,jatom,ncoeff,nshell,proj,Sym)

  use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_phij_nbcoeff'
!End of the abilint section

  implicit none


  integer,intent(in) :: iatcell,ishell,jatom,nshell
  integer,intent(out) :: ncoeff
  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4)
  double precision,intent(out) :: proj(9,9,nshell)


  integer :: ii,jj,isym,LWORK,INFO,ieig,const_tot,jsym
  integer :: kk,ncount,icoeff,jatcell,mu,nu
  integer, allocatable :: IPIV(:)
  integer :: iconst(Sym%nsym)
  double precision :: prod_scal,drandom
  double precision :: eigvec(3,3),vect1(3),vect_trial(3),eigval(3)
  double precision :: WR(3),WI(3),VL(3,3),VR(3,3)
  double precision, allocatable :: WORK(:)
  double complex :: eigenvectors(3,3),eigenvalues(3),alphaij(Sym%nsym,9,9),pp(3,3),lambda
  double complex, allocatable :: tab_vec(:,:),temp(:,:)
  logical :: ok
  logical :: unchanged(Sym%nsym),KeptInvariant,Reversed

!=====================================================
!TODO: Check that transpose is equal to inverse!!!!!!!
!=====================================================

  if (jatom==iatcell) return

! Si on souhaite supprimer la prise en compte des symetries
!FB  if (.true.) then
!FB    do ii=1,9
!FB      proj(ii,ii,ishell)=1.d0
!FB    end do
!FB    ncoeff=9
!FB    return
!FB  end if  

  alphaij(:,:,:)=czero
  unchanged(:)=.false.
  const_tot=0

! Search of the inversion operation symmetry  

! On cherche R tel que: vect1=R.vect1
  do ii=1,3
    vect1(ii)=distance(iatcell,jatom,ii+1)
  end do
  do isym=1,Sym%nsym
    ok=.false.
!FB    if ((matR_ref(1,1,isym,2)==0).and.(isym.gt.2)) cycle
    iconst(isym)=0
    vect_trial(:)=zero
    do ii=1,3
      do jj=1,3
        vect_trial(ii)=vect_trial(ii)+Sym%matR_ref(ii,jj,isym,1)*vect1(jj)
      end do  
    end do
!   Search if the bond is kept invariant or reversed    
    KeptInvariant=.false.
    Reversed=.false.
    if (abs(vect_trial(1)-vect1(1)).lt.tol8.and.&
&       abs(vect_trial(2)-vect1(2)).lt.tol8.and.&
&       abs(vect_trial(3)-vect1(3)).lt.tol8) KeptInvariant=.true.
    if (abs(vect_trial(1)+vect1(1)).lt.tol8.and.&
&       abs(vect_trial(2)+vect1(2)).lt.tol8.and.&
&       abs(vect_trial(3)+vect1(3)).lt.tol8) Reversed=.true.

!   The bond (1,3) becomes (1,3) or (3,1) according to the kind of symetry
!   (Reversed or Keptinvariant)
    jatcell=mod(jatom-1,InVar%natom_unitcell)+1
    if (Keptinvariant.and.(Sym%indsym(4,isym,iatcell)==iatcell).and.&
&                         (Sym%indsym(4,isym,jatom  )==jatcell)) ok=.true.
    if (Reversed     .and.(Sym%indsym(4,isym,jatom  )==iatcell).and.&
&                         (Sym%indsym(4,isym,iatcell)==jatcell)) ok=.true.

!   We write the conditions of the calculation for this symetry
    if (.not.ok) cycle

!   If the bond is kept invariant or reversed    
    write(16,*) ' '
    write(16,*) 'For shell number=',ishell
    if (KeptInvariant) then
      write(16,*)'===========The bond is kept invariant for isym=',isym
    else        
      write(16,*)'===========The bond is reversed for isym=',isym
    end if        
    write(16,'(3(f16.12,x))') Sym%matR_ref(1,1,isym,1),Sym%matR_ref(1,2,isym,1),Sym%matR_ref(1,3,isym,1)
    write(16,'(3(f16.12,x))') Sym%matR_ref(2,1,isym,1),Sym%matR_ref(2,2,isym,1),Sym%matR_ref(2,3,isym,1)
    write(16,'(3(f16.12,x))') Sym%matR_ref(3,1,isym,1),Sym%matR_ref(3,2,isym,1),Sym%matR_ref(3,3,isym,1)

!   If the transformation matrix keep the bond invariant:
!       Phi_{\alpha\beta}=\sum_{\mu\nu} Phi_{\mu\nu}.S_{\mu\alpha}^T.S_{\nu\beta}^T
!       If lambda and p are the eigenvectors and eigenvalues of the S^T matrix, then:
!       \sum_{\alpha\beta} Phi_{\alpha\beta}.p_{\alpha}^l.p_{\beta}^k 
!     = \sum_{\mu\nu,\alpha\beta} Phi_{\mu\nu}.S_{\mu\alpha}^T.S_{\nu\beta}^T.p_{\alpha}^l.p_{\beta}^k
!     = lambda^l.lambda^k \sum_{\mu\nu} Phi_{\mu\nu}.p_{\mu}^l.p_{\nu}^k
!   So, if lambda_{\alpha}.lambda{\beta} = -1, we must have:
!      \sum_{\alpha\beta} Phi_{\alpha\beta}.p_{\alpha}^l.p_{\beta}^k = 0      
!
!   In the case of the reversed bond, one obtains the following constraint:
!      \sum_{\alpha\beta} Phi_{\alpha\beta}.(lambda^l.lambda^k.p_{\alpha}^l.p_{\beta}^k-p_{\beta}^l.p_{\alpha}^k) = 0 
!   which applies whether lambda^l.lambda^k = \pm 1

    do ii=1,3
      do jj=1,3
        eigvec(ii,jj)=Sym%matR_ref(jj,ii,isym,1)
      end do  
    end do  
    LWORK=4*3
    ABI_MALLOC(WORK,(LWORK)); WORK(:)=zero
!   La matrice de transformation eigvec est reelle et peut etre non-symetrique
    call dgeev( 'N', 'V', 3, eigvec, 3, WR, WI, VL, 3, VR, 3, WORK, LWORK, INFO)
    ABI_FREE(WORK)

!   On reconstruit les parties reelles et complexes des vecteurs/valeurs propres
    jj=0
    do ii=1,3
      eigenvalues(ii)=dcmplx(WR(ii),WI(ii))
      if (WI(ii).ne.zero.and.jj==0) then
        do kk=1,3
          eigenvectors(kk,ii)=dcmplx(VR(kk,ii),VR(kk,ii+1))
        end do
        jj=jj+1
      else if (WI(ii).ne.zero.and.jj==1) then
        do kk=1,3
          eigenvectors(kk,ii)=dcmplx(VR(kk,ii-1),-VR(kk,ii))
        end do
        jj=jj+1
      else
        do kk=1,3
          eigenvectors(kk,ii)=dcmplx(VR(kk,ii),zero)
        end do  
      end if
    end do  

!   On ecrit les valeurs propres et vecteurs propres
!   LES VALEURS PROPRES PEUVENT ETRE COMPLEXES
    do ii=1,3
      write(16,*)'  For eigenvalue number',ii
      write(16,*)'     The eigenvalue is:',eigenvalues(ii)
      write(16,*)'     The eigenvector is:'
      write(16,'(2(f16.12,x))') eigenvectors(1,ii)
      write(16,'(2(f16.12,x))') eigenvectors(2,ii)
      write(16,'(2(f16.12,x))') eigenvectors(3,ii)
      if ((aimag(eigenvalues(1)).ne.0).or.(aimag(eigenvalues(2)).ne.0).or.(aimag(eigenvalues(3)).ne.0)) then
        write(16,*) '  WARNING: THERE IS COMPLEX EIGENVALUES:'
      end if  
    end do  
    write(16,*) ' '

!   On calcule le produit des vecteurs propres: alphaij(i,j)=(p1)_i^T.(p2)_j      
!   Ce produit n'est realise que lorsque le produit des valeurs propres/=1 ou 0
!   On obtient n vecteurs a 9 coeff (definis dans R^9).
!   L'espace de recherche des solutions independantes (R^(9-n)) est orthogonal 
!   a ces n vecteurs.
    do ii=1,3
      do jj=1,3
        lambda=eigenvalues(ii)*eigenvalues(jj)
        if (((abs(real(lambda)-1.d0).lt.1.d-6).and.(abs(aimag(lambda)).lt.1.d-6).and.KeptInvariant).or.&
!FB&           ((aimag(eigenvalues(1)).ne.0).or.(aimag(eigenvalues(2)).ne.0).or.(aimag(eigenvalues(3)).ne.0)).or.&
&           ((abs(real(lambda)-1.d0).lt.1.d-6).and.(abs(aimag(lambda)).lt.1.d-6).and.Reversed.and.(ii==jj))) then
          cycle
        else  
          unchanged(isym)=.true.
          iconst(isym)=iconst(isym)+1
          const_tot=const_tot+1
          write(16,*)'  The product of eigenvalues',ii,jj
          write(16,*)'  is equal to ',lambda  
          do mu=1,3
            do nu=1,3
              pp(mu,nu)=eigenvectors(mu,ii)*eigenvectors(nu,jj)
            end do
          end do  
          do mu=1,3
            do nu=1,3
              if (KeptInvariant) then
                alphaij(isym,(mu-1)*3+nu,iconst(isym))=pp(mu,nu)
              else if (Reversed) then
                alphaij(isym,(mu-1)*3+nu,iconst(isym))=lambda*pp(mu,nu)-pp(nu,mu)
              else
                write(16,*)'  This symetry cannot be neither Keptinvariant nor Reversed'
              end if
            end do  
          end do  
          write(16,*)'  Real & imaginary parts of the eigenvectors product:'
          write(16,'(9(f16.12,x))')  real(alphaij(isym,:,iconst(isym)))
          write(16,'(9(f16.12,x))') aimag(alphaij(isym,:,iconst(isym)))
        end if !epsilon_i*epsilon_j=1 
      end do !jj
    end do !ii

!   WARNING: There are some minimum and maximum of constraints
    if ((iconst(isym).eq.1).and.Keptinvariant) then
      write(16,*)'  There is only one constraint'
      stop
    end if
    if ((iconst(isym).eq.2).and.Reversed) then
      write(16,*)'  There is two constraints'
      stop
    end if
    if (iconst(isym).gt.8) then
      write(16,*)'  There are more than 8 constraints'
      stop
    end if
  end do !isym

! In the case where the matrix has 9 elements  
  if (const_tot==0) then 
    write(16,*)'  WARNING: There is no symetry operation leaving the bond unchanged'
    write(16,*) isym,ishell,jatom,iatcell
    do ii=1,9
      proj(ii,ii,ishell)=1.d0
    end do
    ncoeff=9
    return
  end if

! When some constraints have been found
  ncount=const_tot
  write(16,*) ' '
  write(16,*) 'There is a total of ',ncount,' non-independant constraints for this shell'
  ii=0
  ABI_MALLOC(tab_vec,(9,ncount)); tab_vec(:,:)=czero
  do isym=1,Sym%nsym
    if (unchanged(isym)) then
      do jj=1,iconst(isym)
        ii=ii+1
        tab_vec(:,ii)=alphaij(isym,:,jj)
      end do  
    end if  
  end do
  if (ii.ne.ncount) then
    write(16,*) ii,' non equal to ',ncount
    stop
  end if  
  do ii=1,9
    do jj=1,ncount
      if (abs( real(tab_vec(ii,jj))).lt.1.d-8) tab_vec(ii,jj)=dcmplx(zero,aimag(tab_vec(ii,jj)))
      if (abs(aimag(tab_vec(ii,jj))).lt.1.d-8) tab_vec(ii,jj)=dcmplx( real(tab_vec(ii,jj)),zero)
    end do
  end do

! L'ensemble des vecteurs reduisants l'espace de R^9 a R^n ne forment pas une base
! independante. Il faut donc trouver les vecteurs independants.
! --> Orthogonalisation de Gram-Schmidt
  do kk=2,ncount
    do jj=1,kk-1
      prod_scal=sum( real(tab_vec(:,jj))* real(tab_vec(:,jj))+aimag(tab_vec(:,jj))*aimag(tab_vec(:,jj)))
      if (abs(prod_scal).gt.tol8) then
        tab_vec(:,kk)=tab_vec(:,kk)-sum(tab_vec(:,kk)*conjg(tab_vec(:,jj)))/dcmplx(prod_scal,zero)*tab_vec(:,jj)
        do ii=1,9
          if (abs( real(tab_vec(ii,kk))).lt.1.d-8) tab_vec(ii,kk)=dcmplx(zero,aimag(tab_vec(ii,kk)))
          if (abs(aimag(tab_vec(ii,kk))).lt.1.d-8) tab_vec(ii,kk)=dcmplx( real(tab_vec(ii,kk)),zero)
        end do
      end if  
    end do
  end do

! On stocke les vecteurs non-nuls
  ABI_MALLOC(temp,(9,ncount)); temp(:,:)=czero
  ii=0
  do kk=1,ncount
    prod_scal=sum( real(tab_vec(:,kk))* real(tab_vec(:,kk))+aimag(tab_vec(:,kk))*aimag(tab_vec(:,kk)))
    if (abs(prod_scal).gt.tol8) then
      ii=ii+1
      temp(:,ii)=tab_vec(:,kk)/dsqrt(prod_scal)
    end if  
  end do  
  ncount=ii
  ABI_FREE(tab_vec)

! On ecrit les vecteurs non-nuls
  write(16,*) ' '
  write(16,*) '  ========The final set of vectors is:'
  do kk=1,ncount
    write(16,'(9(f16.12,x))')  real(temp(:,kk))
    write(16,'(9(f16.12,x))') aimag(temp(:,kk))
  end do  
  write(16,*) '  =======Au total, il y a ',ncount,' vecteurs independants'
  if (ncount.gt.8) then
    write(16,*) '  ERROR : There is too many independant vectors'
    stop
  end if

! On cherche les (9-ncount) vecteurs orthogonaux aux vecteurs non-nuls
! --> Orthogonalisation de Gram-Schmidt
  ABI_MALLOC(tab_vec,(9,9)); tab_vec(:,:)=czero
  do kk=1,9
    if (kk.le.ncount) then
      tab_vec(:,kk)=temp(:,kk)
    else  
      do jj=1,9
        call random_number(drandom)
        tab_vec(jj,kk)=dcmplx(drandom,zero)
      end do  
      do jj=1,kk-1
        prod_scal=sum( real(tab_vec(:,jj))* real(tab_vec(:,jj))+aimag(tab_vec(:,jj))*aimag(tab_vec(:,jj)))
        if (abs(prod_scal).gt.tol8) then
          tab_vec(:,kk)=tab_vec(:,kk)-sum(tab_vec(:,kk)*conjg(tab_vec(:,jj)))/prod_scal*tab_vec(:,jj)
          do ii=1,9
            if (abs( real(tab_vec(ii,kk))).lt.1.d-8) tab_vec(ii,kk)=dcmplx(zero,aimag(tab_vec(ii,kk)))
            if (abs(aimag(tab_vec(ii,kk))).lt.1.d-8) tab_vec(ii,kk)=dcmplx( real(tab_vec(ii,kk)),zero)
          end do
        end if  
        prod_scal=sum( real(tab_vec(:,kk))* real(tab_vec(:,kk))+aimag(tab_vec(:,kk))*aimag(tab_vec(:,kk)))
        tab_vec(:,kk)=tab_vec(:,kk)/dsqrt(prod_scal)
      end do
    end if
  end do
  ABI_FREE(temp)

! On ecrit les vecteurs non-nuls
  write(16,*) ' '
  write(16,*) '  ========The orthogonal set of vectors is:'
  do kk=ncount+1,9
    write(16,'(9(f16.12,x))')  real(tab_vec(:,kk))
    write(16,'(9(f16.12,x))') aimag(tab_vec(:,kk))
    if ((abs(aimag(tab_vec(1,kk))).gt.1.d-6).or.&
&       (abs(aimag(tab_vec(1,kk))).gt.1.d-6).or.&
&       (abs(aimag(tab_vec(1,kk))).gt.1.d-6)) then
      write(16,*) ' ERROR : the constraint has an imaginary part'
      stop
    end if
  end do  
  ncoeff=9-ncount
  write(16,*) '  =======Au total, il y a ',ncoeff,' coefficients'

! On copie tab_vec dans proj
  do icoeff=1,ncoeff
    proj(:,icoeff,ishell)=tab_vec(:,ncount+icoeff)
  end do
  ABI_FREE(tab_vec)

end subroutine calc_phij_nbcoeff

!=====================================================================================================
subroutine build_phij(distance,InVar,bond_ref,nshell,ntotcoeff,proj,Phij_coeff,Phij_NN,shell,Sym)

  use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'build_phij'
!End of the abilint section

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  integer,intent(in) :: nshell,ntotcoeff
  integer,intent(in) :: bond_ref(InVar%natom,InVar%natom,3)
  integer,intent(in) :: shell(nshell,4)
  double precision,intent(in) :: distance(InVar%natom,InVar%natom,4),proj(9,9,nshell)
  double precision,intent(in) :: Phij_coeff(ntotcoeff,1)
  double precision,intent(out) :: Phij_NN(3*InVar%natom,3*InVar%natom)

  integer :: iatcell,ishell,jshell,isym,iatom,jatom,katom,latom,ncoeff,ncoeff_prev
  integer :: ii,jj,kk,ll,this_shell,kshell,jatcell
  double precision :: sum1,sum2,max_bound,min_bound,dist_shell,delta
  integer,allocatable :: atoms_in_shell(:),tab_shell(:),counter(:,:,:),nb_atoms(:)
  double precision,allocatable :: Phij_33(:,:),Phij_tmp(:,:),Phij_shell(:,:,:),correction(:,:,:)
  double precision :: tmp1(3,3),tmp3(3)
  logical :: ok

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '#### For each shell, list of coefficients (IFC), number of neighbours... ####'
  write(InVar%stdout,*) '#############################################################################'

  ABI_MALLOC(atoms_in_shell,(nshell)) ; atoms_in_shell(:)=0
  ABI_MALLOC(tab_shell,(nshell)); tab_shell(:)=0
!==========================================================================================
!======== 1/ Build the Phij_NN ============================================================
!==========================================================================================
  ABI_MALLOC(Phij_33,(3,3)) ; Phij_33(:,:)=0.d0
  ABI_MALLOC(Phij_tmp,(3,3)) ; Phij_tmp(:,:)=0.d0
  do katom=1,InVar%natom
    do latom=katom+1,InVar%natom
      if ((bond_ref(katom,latom,1).ne.0).and.(bond_ref(katom,latom,2).ne.0)) then 
        Phij_33(:,:)=0.d0
        iatcell=bond_ref(katom,latom,1)
        jatom  =bond_ref(katom,latom,2)
        ishell =bond_ref(katom,latom,3)
        ncoeff     =shell(ishell,1)
        ncoeff_prev=shell(ishell,2)

        do ii=1,9
          jj=mod(ii-1,3)+1
          kk=(ii-jj)/3+1
          Phij_33(kk,jj)=sum(proj(ii,1:ncoeff,ishell)*Phij_coeff(ncoeff_prev+1:ncoeff_prev+ncoeff,1))
        end do  

!       Compute the number of atoms in a shell
        if ((bond_ref(katom,latom,1).eq.katom)) then
          atoms_in_shell(ishell)=atoms_in_shell(ishell)+1
        end if  
        isym=abs(Sym%matR(katom,latom))

!       Transform in the new basis wrt matR_ref
        call DGEMM('N','N',3,3,3,1.d0,Sym%matR_ref(:,:,isym,1),3,Phij_33,3,0.d0,Phij_tmp,3)
        call DGEMM('N','N',3,3,3,1.d0,Phij_tmp,3,Sym%matR_inv(:,:,isym,1),3,0.d0,Phij_33,3)

!       Transpose the 3x3 matrix if required
        if (Sym%matR(katom,latom).lt.0) then
          tmp1(:,:)=Phij_33(:,:)
          Phij_33(1,2)=tmp1(2,1)
          Phij_33(1,3)=tmp1(3,1)
          Phij_33(2,3)=tmp1(3,2)
          Phij_33(2,1)=tmp1(1,2)
          Phij_33(3,1)=tmp1(1,3)
          Phij_33(3,2)=tmp1(2,3)
        end if

!       Symetrization of the Phij_NN matrix
        Phij_NN((katom-1)*3+1:(katom-1)*3+3,3*(latom-1)+1:3*(latom-1)+3)=Phij_33(:,:)
        do ii=1,3
          do jj=1,3
            Phij_NN((latom  -1)*3+ii,3*(katom-1)+jj)=Phij_33(jj,ii)
          end do        
        end do  
      end if  
    end do !latom
!   Regle de somme acoustique
    do jj=1,3
      do kk=1,3
        do latom=1,InVar%natom
          if (latom==katom) cycle
          Phij_NN((katom-1)*3+jj,(katom-1)*3+kk)=Phij_NN((katom-1)*3+jj,3*(katom-1)+kk)&
&                                               -Phij_NN((katom-1)*3+jj,3*(latom-1)+kk)
        enddo
      enddo
    enddo
  end do !katom
  ABI_FREE(Phij_33)
  ABI_FREE(Phij_tmp)

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
    do katom=1,InVar%natom
      Phij_shell(:,:,:)=0.d0
      do kk=1,3
        do latom=1,InVar%natom
          do ll=1,3
            if (katom==latom) cycle
            this_shell=bond_ref(katom,latom,3)
            if (this_shell.eq.0) cycle
            Phij_shell(kk,ll,this_shell)=Phij_shell(kk,ll,this_shell)+Phij_NN((katom-1)*3+kk,(latom-1)*3+ll)
          end do !ll
        end do !latom
      end do !kk
    end do !katom
!   Write the non-symetric contribution
    do ishell=1,nshell
      do kk=1,3
        do ll=1,3
          if (abs(Phij_shell(kk,ll,ishell)-Phij_shell(ll,kk,ishell)).gt.tol8) then
            write(InVar%stdout,'(a,i4,a,i4,x,i4)') '  WARNING: For shell',ishell ,' and directions (kk,ll)=',kk,ll
            write(InVar%stdout,*) '          the shell gives a non symetric contribution to Dij:'
            write(InVar%stdout,'(a,x,f9.6,x,f9.6)') '           the Dij(qpt=0) and Dji(qpt=0)=',Phij_shell(kk,ll,ishell),Phij_shell(ll,kk,ishell)
            write(InVar%stdout,*) '          This could lead to a non-hermitian Dij matrix.'
            write(InVar%stdout,*) ' '
          end if  
        end do !ll
      end do !kk
    end do !ishell  
!   Compute the total non-symetric contribution for each line of the IFC
    do katom=1,InVar%natom
      correction(:,:,:)=zero
      counter(:,:,:)   =zero
      do kk=1,3
        do latom=1,InVar%natom
          do ll=1,3
            if (katom==latom) cycle
            this_shell=bond_ref(katom,latom,3)
            if (this_shell.eq.0) cycle
            if (abs(Phij_shell(kk,ll,this_shell)-Phij_shell(ll,kk,this_shell)).gt.tol8) then
              correction(katom,kk,ll)=correction(katom,kk,ll)+Phij_NN((katom-1)*3+kk,(latom-1)*3+ll)
              counter(katom,kk,ll)   =counter(katom,kk,ll)+1
            end if  
          end do !ll 
        end do !latom 
      end do !kk 
!     Verify that the number of non-symetric contributions is symetric      
      do kk=1,3
        do ll=1,3
          if (counter(katom,kk,ll).ne.counter(katom,ll,kk)) then
            write(InVar%stdout,*) ' BUG: The correction cannot be applied'
            stop
          end if
        end do !ll
      end do !kk 
!     Apply the non-symetric contributions to all non-symetric terms in the IFC      
      do kk=1,3
        do latom=1,InVar%natom
          do ll=1,3
            if (katom==latom) cycle
            this_shell=bond_ref(katom,latom,3)
            if (this_shell.eq.0) cycle
            if (abs(Phij_shell(kk,ll,this_shell)-Phij_shell(ll,kk,this_shell)).gt.tol8) then
              delta=(correction(katom,kk,ll)-correction(katom,ll,kk))/2.d0/counter(katom,kk,ll)
              if (InVar%debug.and.katom.le.InVar%natom_unitcell) then
                write(InVar%stdout,'(a,f9.6,a,3(i4,x))') 'Correction=',delta,' for iatom,jatom,shell=',katom,latom,this_shell
              end if  
              Phij_NN((katom-1)*3+kk,(katom-1)*3+ll)=Phij_NN((katom-1)*3+kk,(katom-1)*3+ll)+delta
              Phij_NN((katom-1)*3+kk,(latom-1)*3+ll)=Phij_NN((katom-1)*3+kk,(latom-1)*3+ll)-delta
            end if !Phij_shell
          end do !ll
        end do !latom
      end do !kk
    end do !katom 
    ABI_FREE(correction)
    ABI_FREE(counter)
    ABI_FREE(Phij_shell)
  end if !Impose_symetry 

!==========================================================================================
!======== 3/ Write the Phij_NN in output ==================================================
!==========================================================================================
  ABI_MALLOC(nb_atoms,(InVar%natom_unitcell)); nb_atoms(:)=zero
! Remove the rounding errors before writing (for non regression testing purposes)
  do ii=1,3*InVar%natom
    do jj=1,3*InVar%natom
      if (abs(Phij_NN(ii,jj)).lt.tol8) Phij_NN(ii,jj)=zero
    end do
  end do  

! Write the IFCs in the data.out file (with others specifications: 
! number of atoms in a shell, distance, Trace...)
  do iatcell=1,InVar%natom_unitcell
    tab_shell(:)=0
    write(InVar%stdout,'(a,i4)') ' ############# List of IFC for the reference atom=',iatcell
!   Sort the IFC with distance in increasing order
    min_bound=-1.d0
    do ishell=1,nshell
      do jshell=1,nshell
        if ((distance(shell(jshell,3),shell(jshell,4),1).ge.min_bound).and.(tab_shell(jshell).eq.0)) then
          max_bound=distance(shell(jshell,3),shell(jshell,4),1)
          this_shell=jshell
        end if  
      end do

      do jshell=1,nshell
        if ((distance(shell(jshell,3),shell(jshell,4),1).lt.max_bound).and.&
&           (distance(shell(jshell,3),shell(jshell,4),1).ge.min_bound).and.&
&            (tab_shell(jshell).eq.0)) then
          max_bound=distance(shell(jshell,3),shell(jshell,4),1)
          this_shell=jshell
        end if
      end do
      tab_shell(this_shell)=1
      min_bound=max_bound
      dist_shell=distance(shell(this_shell,3),shell(this_shell,4),1)
      
!     Write the IFC properly  
      write(InVar%stdout,'(a,i4,a,i4,a,f9.6)') ' ======== NEW SHELL (ishell=',this_shell,&
&           '): There are',atoms_in_shell(this_shell),' atoms on this shell at distance=',dist_shell
      ok=.false.
      do jatom=1,InVar%natom
        if (bond_ref(iatcell,jatom,2).eq.0.d0) cycle
        if (bond_ref(iatcell,jatom,1).eq.0.d0) cycle
        if (bond_ref(iatcell,jatom,3).eq.this_shell) then
          jatcell=mod(jatom-1,InVar%natom_unitcell)+1
          nb_atoms(jatcell)=nb_atoms(jatcell)+1
          ok=.true.
          write(InVar%stdout,'(a,i4,a,i4)') '  For jatom=',jatom,' ,with label=',jatcell
          do ii=1,3
            write(InVar%stdout,'(2x,3(f9.6,x))') Phij_NN((iatcell-1)*3+ii,(jatom-1)*3+1),Phij_NN((iatcell-1)*3+ii,(jatom-1)*3+2),Phij_NN((iatcell-1)*3+ii,(jatom-1)*3+3)
          end do
          write(InVar%stdout,'(a,3(x,f9.6))') '  The components of the vector are:', distance(iatcell,jatom,2:4)
          write(InVar%stdout,'(a,3(x,f9.6))') '  Trace=',Phij_NN((iatcell-1)*3+1,(jatom-1)*3+1)+Phij_NN((iatcell-1)*3+2,(jatom-1)*3+2)+Phij_NN((iatcell-1)*3+3,(jatom-1)*3+3)
          write(InVar%stdout,*) ' '
        end if
      end do
      if (ok) then
        do jatcell=1,InVar%natom_unitcell
          write(InVar%stdout,'(a,i4,a,i4,a)') & 
&               '  At this shell, there are',nb_atoms(jatcell),' atoms labeled',jatcell,' in the sphere.'
        end do
      end if
      if (.not.ok) write(InVar%stdout,*) '  There is no IFC (no neighbour) for this shell'

    end do !ishell 
  end do !iatcell  

  ABI_FREE(atoms_in_shell)
  ABI_FREE(tab_shell)
  ABI_FREE(nb_atoms)

! Write the Phij_unitcell.dat and Phij_NN.dat files
  if (InVar%debug) then
    write(InVar%stdout,'(a)') ' See the Phij*.dat file'
    open(unit=52,file='Phij_unitcell.dat')
    open(unit=55,file='Phij_NN.dat')
    do jatom=1,3*InVar%natom
      if (jatom.le.3*InVar%natom_unitcell) then
        write(52,'(10000(f10.6,x))') Phij_NN(jatom,:)
      end if  
      write(55,'(10000(f10.6,x))') Phij_NN(jatom,:)
    end do  
    close(52)
    close(55)
  end if  

end subroutine build_phij 

!=====================================================================================================
subroutine calc_dij(InVar,Lattice,Phij_NN,Qpt,Rlatt_cart)

  use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'calc_dij'
!End of the abilint section

  implicit none

  type(Input_Variables_type),intent(in) :: InVar
  type(Lattice_Variables_type) :: Lattice
  type(Qpoints_type) :: Qpt
  double precision,intent(in) :: Phij_NN(3*InVar%natom,3*InVar%natom)
  double precision,intent(in) :: Rlatt_cart(3,InVar%natom_unitcell,InVar%natom)

  integer :: LWORK,iqpt,ii,jj,kk,iatom,jatom,iatcell,itypat,iat_mod,INFO,itemp
  double precision :: phase
  double complex :: ctemp
  double precision, allocatable :: RWORK(:),omega(:,:)
  double complex, allocatable :: dij(:,:,:),WORKC(:),mass_mat(:,:)
  logical :: ok


  write(InVar%stdout,*)' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '######################## Dynamical matrix ###################################'
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,'(a)') ' See the dij.dat file'
! Calcul du Dij
  ABI_MALLOC(dij,(3*InVar%natom_unitcell,3*InVar%natom_unitcell,Qpt%nqpt)) ; dij(:,:,:)=zero
  do iqpt=1,Qpt%nqpt 
    do iatcell=1,InVar%natom_unitcell
!     Calcul de la phase et de la matrice dynamique      
      do jatom=1,InVar%natom
        iat_mod=mod(jatom+InVar%natom_unitcell-1,InVar%natom_unitcell)
        phase=0.d0
        do kk=1,3
          phase=phase+2*pi*Rlatt_cart(kk,iatcell,jatom)*Qpt%qpt_cart(kk,iqpt)
        end do
        do ii=1+(iatcell-1)*3,3+(iatcell-1)*3
          do jj=1,3
            dij(ii,jj+iat_mod*3,iqpt)=dij(ii,jj+iat_mod*3,iqpt)+dcmplx(Phij_NN(ii,3*(jatom-1)+jj),0.d0)*exp(dcmplx(0.d0,phase))
          end do !jj
        end do !ii 
        if (InVar%debug.and.iqpt==1) then
          write(InVar%stdout,'(a,i4,a,i4)') 'For iatcell=',iatcell,' and jatom=',jatom
          write(InVar%stdout,'(a11,i4,x,a)') '  Dij(iqpt=',iqpt,', real)='
          write(InVar%stdout,'(100(f10.6,x))') real(dij(1+(iatcell-1)*3,:,iqpt))
          write(InVar%stdout,'(100(f10.6,x))') real(dij(2+(iatcell-1)*3,:,iqpt))
          write(InVar%stdout,'(100(f10.6,x))') real(dij(3+(iatcell-1)*3,:,iqpt))
          write(InVar%stdout,'(a11,i4,x,a)') '  Dij(iqpt=',iqpt,', imag)='
          write(InVar%stdout,'(100(f10.6,x))') imag(dij(1+(iatcell-1)*3,:,iqpt))
          write(InVar%stdout,'(100(f10.6,x))') imag(dij(2+(iatcell-1)*3,:,iqpt))
          write(InVar%stdout,'(100(f10.6,x))') imag(dij(3+(iatcell-1)*3,:,iqpt))
          write(InVar%stdout,'(a,i4,a)') '  '
        end if
      end do !jatom 
    end do !iatcell
  end do !iqpt

! Diagonalization of dynamical matrix
  LWORK=2*3*InVar%natom_unitcell-1
  ABI_MALLOC(omega,(3*InVar%natom_unitcell,Qpt%nqpt)); omega(:,:)=zero
  ABI_MALLOC(WORKC,(LWORK))
  ABI_MALLOC(RWORK,(3*3*InVar%natom_unitcell-2)); WORKC(:)=czero; RWORK(:)=zero
  open(unit=53,file='omega.dat')
  open(unit=52,file='dij.dat')
  open(unit=51,file='eigenvectors.dat')
  itemp=0
  do iqpt=1,Qpt%nqpt
!   The Dij has to be an hermitian matrix
    ok=.false.
    do ii=1,3*InVar%natom_unitcell
      do jj=ii,3*InVar%natom_unitcell
        if ((abs(real(dij(ii,jj,iqpt))-real(dij(jj,ii,iqpt))).gt.tol8).or.(abs(imag(dij(ii,jj,iqpt))+imag(dij(jj,ii,iqpt))).gt.tol8)) then
          ok=.true.
          itemp=itemp+1
          if ((InVar%Impose_symetry.eq.0).or.(InVar%Impose_symetry.eq.2)) then
            write(InVar%stdout,*) 'STOP: The Dij matrix is not hermitian'
            write(InVar%stdout,'(a,x,3(f10.6,x))') 'For qpt=',Qpt%qpt_cart(:,iqpt)*Lattice%acell_unitcell(:)
            write(InVar%stdout,'(a,i4,a)') '  Dij(',iqpt,'real)='
            do iatcell=1,InVar%natom_unitcell
              write(InVar%stdout,'(100(f10.6,x))') real(dij(1+(iatcell-1)*3,:,iqpt))
              write(InVar%stdout,'(100(f10.6,x))') real(dij(2+(iatcell-1)*3,:,iqpt))
              write(InVar%stdout,'(100(f10.6,x))') real(dij(3+(iatcell-1)*3,:,iqpt))
            end do  
            write(InVar%stdout,'(a,i4,a)') '  Dij(',iqpt,'imag)='
            do iatcell=1,InVar%natom_unitcell
              write(InVar%stdout,'(100(f10.6,x))') imag(dij(1+(iatcell-1)*3,:,iqpt))
              write(InVar%stdout,'(100(f10.6,x))') imag(dij(2+(iatcell-1)*3,:,iqpt))
              write(InVar%stdout,'(100(f10.6,x))') imag(dij(3+(iatcell-1)*3,:,iqpt))
            end do  
            stop
          else if ((InVar%Impose_symetry.eq.1).or.(InVar%Impose_symetry.eq.3)) then
            ctemp=(dij(ii,jj,iqpt)+conjg(dij(jj,ii,iqpt)))/2.d0
            dij(ii,jj,iqpt)=ctemp
            dij(jj,ii,iqpt)=conjg(ctemp)
          end if  
        end if  
      end do
    end do

!   Remove the rounding errors before writing (for non regression testing purposes)
    do ii=1,3*InVar%natom_unitcell
      do jj=1,3*InVar%natom_unitcell
        if (abs(real(dij(ii,jj,iqpt))).lt.tol8) dij(ii,jj,iqpt)=dcmplx(zero,imag(dij(ii,jj,iqpt)))
        if (abs(imag(dij(ii,jj,iqpt))).lt.tol8) dij(ii,jj,iqpt)=dcmplx(real(dij(ii,jj,iqpt)),zero)
      end do
    end do  
!   Print the Dij matrix 
    write(52,'(a,x,3(f10.6,x))') 'For qpt=',Qpt%qpt_cart(:,iqpt)*Lattice%acell_unitcell(:)
    write(52,'(a,i4,a)') '  Dij(',iqpt,'real)='
    do iatcell=1,InVar%natom_unitcell
      write(52,'(100(f10.6,x))') real(dij(1+(iatcell-1)*3,:,iqpt))
      write(52,'(100(f10.6,x))') real(dij(2+(iatcell-1)*3,:,iqpt))
      write(52,'(100(f10.6,x))') real(dij(3+(iatcell-1)*3,:,iqpt))
    end do  
    write(52,'(a,i4,a)') '  Dij(',iqpt,'imag)='
    do iatcell=1,InVar%natom_unitcell
      write(52,'(100(f10.6,x))') imag(dij(1+(iatcell-1)*3,:,iqpt))
      write(52,'(100(f10.6,x))') imag(dij(2+(iatcell-1)*3,:,iqpt))
      write(52,'(100(f10.6,x))') imag(dij(3+(iatcell-1)*3,:,iqpt))
    end do  
    write(52,*)' '

!   Diagonalisation of Dij
    ABI_MALLOC(mass_mat,(3*InVar%natom_unitcell,3*InVar%natom_unitcell)); mass_mat(:,:)=dcmplx(0.d0,0.d0)
    do iatom=1,InVar%natom_unitcell
      itypat=InVar%typat_unitcell(iatom)
      do ii=1,3
        mass_mat(3*(iatom-1)+ii,3*(iatom-1)+ii)=dcmplx(InVar%amu(itypat)*1.660538782e-27/9.10938215e-31,0.d0)
      end do
    end do  
    call ZHEGV(1,'V','U',3*InVar%natom_unitcell,dij(:,:,iqpt),3*InVar%natom_unitcell,mass_mat(:,:),3*InVar%natom_unitcell,omega(:,iqpt),WORKC,LWORK,RWORK,INFO)
    ABI_FREE(mass_mat)
    do ii=1,InVar%natom_unitcell
      do jj=1,3
        if (omega((ii-1)*3+jj,iqpt).lt.0.d0) then
          omega((ii-1)*3+jj,iqpt)=-dsqrt(-omega((ii-1)*3+jj,iqpt))*219474.64
        else
          omega((ii-1)*3+jj,iqpt)=dsqrt(omega((ii-1)*3+jj,iqpt))*219474.64
        end if
      end do  
    end do
    if (iqpt.le.Qpt%nqpt) write(53,'(i5,x,100(f15.6,x))') iqpt,(omega(ii,iqpt),ii=1,3*InVar%natom_unitcell)
    write(51,*) 'For iqpt=',iqpt
    do ii=1,3*InVar%natom_unitcell
      write(51,*) 'Mode number',ii,' energy',omega(ii,iqpt)
      write(51,*) '  Real:'
      write(51,*) real(dij(:,ii,iqpt))
      write(51,*) '  Imag:'
      write(51,*) imag(dij(:,ii,iqpt))
    end do  
    write(51,*) ' '
  end do  
  if (itemp.ne.0) then
    write(InVar%stdout,*) 'WARNING: The Dij matrix is not hermitian'
    write(InVar%stdout,*) '  Probably: one shell may not have the whole number of atoms'
    write(InVar%stdout,*) '  The Dij matrix is symetrized'
  end if  
  close(53)
  close(52)
  close(51)
  ABI_FREE(WORKC)
  ABI_FREE(RWORK)
  ABI_FREE(dij)
  ABI_FREE(omega)
end subroutine calc_dij

!=====================================================================================================
end module m_phij
