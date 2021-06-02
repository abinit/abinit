
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_constraints

  use defs_basis
  use m_errors
  use m_profiling_abi
  use m_tdep_readwrite,   only : Input_Variables_type
  use m_tdep_sym,         only : Symetries_Variables_type
  use m_tdep_utils,       only : Coeff_Moore_type
  use m_tdep_shell,       only : Shell_Variables_type
  use m_tdep_psij,        only : Shell_Variables_type, tdep_build_psij333

  implicit none

  public :: tdep_calc_orthonorm
  public :: tdep_calc_constraints
  public :: tdep_check_constraints

contains 

!====================================================================================================
subroutine tdep_calc_orthonorm(dim1,dim2,nindep,vectin,vectout)

  integer, intent(in) :: dim1,dim2
  integer, intent(out) :: nindep
  double precision, intent(inout) :: vectin(dim1,dim2)
  double precision, intent(out) :: vectout(dim1,dim2)

  integer :: ii,jj,kk
  double precision :: prod_scal

  ii=0
  do kk=1,dim2
    if (sum(abs(vectin(:,kk))).gt.tol8) then
      ii=ii+1
      vectout(:,ii)=vectin(:,kk)
    end if 
  end do
  nindep=ii
  return

! Gram-Schmidt orthogonalization
  do kk=2,dim2
    do jj=1,kk-1
      prod_scal=sum(vectin(:,jj)*vectin(:,jj))
      if (abs(prod_scal).gt.tol8) then
        vectin(:,kk)=vectin(:,kk)-sum(vectin(:,kk)*vectin(:,jj))/prod_scal*vectin(:,jj)
      end if  
    end do
  end do

! Store the non-zero vectors and normalize
  ii=0
  do kk=1,dim2
    prod_scal=sum(vectin(:,kk)*vectin(:,kk))
    if (abs(prod_scal).gt.tol8) then
      ii=ii+1
      vectout(:,ii)=vectin(:,kk)/dsqrt(prod_scal)
    end if  
  end do  
  nindep=ii

end subroutine tdep_calc_orthonorm
!====================================================================================================
subroutine tdep_calc_constraints(CoeffMoore,distance,InVar,nshell1at,nshell2at,nshell3at,Sym,&
&                                proj1st,Shell1at,&!optional
&                                proj2nd,Shell2at,&!optional 
&                                proj3rd,Shell3at) !optional

  type(Coeff_Moore_type), intent(inout) :: CoeffMoore
  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(Shell_Variables_type),intent(in), optional :: Shell1at
  type(Shell_Variables_type),intent(in), optional :: Shell2at
  type(Shell_Variables_type),intent(in), optional :: Shell3at
  integer, intent(in) :: nshell1at,nshell2at,nshell3at
  double precision, intent(in), optional :: proj1st(3 ,3 ,nshell1at)
  double precision, intent(in), optional :: proj2nd(9 ,9 ,nshell2at)
  double precision, intent(in), optional :: proj3rd(27,27,nshell3at)
  double precision, intent(in) :: distance(InVar%natom,InVar%natom,4)

  integer :: ishell,ncoeff,ncoeff_prev,iatom,jatom,iatshell !katom,
  integer :: icoeff,iconst,nconst_loc,iconst_loc,iconst_new,isym,ntotcoeff,iat_mod
  integer :: mu,nu,xi,alpha,beta,gama,lambda,trans,natom_unitcell,natom !,ii,trans1,trans2,trans3,trans4,trans5
  double precision :: terme,temp !,terme1,terme2,terme3,terme4,terme5
  double precision :: Levi_Civita(3,3,3)
  double precision, allocatable :: SS_ref(:,:,:,:,:),dlevi(:,:,:,:)
  double precision, allocatable :: vectin(:,:),vectout(:,:)
  double precision, allocatable :: SSS_ref(:,:,:,:,:,:)
  double precision, allocatable :: const_rot1st(:,:)
  double precision, allocatable :: const_rot2nd(:,:,:,:)
  double precision, allocatable :: const_dynmat(:,:,:,:,:)
  double precision, allocatable :: const_huang(:,:,:,:,:)
  !double precision, allocatable :: const_asr3rd(:,:,:,:,:,:,:)
!FB  double precision, allocatable :: const_rot3rd(:,:,:,:,:,:)
  !double precision, allocatable :: const_rot3rd(:,:,:,:,:,:,:,:)

  ABI_UNUSED(shell3at%nshell)

  natom_unitcell=InVar%natom_unitcell
  natom         =InVar%natom

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '########################## Add the constraints ##############################'

! For each couple of atoms, transform the Phij (3x3) ifc matrix using the symetry operation (S)
  if (present(proj2nd)) then
    ABI_MALLOC(SS_ref,(3,9,3,Sym%nsym,2)); SS_ref(:,:,:,:,:)=zero
    do isym=1,Sym%nsym
      do alpha=1,3
        do mu=1,3
          do beta=1,3
            do nu=1,3
              temp=Sym%S_ref(alpha,mu,isym,1)*Sym%S_ref(beta,nu,isym,1)
              SS_ref(alpha,nu+(mu-1)*3,beta,isym,1)=temp
              SS_ref(alpha,mu+(nu-1)*3,beta,isym,2)=temp
            end do
          end do
        end do
      end do
    end do  
  end if 
! For each couple of atoms, transform the Psij (3x3x3) ifc matrix using the symetry operation (S)
  if (present(proj3rd)) then
    ABI_MALLOC(SSS_ref,(3,27,3,3,Sym%nsym,6)); SSS_ref(:,:,:,:,:,:)=zero
    do isym=1,Sym%nsym
      do alpha=1,3
        do mu=1,3
          do beta=1,3
            do nu=1,3
              do gama=1,3
                do xi=1,3
                  temp=Sym%S_ref(alpha,mu,isym,1)*Sym%S_ref(beta,nu,isym,1)*Sym%S_ref(gama,xi,isym,1)
                  SSS_ref(alpha,xi+(nu-1)*3+(mu-1)*9,beta,gama,isym,1)=temp !\Psi_efg
                  SSS_ref(alpha,nu+(xi-1)*3+(mu-1)*9,beta,gama,isym,2)=temp !\Psi_egf
                  SSS_ref(alpha,xi+(mu-1)*3+(nu-1)*9,beta,gama,isym,3)=temp !\Psi_feg
                  SSS_ref(alpha,mu+(xi-1)*3+(nu-1)*9,beta,gama,isym,4)=temp !\Psi_fge
                  SSS_ref(alpha,nu+(mu-1)*3+(xi-1)*9,beta,gama,isym,5)=temp !\Psi_gef
                  SSS_ref(alpha,mu+(nu-1)*3+(xi-1)*9,beta,gama,isym,6)=temp !\Psi_gfe
                end do
              end do
            end do
          end do
        end do
      end do
    end do  
  end if 
        
! Compute the invariance under an arbitrary rotation of the system
  ABI_MALLOC(dlevi,(natom,natom,3,3)) ; dlevi(:,:,:,:)=zero
  Levi_Civita(:,:,:)=zero
  Levi_Civita(1,2,3)=+1 ; Levi_Civita(2,3,1)=+1 ; Levi_Civita(3,1,2)=+1
  Levi_Civita(3,2,1)=-1 ; Levi_Civita(1,3,2)=-1 ; Levi_Civita(2,1,3)=-1
  do beta=1,3
    do nu=1,3
      do iatom=1,natom
        do jatom=1,natom
          do gama=1,3
            dlevi(iatom,jatom,beta,nu)=dlevi(iatom,jatom,beta,nu)+&
&             distance(iatom,jatom,gama+1)*Levi_Civita(beta,gama,nu)
          end do
        end do
      end do
    end do
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Compute the constraints !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  open(unit=16,file=trim(InVar%output_prefix)//'constraints.dat')
  ntotcoeff=CoeffMoore%ntotcoeff
  if (present(proj2nd)) then
    ABI_MALLOC(const_rot1st,(3,                                ntotcoeff)); const_rot1st(:,:)        =zero
    ABI_MALLOC(const_rot2nd,(3,3,               natom_unitcell,ntotcoeff)); const_rot2nd(:,:,:,:)    =zero
    ABI_MALLOC(const_dynmat,(3,3,natom_unitcell,natom_unitcell,ntotcoeff)); const_dynmat(:,:,:,:,:)  =zero
    ABI_MALLOC(const_huang ,(3,3,3,3                          ,ntotcoeff)); const_huang(:,:,:,:,:)   =zero
  end if  
!  if (present(proj3rd)) then
!!FB    ABI_MALLOC(const_rot3rd,(3,3,3,       natom_unitcell,natom,ntotcoeff)); const_rot3rd(:,:,:,:,:,:)  =zero
!    ABI_MALLOC(const_rot3rd,(3,3,3,3,3,     natom_unitcell,natom,ntotcoeff)); const_rot3rd(:,:,:,:,:,:,:,:)=zero
!    ABI_MALLOC(const_asr3rd,(8,3,3,3,     natom_unitcell,natom,ntotcoeff)); const_asr3rd(:,:,:,:,:,:,:)=zero
!  end if

! First order only
  do ishell=1,Shell1at%nshell
    write(16,*) ' Compute the 1st order'
    if (Shell1at%neighbours(1,ishell)%n_interactions.eq.0) cycle
    do iatshell=1,Shell1at%neighbours(1,ishell)%n_interactions
      iatom=Shell1at%neighbours(1,ishell)%atomj_in_shell(iatshell) 
      if (iatom.gt.natom_unitcell) exit
      if (iatom.eq.1) cycle
      isym=Shell1at%neighbours(1,ishell)%sym_in_shell(iatshell)
      ncoeff     =Shell1at%ncoeff(ishell)
      ncoeff_prev=Shell1at%ncoeff_prev(ishell)
      do nu=1,3
        do alpha=1,3
          do icoeff=1,ncoeff
!           1/ Rotational invariances (1st order)
            terme=sum(Sym%S_ref(alpha,:,isym,1)*proj1st(:,icoeff,ishell))*dlevi(1,iatom,alpha,nu)
            const_rot1st(nu,icoeff+ncoeff_prev)= &
&           const_rot1st(nu,icoeff+ncoeff_prev)+terme

!           2/ Rotational invariances (for the 2nd order)
            do beta=1,3
              terme=sum(Sym%S_ref(beta,:,isym,1)*proj1st(:,icoeff,ishell))*Levi_Civita(beta,alpha,nu)
              const_rot2nd(alpha,nu,iatom,icoeff+ncoeff_prev)=&
&             const_rot2nd(alpha,nu,iatom,icoeff+ncoeff_prev)+terme                   
              const_rot2nd(alpha,nu,1,icoeff+ncoeff_prev)=&
&             const_rot2nd(alpha,nu,1,icoeff+ncoeff_prev)-terme                   
            end do
          end do    
        end do    
      end do  
    end do !iatshell
  end do !ishell

! First + second order
  if (present(proj2nd)) then
    write(16,*) ' Compute the 2nd order'
    do ishell=1,nshell2at
      do iatom=1,natom_unitcell
        if (Shell2at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell2at%neighbours(iatom,ishell)%n_interactions
          jatom=Shell2at%neighbours(iatom,ishell)%atomj_in_shell(iatshell) 
          if (iatom==jatom) cycle
          isym=Shell2at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
          trans=Shell2at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
          ncoeff     =Shell2at%ncoeff(ishell)
          ncoeff_prev=Shell2at%ncoeff_prev(ishell)+CoeffMoore%ncoeff1st
          iat_mod=mod(jatom+natom_unitcell-1,natom_unitcell)+1
!         1/ Rotational invariances (2nd order). Number of constraints = natom_unitcell*3**2
          do nu=1,3
            do alpha=1,3
              do beta=1,3
                do icoeff=1,ncoeff
                  terme=sum(SS_ref(alpha,:,beta,isym,trans)*proj2nd(:,icoeff,ishell))*dlevi(iatom,jatom,beta,nu)
                  const_rot2nd(alpha,nu,iatom,icoeff+ncoeff_prev)=&
&                 const_rot2nd(alpha,nu,iatom,icoeff+ncoeff_prev)+terme                   
                end do
              end do
            end do    
          end do  
!         2/ Enforce the symetry of the dynamical matrix. Number of constraints = (3*natom_unitcell)**2
!            Note that we are unable to enforce the symetry when iatom=jatom (We have to write the quations)
          do alpha=1,3
            do beta=1,3
              do icoeff=1,ncoeff
                terme=sum(SS_ref(alpha,:,beta,isym,trans)*proj2nd(:,icoeff,ishell))-&
&                     sum(SS_ref(beta,:,alpha,isym,trans)*proj2nd(:,icoeff,ishell))
                const_dynmat(alpha,beta,iatom,iat_mod,icoeff+ncoeff_prev)=&
&               const_dynmat(alpha,beta,iatom,iat_mod,icoeff+ncoeff_prev)+terme                   
              end do
            end do    
          end do  
!         3/ Huang invariances. Number of constraints = 3**4
          do alpha=1,3
            do beta=1,3
              do gama=1,3
                do lambda=1,3
                  do icoeff=1,ncoeff
                    terme=sum(SS_ref(alpha,:,beta,isym,trans)*proj2nd(:,icoeff,ishell))*&
&                             distance(iatom,jatom,gama+1)*&
&                             distance(iatom,jatom,lambda+1)-&
&                         sum(SS_ref(gama,:,lambda,isym,trans)*proj2nd(:,icoeff,ishell))*&
&                             distance(iatom,jatom,alpha+1)*&
&                             distance(iatom,jatom,beta+1)
                    const_huang(alpha,beta,gama,lambda,icoeff+ncoeff_prev)=&
&                   const_huang(alpha,beta,gama,lambda,icoeff+ncoeff_prev)+terme                   
                  end do
                end do
              end do
            end do
          end do
!         4/ Rotational invariances (for the 3nd order). Number of constraints = natom_unitcell*natom*3**3
!FB          if (present(proj3rd).and.(distance(iatom,jatom,1).lt.InVar%Rcut3)) then
!          if (present(proj3rd)) then
!            do alpha=1,3
!              do beta=1,3
!                do gama=1,3
!                  do lambda=1,3
!                    do icoeff=1,ncoeff
!                      terme1=zero ; terme2=zero ; terme3=zero ; terme4=zero ;
!                      if (alpha.eq.lambda) terme1=sum(SS_ref(gama  ,:,beta  ,isym,trans)*proj2nd(:,icoeff,ishell))
!                      if (beta.eq.lambda)  terme2=sum(SS_ref(alpha ,:,gama  ,isym,trans)*proj2nd(:,icoeff,ishell))
!                      if (alpha.eq.gama)   terme3=sum(SS_ref(lambda,:,beta  ,isym,trans)*proj2nd(:,icoeff,ishell))
!                      if (beta.eq.gama)    terme4=sum(SS_ref(alpha ,:,lambda,isym,trans)*proj2nd(:,icoeff,ishell))
!                      if (distance(iatom,jatom,1).lt.InVar%Rcut3) then
!                        const_rot3rd(:,alpha,beta,gama,lambda,iatom,jatom,icoeff+ncoeff_prev)=&
!&                       const_rot3rd(:,alpha,beta,gama,lambda,iatom,jatom,icoeff+ncoeff_prev)+terme1+terme2-terme3-terme4
!                      end if
!                      const_rot3rd(:,alpha,beta,gama,lambda,iatom,iatom,icoeff+ncoeff_prev)=&
!&                     const_rot3rd(:,alpha,beta,gama,lambda,iatom,iatom,icoeff+ncoeff_prev)-terme1-terme2+terme3+terme4
!                    end do
!                  end do
!                end do
!              end do    
!            end do  
!!FB            do nu=1,3
!!FB              do alpha=1,3
!!FB                do beta=1,3
!!FB                  do gama=1,3
!!FB                    do icoeff=1,ncoeff
!!FB                      terme1=sum(SS_ref(alpha,:,beta,isym,trans)*proj2nd(:,icoeff,ishell))*Levi_Civita(gama,alpha,nu)
!!FB                      terme2=sum(SS_ref(alpha,:,gama,isym,trans)*proj2nd(:,icoeff,ishell))*Levi_Civita(gama,beta,nu)
!!FB!FB                      if (iatom.ne.jatom) then
!!FB                        const_rot3rd(alpha,beta,gama,nu,iatom,jatom,icoeff+ncoeff_prev)=&
!!FB&                       const_rot3rd(alpha,beta,gama,nu,iatom,jatom,icoeff+ncoeff_prev)+terme1+terme2                   
!!FB!FB                      else if (iatom.eq.jatom) then
!!FB!FB                        const_rot3rd(alpha,beta,nu,iatom,iatom,icoeff+ncoeff_prev)=&
!!FB!FB&                       const_rot3rd(alpha,beta,nu,iatom,iatom,icoeff+ncoeff_prev)-terme1-terme2                   
!!FB!FB                      end if
!!FB                    end do
!!FB                  end do
!!FB                end do
!!FB              end do    
!!FB            end do  
!          end if
        end do !iatshell
      end do !iatom
    end do !ishell
  end if !Order=2

! Second + third order
!  if (present(proj3rd)) then
!    write(16,*) ' Compute the 3rd order'
!    do ishell=1,nshell3at
!      do iatom=1,natom_unitcell
!        if (Shell3at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
!        do iatshell=1,Shell3at%neighbours(iatom,ishell)%n_interactions
!          jatom=Shell3at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
!          katom=Shell3at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
!          if ((iatom.eq.jatom).and.(jatom.eq.katom)) cycle
!          isym =Shell3at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
!          trans=Shell3at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
!          ncoeff     =Shell3at%ncoeff(ishell)
!          ncoeff_prev=Shell3at%ncoeff_prev(ishell)+CoeffMoore%ncoeff2nd+CoeffMoore%ncoeff1st
!!         1/ Acoustic sum rules (3rd order). Number of constraints = 6*natom_unitcell*natom*2*3**3
!          if (trans==1) then ; trans1=2 ; trans2=6 ; trans3=4 ; trans4=5 ; trans5=3 ; end if
!          if (trans==2) then ; trans1=1 ; trans2=4 ; trans3=6 ; trans4=3 ; trans5=5 ; end if
!          if (trans==3) then ; trans1=4 ; trans2=5 ; trans3=2 ; trans4=6 ; trans5=1 ; end if
!          if (trans==4) then ; trans1=3 ; trans2=2 ; trans3=5 ; trans4=1 ; trans5=6 ; end if
!          if (trans==5) then ; trans1=6 ; trans2=3 ; trans3=1 ; trans4=4 ; trans5=2 ; end if
!          if (trans==6) then ; trans1=5 ; trans2=1 ; trans3=3 ; trans4=2 ; trans5=4 ; end if
!          do alpha=1,3
!            do beta=1,3
!              do gama=1,3
!!FB                write(6,'(a,1x,8(i3,1x))') '    --->',alpha,beta,gama,iatom,jatom,katom,ishell,trans
!                do icoeff=1,ncoeff
!                  terme =sum(SSS_ref(alpha,:,beta,gama,isym,trans )*proj3rd(:,icoeff,ishell))
!                  terme1=sum(SSS_ref(alpha,:,beta,gama,isym,trans1)*proj3rd(:,icoeff,ishell))
!                  terme2=sum(SSS_ref(alpha,:,beta,gama,isym,trans2)*proj3rd(:,icoeff,ishell))
!!                 First acoustic sum rules (i.ne.j)
!                  if (iatom.ne.jatom) then
!                    const_asr3rd(1,alpha,beta,gama,iatom,jatom,icoeff+ncoeff_prev)=&
!&                   const_asr3rd(1,alpha,beta,gama,iatom,jatom,icoeff+ncoeff_prev)+terme
!                    const_asr3rd(2,alpha,beta,gama,iatom,jatom,icoeff+ncoeff_prev)=&
!&                   const_asr3rd(2,alpha,beta,gama,iatom,jatom,icoeff+ncoeff_prev)+terme1
!                    const_asr3rd(3,alpha,beta,gama,iatom,jatom,icoeff+ncoeff_prev)=&
!&                   const_asr3rd(3,alpha,beta,gama,iatom,jatom,icoeff+ncoeff_prev)+terme2
!                  end if
!!FB                  if (iatom.ne.katom) then
!!FB                    const_asr3rd(9,alpha,beta,gama,iatom,katom,icoeff+ncoeff_prev)=&
!!FB&                   const_asr3rd(9,alpha,beta,gama,iatom,katom,icoeff+ncoeff_prev)+terme
!!FB                  end if
!                  if ((iatom.eq.jatom).and.(iatom.ne.katom)) then
!                    terme3=sum(SSS_ref(alpha,:,beta,gama,isym,trans3)*proj3rd(:,icoeff,ishell))
!                    terme4=sum(SSS_ref(alpha,:,beta,gama,isym,trans4)*proj3rd(:,icoeff,ishell))
!                    terme5=sum(SSS_ref(alpha,:,beta,gama,isym,trans5)*proj3rd(:,icoeff,ishell))
!
!                    const_asr3rd(4,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)=&
!&                   const_asr3rd(4,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)+terme
!                    const_asr3rd(4,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)=&
!&                   const_asr3rd(4,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)-terme1
!                    const_asr3rd(5,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)=&
!&                   const_asr3rd(5,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)+terme
!                    const_asr3rd(5,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)=&
!&                   const_asr3rd(5,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)-terme2
!                    const_asr3rd(6,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)=&
!&                   const_asr3rd(6,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)+terme
!                    const_asr3rd(6,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)=&
!&                   const_asr3rd(6,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)-terme3
!                    const_asr3rd(7,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)=&
!&                   const_asr3rd(7,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)+terme
!                    const_asr3rd(7,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)=&
!&                   const_asr3rd(7,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)-terme4
!                    const_asr3rd(8,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)=&
!&                   const_asr3rd(8,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)+terme
!                    const_asr3rd(8,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)=&
!&                   const_asr3rd(8,alpha,beta,gama,iatom,iatom,icoeff+ncoeff_prev)-terme5
!                  end if
!                end do
!              end do
!            end do
!          end do
!!         2/ Rotational invariances (3rd order). Number of constraints = natom_unitcell*natom*2*3**3
!          if (katom==iatom) cycle
!            do alpha=1,3
!              do beta=1,3
!                do gama=1,3
!                  do lambda=1,3
!                    do icoeff=1,ncoeff
!                      terme1=sum(SSS_ref(alpha,:,beta,gama  ,isym,trans )*proj3rd(:,icoeff,ishell))*distance(iatom,katom,lambda+1)
!                      terme2=sum(SSS_ref(alpha,:,beta,lambda,isym,trans )*proj3rd(:,icoeff,ishell))*distance(iatom,katom,gama+1)
!                      const_rot3rd(1,alpha,beta,gama,lambda,iatom,jatom,icoeff+ncoeff_prev)=&
!&                     const_rot3rd(1,alpha,beta,gama,lambda,iatom,jatom,icoeff+ncoeff_prev)+terme1-terme2
!                      terme1=sum(SSS_ref(alpha,:,beta,gama  ,isym,trans1)*proj3rd(:,icoeff,ishell))*distance(iatom,katom,lambda+1)
!                      terme2=sum(SSS_ref(alpha,:,beta,lambda,isym,trans1)*proj3rd(:,icoeff,ishell))*distance(iatom,katom,gama+1)
!                      const_rot3rd(2,alpha,beta,gama,lambda,iatom,jatom,icoeff+ncoeff_prev)=&
!&                     const_rot3rd(2,alpha,beta,gama,lambda,iatom,jatom,icoeff+ncoeff_prev)+terme1-terme2
!                      terme1=sum(SSS_ref(alpha,:,beta,gama  ,isym,trans2)*proj3rd(:,icoeff,ishell))*distance(iatom,katom,lambda+1)
!                      terme2=sum(SSS_ref(alpha,:,beta,lambda,isym,trans2)*proj3rd(:,icoeff,ishell))*distance(iatom,katom,gama+1)
!                      const_rot3rd(3,alpha,beta,gama,lambda,iatom,jatom,icoeff+ncoeff_prev)=&
!&                     const_rot3rd(3,alpha,beta,gama,lambda,iatom,jatom,icoeff+ncoeff_prev)+terme1-terme2
!                    end do
!                  end do
!                end do
!              end do    
!            end do  
!!FB          do nu=1,3
!!FB            do alpha=1,3
!!FB              do beta=1,3
!!FB                do gama=1,3
!!FB                  do icoeff=1,ncoeff
!!FB                    terme=sum(SSS_ref(alpha,:,beta,gama,isym,trans)*proj3rd(:,icoeff,ishell))*dlevi(iatom,katom,gama,nu)
!!FB                    if (iatom.ne.jatom) then
!!FB                      const_rot3rd(alpha,beta,gama,nu,iatom,jatom,icoeff+ncoeff_prev)=&
!!FB&                     const_rot3rd(alpha,beta,gama,nu,iatom,jatom,icoeff+ncoeff_prev)+terme
!!FB!FB                    else if (iatom.eq.jatom) then
!!FB!FB                      const_rot3rd(alpha,beta,nu,iatom,iatom,icoeff+ncoeff_prev)=&
!!FB!FB&                     const_rot3rd(alpha,beta,nu,iatom,iatom,icoeff+ncoeff_prev)-terme
!!FB                    end if 
!!FB                  end do
!!FB                end do
!!FB              end do
!!FB            end do    
!!FB          end do  
!        end do !iatshell
!      end do !iatom   
!    end do !ishell   
!  end if

  if (present(proj2nd)) ABI_FREE(SS_ref)
  !if (present(proj3rd)) ABI_FREE(SSS_ref)
  ABI_FREE(dlevi)

! Reduce the number of constraints by selecting the non-zero equations
  iconst_new=0
  if (present(proj2nd)) then
!   1/ For Rotational invariances (1st order)
    write(16,*) ' ======== Constraints at the 1st order (Rotational Invariances)='
    iconst=0
    ABI_MALLOC(vectin ,(ntotcoeff,CoeffMoore%nconst_1st)) ; vectin(:,:)=zero
    ABI_MALLOC(vectout,(ntotcoeff,CoeffMoore%nconst_1st)) ; vectout(:,:)=zero
    do nu=1,3
      iconst=iconst+1 
      vectin(:,iconst)=const_rot1st(nu,:)
    end do
    call tdep_calc_orthonorm(ntotcoeff,CoeffMoore%nconst_1st,nconst_loc,vectin,vectout)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        write(16,'(500(f10.5,1x))') vectout(:,iconst_loc)
        CoeffMoore%fcoeff(3*natom*InVar%nstep+iconst_new,:)= &
&       CoeffMoore%fcoeff(3*natom*InVar%nstep+iconst_new,:)+ &
&       (3*natom*InVar%nstep)*vectout(:,iconst_loc)
      end do
    end if 
    ABI_FREE(vectin)
    ABI_FREE(vectout)
    write(16,*) ' Number of constraints for the 1st order (Rotational Invariances)=',nconst_loc

!   2/ For Rotational invariances (2nd order)
    write(16,*) ' ======== Constraints at the 2nd order (Rotational Invariances) ='
    iconst=0
    ABI_MALLOC(vectin ,(ntotcoeff,CoeffMoore%nconst_rot2nd)) ; vectin(:,:)=zero
    ABI_MALLOC(vectout,(ntotcoeff,CoeffMoore%nconst_rot2nd)) ; vectout(:,:)=zero
    do iatom=1,natom_unitcell
      do alpha=1,3
        do nu=1,3
          iconst=iconst+1  
          vectin(:,iconst)=const_rot2nd(alpha,nu,iatom,:)
        end do
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,CoeffMoore%nconst_rot2nd,nconst_loc,vectin,vectout)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        write(16,'(500(f10.5,1x))') vectout(:,iconst_loc)
        CoeffMoore%fcoeff(3*natom*InVar%nstep+iconst_new,:)= &
&       CoeffMoore%fcoeff(3*natom*InVar%nstep+iconst_new,:)+ &
&       (3*natom*InVar%nstep)*vectout(:,iconst_loc)
      end do  
    end if
    ABI_FREE(vectin)
    ABI_FREE(vectout)
    write(16,*) ' Number of constraints for the 2nd order (Rotational Invariances)=',nconst_loc

!   3/ For symetry of the dynamical matrix
    write(16,*) ' ======== Constraints at the 2nd order (Dynamical Matrix) ='
    iconst=0
    ABI_MALLOC(vectin ,(ntotcoeff,CoeffMoore%nconst_dynmat)) ; vectin(:,:)=zero
    ABI_MALLOC(vectout,(ntotcoeff,CoeffMoore%nconst_dynmat)) ; vectout(:,:)=zero
    do iatom=1,natom_unitcell
      do jatom=1,natom_unitcell
        do alpha=1,3
          do beta=1,3
            iconst=iconst+1
            vectin(:,iconst)=const_dynmat(alpha,beta,iatom,jatom,:)
          end do
        end do
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,CoeffMoore%nconst_dynmat,nconst_loc,vectin,vectout)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        write(16,'(500(f10.5,1x))') vectout(:,iconst_loc)
        CoeffMoore%fcoeff(3*natom*InVar%nstep+iconst_new,:)= &
&       CoeffMoore%fcoeff(3*natom*InVar%nstep+iconst_new,:)+ &
&       (3*natom*InVar%nstep)*vectout(:,iconst_loc)
      end do  
    end if
    ABI_FREE(vectin)
    ABI_FREE(vectout)
    write(16,*) ' Number of constraints at the 2nd order (Dynamical Matrix)=',nconst_loc

!   4/ For Huang invariances 
    write(16,*) ' ======== Constraints at the 2nd order (Huang) ='
    iconst=0
    ABI_MALLOC(vectin ,(ntotcoeff,CoeffMoore%nconst_huang)) ; vectin(:,:)=zero
    ABI_MALLOC(vectout,(ntotcoeff,CoeffMoore%nconst_huang)) ; vectout(:,:)=zero
    do alpha=1,3
      do beta=1,3
        do gama=1,3
          do lambda=1,3
            iconst=iconst+1
            vectin(:,iconst)=const_huang(alpha,beta,gama,lambda,:)
          end do
        end do
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,CoeffMoore%nconst_huang,nconst_loc,vectin,vectout)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        write(16,'(500(f10.51,1x))') vectout(:,iconst_loc)
        CoeffMoore%fcoeff(3*natom*InVar%nstep+iconst_new,:)= &
&       CoeffMoore%fcoeff(3*natom*InVar%nstep+iconst_new,:)+ &
&       (3*natom*InVar%nstep)*vectout(:,iconst_loc)
      end do  
    end if
    ABI_FREE(vectin)
    ABI_FREE(vectout)
    write(16,*) ' Number of constraints at the 2nd order (Huang)=',nconst_loc
  end if  

!  if (present(proj3rd)) then
!!   1/ For acoustic sum rules (3rd order)
!    write(16,*) ' ======== Constraints at the 3rd order (Acoustic sum rules) ='
!    iconst=0
!    ABI_MALLOC(vectin ,(ntotcoeff,CoeffMoore%nconst_asr3rd)) ; vectin(:,:)=zero
!    ABI_MALLOC(vectout,(ntotcoeff,CoeffMoore%nconst_asr3rd)) ; vectout(:,:)=zero
!    do iatom=1,natom_unitcell
!      do jatom=1,natom
!        do alpha=1,3
!          do beta=1,3
!            do gama=1,3 
!              do ii=1,8
!                iconst=iconst+1
!                vectin(:,iconst)=const_asr3rd(ii,alpha,beta,gama,iatom,jatom,:)
!              end do 
!            end do
!          end do
!        end do
!      end do
!    end do
!    call tdep_calc_orthonorm(ntotcoeff,CoeffMoore%nconst_asr3rd,nconst_loc,vectin,vectout)
!    if (nconst_loc.ne.0) then
!      do iconst_loc=1,nconst_loc
!        iconst_new=iconst_new+1
!        write(16,'(500(f10.5,x))') vectout(:,iconst_loc)
!        CoeffMoore%fcoeff(3*natom*InVar%nstep+iconst_new,:)= &
!&       CoeffMoore%fcoeff(3*natom*InVar%nstep+iconst_new,:)+ &
!&       (3*natom*InVar%nstep)*vectout(:,iconst_loc)
!      end do  
!    end if
!    ABI_FREE(vectin)
!    ABI_FREE(vectout)
!    write(16,*) ' Number of constraints at the 3rd order (Acoustic sum rules)=',nconst_loc
!
!!   2/ For Rotational invariances (3rd order)
!    write(16,*) ' ======== Constraints at the 3rd order (Rotational Invariances) ='
!    iconst=0
!    ABI_MALLOC(vectin ,(ntotcoeff,CoeffMoore%nconst_rot3rd)) ; vectin(:,:)=zero
!    ABI_MALLOC(vectout,(ntotcoeff,CoeffMoore%nconst_rot3rd)) ; vectout(:,:)=zero
!    do iatom=1,natom_unitcell
!      do jatom=1,natom
!        do alpha=1,3
!          do beta=1,3
!            do gama=1,3
!              do lambda=1,3
!                do ii=1,3
!                  iconst=iconst+1
!                  vectin(:,iconst)=const_rot3rd(ii,alpha,beta,gama,lambda,iatom,jatom,:)
!                end do  
!              end do
!            end do
!          end do
!        end do
!!FB        do alpha=1,3
!!FB          do beta=1,3
!!FB            do nu=1,3 
!!FB              iconst=iconst+1
!!FB              vectin(:,iconst)=const_rot3rd(alpha,beta,nu,iatom,jatom,:)
!!FB            end do
!!FB          end do
!!FB        end do
!      end do
!    end do
!    call tdep_calc_orthonorm(ntotcoeff,CoeffMoore%nconst_rot3rd,nconst_loc,vectin,vectout)
!    if (nconst_loc.ne.0) then
!      do iconst_loc=1,nconst_loc
!        iconst_new=iconst_new+1
!        write(16,'(500(f10.5,x))') vectout(:,iconst_loc)
!        CoeffMoore%fcoeff(3*natom*InVar%nstep+iconst_new,:)= &
!&       CoeffMoore%fcoeff(3*natom*InVar%nstep+iconst_new,:)+ &
!&       (3*natom*InVar%nstep)*vectout(:,iconst_loc)
!      end do  
!    end if
!    ABI_FREE(vectin)
!    ABI_FREE(vectout)
!    write(16,*) ' Number of constraints at the 3rd order (Rotational Invariances)=',nconst_loc
!
!  end if
  write(16,*) '================================================================='
  write(16,*) ' Total number of constraints =',iconst_new

! Deallocate constraints  
  CoeffMoore%ntotconst=iconst_new
  if (present(proj2nd)) then
    ABI_FREE(const_rot1st)
    ABI_FREE(const_rot2nd)
    ABI_FREE(const_dynmat)
    ABI_FREE(const_huang)
  end if  
!  if (present(proj3rd)) then
!    ABI_FREE(const_rot3rd)
!    ABI_FREE(const_asr3rd)
!  end if  
  close(16)

end subroutine tdep_calc_constraints

!====================================================================================================
 subroutine tdep_check_constraints(distance,InVar,Phij_NN,Pij_N,nshell3at,&
&                 ftot3,Psij_ref,Shell3at,Sym,ucart) !optional 

  type(Input_Variables_type),intent(in) :: InVar
  integer, intent(in)  :: nshell3at
  double precision, intent(in)  :: Phij_NN(3*InVar%natom,3*InVar%natom)
  double precision, intent(in)  :: Pij_N(3*InVar%natom)
  double precision, intent(in) :: distance(InVar%natom,InVar%natom,4)
  type(Shell_Variables_type),intent(in), optional :: Shell3at
  type(Symetries_Variables_type),intent(in), optional :: Sym
  double precision, intent(in), optional  :: Psij_ref(3,3,3,nshell3at)
  double precision, intent(inout), optional :: ftot3(3*InVar%natom,InVar%nstep)
  double precision, intent(in), optional  :: ucart(3,InVar%natom,InVar%nstep)
  
  integer :: ii,jj,kk,iatom,jatom,katom,isym,trans
  integer :: alpha,beta,gama,lambda
  integer :: ishell,iatshell,istep
  double precision :: norm1
  double precision :: Kroenecker(3,3),Psij_333(3,3,3)
  double precision, allocatable :: asr3(:,:,:,:,:),rot3(:,:,:,:,:)
  integer :: ierr

  ierr = 0;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!! Compute the invariance under an arbitrary rotation of the system !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Kroenecker(:,:)=zero
  Kroenecker(1,1)=1 ; Kroenecker(2,2)=1 ; Kroenecker(3,3)=1 

! Invariance under an arbitrary rotation (first order)'
  do alpha=1,3
    do beta=1,3
      norm1=zero
      do iatom=1,InVar%natom
        norm1=norm1+Pij_N(3*(iatom-1)+alpha)*distance(1,iatom,beta +1)-&
&                   Pij_N(3*(iatom-1)+beta )*distance(1,iatom,alpha+1)
      end do !iatom
      if (abs(norm1).gt.tol8) then
        write(std_out,*) ' BUG : invariance under arbitrary rotation is not fulfilled (order 1)'
        write(std_out,'(a,2(i3,1x),1(e17.10,1x))') 'alpha,beta,Sum=',alpha,beta,norm1
      end if 
    end do !beta
  end do !alpha

! Invariance under an arbitrary rotation (first and second order)'
  do iatom=1,InVar%natom
    do alpha=1,3
      do beta=1,3
        do gama=1,3
          norm1=zero
          do jatom=1,InVar%natom
            norm1=norm1+Phij_NN(3*(iatom-1)+alpha,3*(jatom-1)+beta)*distance(iatom,jatom,gama+1)-&
&                       Phij_NN(3*(iatom-1)+alpha,3*(jatom-1)+gama)*distance(iatom,jatom,beta+1)
          end do !jatom
          norm1=norm1+Pij_N(3*(iatom-1)+beta)*Kroenecker(alpha,gama)&
&                    -Pij_N(3*(iatom-1)+gama)*Kroenecker(alpha,beta)
          if (abs(norm1).gt.tol8) then
            write(std_out,*) ' BUG : invariance under arbitrary rotation is not fulfilled (order 2)'
            write(std_out,'(a,4(i3,1x),1(e17.10,1x))') 'iatom,alpha,beta,gama,Sum =',iatom,alpha,beta,gama,norm1
          end if 
        end do !gama
      end do !beta
    end do !alpha
  end do !iatom  

! THIRD ORDER
  if (present(Psij_ref)) then
    ABI_MALLOC(rot3,(Invar%natom,3,3,3,3)) ; rot3(:,:,:,:,:)=0.d0
    ABI_MALLOC(asr3,(2,Invar%natom,3,3,3)) ; asr3(:,:,:,:,:)=0.d0
    do iatom=1,InVar%natom
      rot3(:,:,:,:,:)=0.d0
      asr3(:,:,:,:,:)=0.d0
      do jatom=1,InVar%natom
        if (distance(iatom,jatom,1).gt.InVar%Rcut3) cycle
!       Compute the rotational invariance (third order)
        do alpha=1,3
          do beta=1,3
            do gama=1,3
              do lambda=1,3
                rot3(jatom,alpha,beta,gama,lambda)=rot3(jatom,alpha,beta,gama,lambda)+&
&                    Phij_NN(3*(iatom-1)+gama  ,3*(jatom-1)+beta  )*Kroenecker(alpha,lambda)+&
&                    Phij_NN(3*(iatom-1)+alpha ,3*(jatom-1)+gama  )*Kroenecker(beta,lambda)-&
&                    Phij_NN(3*(iatom-1)+lambda,3*(jatom-1)+beta  )*Kroenecker(alpha,gama)-&
&                    Phij_NN(3*(iatom-1)+alpha ,3*(jatom-1)+lambda)*Kroenecker(beta,gama)
             if (iatom.le.7) then
!FB	          write(6,'(a,1x,6(i3,1x),1(e17.10,1x))') 'iatom,jatom,alpha,beta,gama,lambda,rot3',&
!FB&		    iatom,jatom,alpha,beta,gama,lambda,&
!FB&                    Phij_NN(3*(iatom-1)+gama  ,3*(jatom-1)+beta)*Kroenecker(alpha,lambda)+&
!FB&                    Phij_NN(3*(iatom-1)+alpha ,3*(jatom-1)+gama)*Kroenecker(beta,lambda)-&
!FB&                    Phij_NN(3*(iatom-1)+lambda,3*(jatom-1)+beta)*Kroenecker(alpha,gama)-&
!FB&                    Phij_NN(3*(iatom-1)+alpha ,3*(jatom-1)+lambda)*Kroenecker(beta,gama)
             end if
              end do !lambda  
            end do !gama
          end do !beta
        end do !alpha
      end do !jatom
      do ishell=1,nshell3at
!       Build the 3x3x3 IFC of an atom in this shell    
        if (Shell3at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell3at%neighbours(iatom,ishell)%n_interactions
          jatom=Shell3at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
          katom=Shell3at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
          isym =Shell3at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
          trans=Shell3at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
          call tdep_build_psij333(isym,InVar,Psij_ref(:,:,:,ishell),Psij_333,Sym,trans) 
!         Calculation of the force components (third order)
          do istep=1,InVar%nstep
            do ii=1,3
              do jj=1,3
                do kk=1,3
                  ftot3(3*(iatom-1)+ii,istep)=ftot3(3*(iatom-1)+ii,istep)+&
&                      Psij_333(ii,jj,kk)*ucart(jj,jatom,istep)*ucart(kk,katom,istep)
                end do  
              end do  
            end do  
          end do
!         Compute the first ASR : sum_k Psi_ijk=0 
!              --> Psi_iji+sum_{k.ne.i} Psi_ijk=0
!              --> if i.eq.j Psi_iii+sum_{k.ne.i} Psi_iik
          asr3(1,jatom,:,:,:)=asr3(1,jatom,:,:,:)+Psij_333(:,:,:)
!         Compute the second ASR : sum_j Psi_ijk=0 
          asr3(2,katom,:,:,:)=asr3(2,katom,:,:,:)+Psij_333(:,:,:)
!         Compute the rotational invariance (third order)
      if (iatom.eq.1.and.jatom.eq.7) then
            write(std_out,'(a,6(i5,1x),2(e12.6,1x))') &
&          'ishell,iatom,jatom,katom,isym,trans,Psi=', &
&           ishell,iatom,jatom,katom,isym,trans,Psij_333(1,3,3),asr3(1,jatom,1,3,3)
      end if
          do alpha=1,3
            do beta=1,3
              do gama=1,3
                do lambda=1,3
                  rot3(jatom,alpha,beta,gama,lambda)=rot3(jatom,alpha,beta,gama,lambda)+&
&                      Psij_333(alpha,beta,gama  )*distance(iatom,katom,lambda+1)-&
&                      Psij_333(alpha,beta,lambda)*distance(iatom,katom,gama  +1)
              if (iatom.le.7) then
!FB	            write(6,'(a,1x,6(i3,1x),1(e17.10,1x))') 'iatom,jatom,alpha,beta,gama,lambda,rot3',&
!FB&		      iatom,jatom,alpha,beta,gama,lambda,&
!FB&                      Psij_333(alpha,beta,gama)*distance(iatom,katom,lambda+1)-&
!FB&                      Psij_333(alpha,beta,lambda)*distance(iatom,katom,gama+1)
!FB	            write(6,'(a,1x,4(e17.10,1x))') &
!FB&		    'Psij_333(alpha,beta,gama  ),distance(iatom,katom,lambda+1),&
!FB&                    Psij_333(alpha,beta,lambda),distance(iatom,katom,gama+1  )',&
!FB&                    Psij_333(alpha,beta,gama  ),distance(iatom,katom,lambda+1),&
!FB&                    Psij_333(alpha,beta,lambda),distance(iatom,katom,gama+1  )
               end if
                end do  
              end do  
            end do  
          end do  
        end do !iatshell
      end do !ishell  
!     Check the acoustic sum rules (third order)
      do ii=1,3
        do jj=1,3
          do kk=1,3
!           Check the first acoustic sum rule
            do jatom=1,InVar%natom
              if (abs(asr3(1,jatom,ii,jj,kk)).gt.tol8) then
                write(std_out,'(a,1x,5(i3,1x),1(e17.10,1x))') ' BUG --->',ii,jj,kk,iatom,jatom,asr3(1,jatom,ii,jj,kk)
                ABI_ERROR('The acoustic sum rule is not fulfilled at the 3rd order (3rd dim)')
              end if
            end do !jatom  
!           Check the second acoustic sum rule
            do katom=1,InVar%natom
              if (abs(asr3(2,katom,ii,jj,kk)).gt.tol8) then
                write(std_out,'(a,1x,5(i3,1x),1(e17.10,1x))') ' BUG --->',ii,jj,kk,iatom,katom,asr3(2,katom,ii,jj,kk)
                ABI_ERROR('The acoustic sum rule is not fulfilled at the 3rd order (2nd dim)')
              end if
            end do !jatom  
          end do !kk
        end do !jj
      end do !ii  
!     Check the rotational invariance (third order)
      do jatom=1,InVar%natom
        do alpha=1,3
          do beta=1,3
            do gama=1,3
              do lambda=1,3
                if (abs(rot3(jatom,alpha,beta,gama,lambda)).gt.tol8) then
                  write(std_out,'(a,6(i3,1x),1(e17.10,1x))') ' BUG ---> iatom,jatom,alpha,beta,gama,lambda,norm =',&
&                                iatom,jatom,alpha,beta,gama,lambda,rot3(jatom,alpha,beta,gama,lambda)
                  ABI_ERROR('The invariance under arbitrary rotation is not fulfilled (order 3)')
                end if  
              end do !lambda 
            end do !gama 
          end do !beta 
        end do !alpha 
      end do !jatom
    end do !iatom  
    ABI_FREE(asr3)
    ABI_FREE(rot3)
  end if !Psij_ref

 end subroutine tdep_check_constraints
!====================================================================================================

end module m_tdep_constraints
