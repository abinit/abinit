
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_constraints

  use defs_basis
  use m_errors
  use m_profiling_abi
  use m_xmpi
  use m_abi_linalg,       only : abi_xorthonormalize  
  use m_tdep_readwrite,   only : Input_Variables_type, MPI_enreg_type
  use m_tdep_sym,         only : Symetries_Variables_type
  use m_tdep_utils,       only : Coeff_Moore_type, Constraints_Variables_type
  use m_tdep_shell,       only : Shell_Variables_type
  use m_tdep_phi3,        only : tdep_build_phi3_333
  use m_tdep_phi4,        only : tdep_build_phi4_3333

  implicit none

  public :: tdep_calc_orthonorm
  public :: tdep_calc_constraints
  public :: tdep_check_constraints

contains 

!====================================================================================================
subroutine tdep_calc_orthonorm(dim1,dim2,nindep,vect)

  integer, intent(in) :: dim1,dim2
  integer, intent(out) :: nindep
  double precision, intent(inout) :: vect(dim1,dim2)

  integer :: ii,jj,kk
  double precision :: prod_scal

! Filter non-zero vectors
  ii=0
  do kk=1,dim2
    if (sum(abs(vect(:,kk))).gt.tol8) then
      ii=ii+1
      vect(:,ii)=vect(:,kk)
    end if 
  end do
  nindep=ii

! Gram-Schmidt orthogonalization
  do kk=2,nindep
    do jj=1,kk-1
      prod_scal=sum(vect(:,jj)*vect(:,jj))
      if (abs(prod_scal).gt.tol8) then
        vect(:,kk)=vect(:,kk)-sum(vect(:,kk)*vect(:,jj))/prod_scal*vect(:,jj)
      end if  
    end do
  end do

! Store the non-zero vectors and normalize
  ii=0
  do kk=1,nindep
    prod_scal=sum(vect(:,kk)*vect(:,kk))
    if (abs(prod_scal).gt.tol8) then
      ii=ii+1
      vect(:,ii)=vect(:,kk)/dsqrt(prod_scal)
    end if  
  end do  
  do kk=nindep+1,dim2
    vect(:,kk)=zero
  end do  
  nindep=ii

end subroutine tdep_calc_orthonorm

!====================================================================================================
subroutine tdep_calc_constraints(CoeffMoore,distance,InVar,MPIdata,nshell1at,nshell2at,nshell3at,nshell4at,&
&                                proj1st,proj2nd,proj3rd,proj4th,Shell1at,Shell2at,Shell3at,Shell4at,Sym) 

  type(Coeff_Moore_type), intent(inout) :: CoeffMoore
  type(Input_Variables_type),intent(in) :: InVar
  type(Symetries_Variables_type),intent(in) :: Sym
  type(MPI_enreg_type), intent(in) :: MPIdata
  type(Shell_Variables_type),intent(in) :: Shell1at
  type(Shell_Variables_type),intent(in) :: Shell2at
  type(Shell_Variables_type),intent(in) :: Shell3at
  type(Shell_Variables_type),intent(in) :: Shell4at
  integer, intent(in) :: nshell1at,nshell2at,nshell3at,nshell4at
  double precision, intent(in) :: proj1st(3 ,3 ,nshell1at)
  double precision, intent(in) :: proj2nd(9 ,9 ,nshell2at)
  double precision, intent(in) :: proj3rd(27,27,nshell3at)
  double precision, intent(in) :: proj4th(81,81,nshell4at)
  double precision, intent(in) :: distance(InVar%natom,InVar%natom,4)

  integer :: ishell,ncoeff,ncoeff_prev,iatom,jatom,katom,latom,iatshell,YES,ishell2at,counter
  integer :: icoeff,iconst,nconst,nconst_loc,iconst_loc,iconst_new,isym,itrans,ntotcoeff,iat_mod
  integer :: mu,nu,xi,zeta,alpha,beta,gama,delta,lambda,natom_unitcell,natom,ii
  double precision :: terme,temp,terme1,terme2,terme3,terme4,terme5
  double precision, allocatable :: SS_ref(:,:,:,:,:)
  double precision, allocatable :: vect(:,:)
  double precision, allocatable :: const_rot1st(:,:,:)
  double precision, allocatable :: const_rot2nd(:,:,:,:,:)
  double precision, allocatable :: const_dynmat(:,:,:,:,:)
  double precision, allocatable :: const_huang(:,:,:,:,:)
  double precision, allocatable :: const_asr4th(:,:,:,:,:,:)
  double precision, allocatable :: const_rot4th(:,:,:,:,:,:,:)
  double precision, allocatable :: cf_ovlp(:,:)
  type(Constraints_Variables_type) :: Const3,Const4
  logical :: order2,order3,order4

  natom_unitcell=InVar%natom_unitcell
  natom         =InVar%natom

  order2 = .false.
  order3 = .false.
  order4 = .false.
  if (InVar%Order.ge.2) order2=.true.
  if (InVar%Order.ge.3) order3=.true.
!FB4th  if (InVar%Order.ge.4) order4=.true.
  if (InVar%Order.ge.4) order4=.false.

  write(InVar%stdout,*) ' '
  write(InVar%stdout,*) '#############################################################################'
  write(InVar%stdout,*) '###################### Compute the constraints ##############################'

! For each couple of atoms, transform the Phi2 (3x3) ifc matrix using the symetry operation (S)
  if (order2.or.order3) then
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
! For each couple of atoms, transform the Phi3 (3x3x3) ifc matrix using the symetry operation (S)
  if (order3.or.order4) then
    ABI_MALLOC(Const3%Sprod,(Sym%nsym,6))
    do isym=1,Sym%nsym
      do itrans=1,6
        ABI_MALLOC(Const3%Sprod(isym,itrans)%SSS,(3,27,3,3)); Const3%Sprod(isym,itrans)%SSS(:,:,:,:)=zero
      end do  
    end do  
    do isym=1,Sym%nsym
      do alpha=1,3
        do mu=1,3
          do beta=1,3
            do nu=1,3
              do gama=1,3
                do xi=1,3
                  temp=Sym%S_ref(alpha,mu,isym,1)*Sym%S_ref(beta,nu,isym,1)*Sym%S_ref(gama,xi,isym,1)
                  Const3%Sprod(isym,1)%SSS(alpha,xi+(nu-1)*3+(mu-1)*9,beta ,gama) =temp !\Psi_efg
                  Const3%Sprod(isym,2)%SSS(alpha,xi+(nu-1)*3+(mu-1)*9,gama ,beta) =temp !\Psi_egf
                  Const3%Sprod(isym,3)%SSS(beta ,xi+(nu-1)*3+(mu-1)*9,alpha,gama) =temp !\Psi_feg
                  Const3%Sprod(isym,4)%SSS(beta ,xi+(nu-1)*3+(mu-1)*9,gama ,alpha)=temp !\Psi_fge
                  Const3%Sprod(isym,5)%SSS(gama ,xi+(nu-1)*3+(mu-1)*9,alpha,beta) =temp !\Psi_gef
                  Const3%Sprod(isym,6)%SSS(gama ,xi+(nu-1)*3+(mu-1)*9,beta ,alpha)=temp !\Psi_gfe
                end do
              end do
            end do
            end do
        end do
      end do
    end do  
  end if 
        
! For each couple of atoms, transform the Phi4 (3x3x3x3) ifc matrix using the symetry operation (S)
  if (order4) then
    ABI_MALLOC(Const4%Sprod,(Sym%nsym,24))
    do isym=1,Sym%nsym
      do itrans=1,24
        ABI_MALLOC(Const4%Sprod(isym,itrans)%SSSS,(3,81,3,3,3)); Const4%Sprod(isym,itrans)%SSSS(:,:,:,:,:)=zero
      end do  
    end do  
    do isym=1,Sym%nsym
      do alpha=1,3
        do mu=1,3
          do beta=1,3
            do nu=1,3
              do gama=1,3
                do xi=1,3
                  do delta=1,3
                    do zeta=1,3
                      counter=zeta+(xi-1)*3+(nu-1)*9+(mu-1)*27
                      temp=Sym%S_ref(alpha,mu,isym,1)*Sym%S_ref(beta  ,nu ,isym,1)*&
&                          Sym%S_ref(gama,xi ,isym,1)*Sym%S_ref(delta,zeta,isym,1)
                      Const4%Sprod(isym,1 )%SSSS(alpha,counter,beta,gama,delta)=temp !\Psi_efgh
                      Const4%Sprod(isym,2 )%SSSS(alpha,counter,gama,beta,delta)=temp !\Psi_egfh
                      Const4%Sprod(isym,3 )%SSSS(beta,counter,alpha,gama,delta)=temp !\Psi_fegh
                      Const4%Sprod(isym,4 )%SSSS(beta,counter,gama,alpha,delta)=temp !\Psi_fgeh
                      Const4%Sprod(isym,5 )%SSSS(gama,counter,alpha,beta,delta)=temp !\Psi_gefh
                      Const4%Sprod(isym,6 )%SSSS(gama,counter,beta,alpha,delta)=temp !\Psi_gfeh
  
                      Const4%Sprod(isym,7 )%SSSS(alpha,counter,beta,delta,gama)=temp !\Psi_efhg
                      Const4%Sprod(isym,8 )%SSSS(alpha,counter,gama,delta,beta)=temp !\Psi_eghf
                      Const4%Sprod(isym,9 )%SSSS(beta,counter,alpha,delta,gama)=temp !\Psi_fehg
                      Const4%Sprod(isym,10)%SSSS(beta,counter,gama,delta,alpha)=temp !\Psi_fghe
                      Const4%Sprod(isym,11)%SSSS(gama,counter,alpha,delta,beta)=temp !\Psi_gehf
                      Const4%Sprod(isym,12)%SSSS(gama,counter,beta,delta,alpha)=temp !\Psi_gfhe
  
                      Const4%Sprod(isym,13)%SSSS(alpha,counter,delta,beta,gama)=temp !\Psi_ehfg
                      Const4%Sprod(isym,14)%SSSS(alpha,counter,delta,gama,beta)=temp !\Psi_ehgf
                      Const4%Sprod(isym,15)%SSSS(beta,counter,delta,alpha,gama)=temp !\Psi_fheg
                      Const4%Sprod(isym,16)%SSSS(beta,counter,delta,gama,alpha)=temp !\Psi_fhge
                      Const4%Sprod(isym,17)%SSSS(gama,counter,delta,alpha,beta)=temp !\Psi_ghef
                      Const4%Sprod(isym,18)%SSSS(gama,counter,delta,beta,alpha)=temp !\Psi_ghfe
  
                      Const4%Sprod(isym,19)%SSSS(delta,counter,alpha,beta,gama)=temp !\Psi_hefg
                      Const4%Sprod(isym,20)%SSSS(delta,counter,alpha,gama,beta)=temp !\Psi_hegf
                      Const4%Sprod(isym,21)%SSSS(delta,counter,beta,alpha,gama)=temp !\Psi_hfeg
                      Const4%Sprod(isym,22)%SSSS(delta,counter,beta,gama,alpha)=temp !\Psi_hfge
                      Const4%Sprod(isym,23)%SSSS(delta,counter,gama,alpha,beta)=temp !\Psi_hgef
                      Const4%Sprod(isym,24)%SSSS(delta,counter,gama,beta,alpha)=temp !\Psi_hgfe
  
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do  
  end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Compute the constraints !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ntotcoeff=CoeffMoore%ntotcoeff
! First order only
  write(InVar%stdout,*) '########################## At the 1st order #################################'
  if (order2) then
    ABI_MALLOC(const_rot1st,(3,3,                              ntotcoeff)); const_rot1st(:,:,:)      =zero
    ABI_MALLOC(const_rot2nd,(3,3,3,             natom_unitcell,ntotcoeff)); const_rot2nd(:,:,:,:,:)  =zero
    ABI_MALLOC(const_dynmat,(3,3,natom_unitcell,natom_unitcell,ntotcoeff)); const_dynmat(:,:,:,:,:)  =zero
    ABI_MALLOC(const_huang ,(3,3,3,3                          ,ntotcoeff)); const_huang(:,:,:,:,:)   =zero
    do ishell=1,Shell1at%nshell
      if (Shell1at%neighbours(1,ishell)%n_interactions.eq.0) cycle
      do iatshell=1,Shell1at%neighbours(1,ishell)%n_interactions
        iatom=Shell1at%neighbours(1,ishell)%atomj_in_shell(iatshell) 
        if (iatom.ge.natom_unitcell) cycle
        if (iatom.eq.1) cycle
        isym=Shell1at%neighbours(1,ishell)%sym_in_shell(iatshell)
        ncoeff     =Shell1at%ncoeff(ishell)
        ncoeff_prev=Shell1at%ncoeff_prev(ishell)
        do alpha=1,3
          do beta=1,3
            do icoeff=1,ncoeff
!             1/ Rotational invariances (1st order)
              terme1=sum(Sym%S_ref(alpha,:,isym,1)*proj1st(:,icoeff,ishell))*distance(1,iatom,beta +1)
              terme2=sum(Sym%S_ref(beta ,:,isym,1)*proj1st(:,icoeff,ishell))*distance(1,iatom,alpha+1)
              const_rot1st(alpha,beta,icoeff+ncoeff_prev)= &
&             const_rot1st(alpha,beta,icoeff+ncoeff_prev)+terme1-terme2
  
!             2/ Rotational invariances (for the 2nd order)
              do gama=1,3
                terme1=zero ; terme2=zero 
                if (alpha.eq.gama) terme1=sum(Sym%S_ref(beta,:,isym,1)*proj1st(:,icoeff,ishell))
                if (alpha.eq.beta) terme2=sum(Sym%S_ref(gama,:,isym,1)*proj1st(:,icoeff,ishell))
                const_rot2nd(alpha,beta,gama,iatom,icoeff+ncoeff_prev)=&
&               const_rot2nd(alpha,beta,gama,iatom,icoeff+ncoeff_prev)+terme1-terme2
                const_rot2nd(alpha,beta,gama,1,icoeff+ncoeff_prev)=&
&               const_rot2nd(alpha,beta,gama,1,icoeff+ncoeff_prev)-terme1+terme2
              end do
            end do    
          end do    
        end do  
      end do !iatshell
    end do !ishell

!   First + second order
    write(InVar%stdout,*) '########################## At the 2nd order #################################'
    do ishell=1,nshell2at
      do iatom=1,natom_unitcell
        if (Shell2at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell2at%neighbours(iatom,ishell)%n_interactions
          jatom=Shell2at%neighbours(iatom,ishell)%atomj_in_shell(iatshell) 
          if (iatom==jatom) cycle
          isym=Shell2at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
          itrans=Shell2at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
          ncoeff     =Shell2at%ncoeff(ishell)
          ncoeff_prev=Shell2at%ncoeff_prev(ishell)+CoeffMoore%ncoeff1st
          iat_mod=mod(jatom+natom_unitcell-1,natom_unitcell)+1
!         1/ Rotational invariances (2nd order). Number of constraints = natom_unitcell*3**2
          do alpha=1,3
            do beta=1,3
              do gama=1,3
                do icoeff=1,ncoeff
                  terme1=sum(SS_ref(alpha,:,beta,isym,itrans)*proj2nd(:,icoeff,ishell))*distance(iatom,jatom,gama+1)
                  terme2=sum(SS_ref(alpha,:,gama,isym,itrans)*proj2nd(:,icoeff,ishell))*distance(iatom,jatom,beta+1)
                  const_rot2nd(alpha,beta,gama,iatom,icoeff+ncoeff_prev)=&
&                 const_rot2nd(alpha,beta,gama,iatom,icoeff+ncoeff_prev)+terme1-terme2                  
                end do
              end do
            end do    
          end do  
!         2/ Enforce the symetry of the dynamical matrix. Number of constraints = (3*natom_unitcell)**2
!            Note that we are unable to enforce the symetry when iatom=jatom (We have to write the equations)
          do alpha=1,3
            do beta=1,3
              do icoeff=1,ncoeff
                terme=sum(SS_ref(alpha,:,beta,isym,itrans)*proj2nd(:,icoeff,ishell))-&
&                     sum(SS_ref(beta,:,alpha,isym,itrans)*proj2nd(:,icoeff,ishell))
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
                    terme=sum(SS_ref(alpha,:,beta,isym,itrans)*proj2nd(:,icoeff,ishell))*&
&                             distance(iatom,jatom,gama+1)*&
&                             distance(iatom,jatom,lambda+1)-&
&                         sum(SS_ref(gama,:,lambda,isym,itrans)*proj2nd(:,icoeff,ishell))*&
&                             distance(iatom,jatom,alpha+1)*&
&                             distance(iatom,jatom,beta+1)
                    const_huang(alpha,beta,gama,lambda,icoeff+ncoeff_prev)=&
&                   const_huang(alpha,beta,gama,lambda,icoeff+ncoeff_prev)+terme                   
                  end do
                end do
              end do
            end do
          end do
        end do !iatshell
      end do !iatom
    end do !ishell
  end if !Order=1,2

! Third order
  if (order3) then
    write(InVar%stdout,*) '########################## At the 3rd order #################################'
    ABI_MALLOC(Const3%AsrRot3,(natom_unitcell,natom,ntotcoeff))
    do iatom=1,natom_unitcell
      do jatom=1,natom
        do icoeff=1,ntotcoeff
          ABI_MALLOC(Const3%AsrRot3(iatom,jatom,icoeff)%ABG,   (3,3,3)); Const3%AsrRot3(iatom,jatom,icoeff)%ABG(:,:,:)   =zero
          ABI_MALLOC(Const3%AsrRot3(iatom,jatom,icoeff)%ABGD,(3,3,3,3)); Const3%AsrRot3(iatom,jatom,icoeff)%ABGD(:,:,:,:)=zero
        end do  
      end do  
    end do  
    do ishell=1,nshell2at
      do iatom=1,natom_unitcell
        if (Shell2at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell2at%neighbours(iatom,ishell)%n_interactions
          jatom=Shell2at%neighbours(iatom,ishell)%atomj_in_shell(iatshell) 
          if (iatom==jatom) cycle
          isym=Shell2at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
          itrans=Shell2at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
          ncoeff     =Shell2at%ncoeff(ishell)
          ncoeff_prev=Shell2at%ncoeff_prev(ishell)+CoeffMoore%ncoeff1st
          iat_mod=mod(jatom+natom_unitcell-1,natom_unitcell)+1
!         1/ Rotational invariances (coming from the 2nd order). Number of constraints = natom_unitcell*natom*3**3
          if (InVar%Order.ge.3) then
            do alpha=1,3
              do beta=1,3
                do gama=1,3
                  do lambda=1,3
                    do icoeff=1,ncoeff
                      terme1=zero ; terme2=zero ; terme3=zero ; terme4=zero ;
                      if (alpha.eq.lambda) terme1=sum(SS_ref(gama  ,:,beta  ,isym,itrans)*proj2nd(:,icoeff,ishell))
                      if (beta.eq.lambda)  terme2=sum(SS_ref(alpha ,:,gama  ,isym,itrans)*proj2nd(:,icoeff,ishell))
                      if (alpha.eq.gama)   terme3=sum(SS_ref(lambda,:,beta  ,isym,itrans)*proj2nd(:,icoeff,ishell))
                      if (beta.eq.gama)    terme4=sum(SS_ref(alpha ,:,lambda,isym,itrans)*proj2nd(:,icoeff,ishell))
                      if (distance(iatom,jatom,1).lt.InVar%Rcut3) then
                        Const3%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)=&
&                       Const3%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)+terme1+terme2-terme3-terme4
                      end if
                      Const3%AsrRot3(iatom,iatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)=&
&                     Const3%AsrRot3(iatom,iatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)-terme1-terme2+terme3+terme4
                    end do
                  end do
                end do
              end do    
            end do  
          end if !proj3rd
        end do !iatshell
      end do !iatom
    end do !ishell
    do ishell=1,nshell3at
      do iatom=1,natom_unitcell
        if (Shell3at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell3at%neighbours(iatom,ishell)%n_interactions
          jatom=Shell3at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
          katom=Shell3at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
          isym =Shell3at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
          itrans=Shell3at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
          ncoeff     =Shell3at%ncoeff(ishell)
          ncoeff_prev=Shell3at%ncoeff_prev(ishell)+CoeffMoore%ncoeff2nd+CoeffMoore%ncoeff1st
!         2/ Acoustic sum rules (3rd order). Number of constraints = natom_unitcell*natom*3**3
          do alpha=1,3
            do beta=1,3
              do gama=1,3
                do icoeff=1,ncoeff
                  terme =sum(Const3%Sprod(isym,itrans)%SSS(alpha,:,beta,gama)*proj3rd(:,icoeff,ishell))
                    Const3%AsrRot3(iatom,katom,icoeff+ncoeff_prev)%ABG(alpha,beta,gama)=&
&                   Const3%AsrRot3(iatom,katom,icoeff+ncoeff_prev)%ABG(alpha,beta,gama)+terme
                end do
              end do
            end do
          end do
!         2/ Rotational invariances (coming from the 3rd order). Number of constraints = natom_unitcell*natom*3**4
          do alpha=1,3
            do beta=1,3
              do gama=1,3
                do lambda=1,3
                  do icoeff=1,ncoeff
                    terme1=sum(Const3%Sprod(isym,itrans)%SSS(alpha,:,beta,gama  )*proj3rd(:,icoeff,ishell))*distance(iatom,katom,lambda+1)
                    terme2=sum(Const3%Sprod(isym,itrans)%SSS(alpha,:,beta,lambda)*proj3rd(:,icoeff,ishell))*distance(iatom,katom,gama+1)
                    Const3%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)=&
&                   Const3%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)+terme1-terme2
                  end do
                end do
              end do
            end do    
          end do  
        end do !iatshell
      end do !iatom   
    end do !ishell   
  end if !Order=3

! Fourth order
  if (order4) then
    write(InVar%stdout,*) '########################## At the 4th order #################################'
    ABI_MALLOC(Const4%AsrRot4,(natom_unitcell,natom,natom,ntotcoeff))
    do iatom=1,natom_unitcell
      do jatom=1,natom
        do katom=1,natom
          do icoeff=1,ntotcoeff
            ABI_MALLOC(Const4%AsrRot4(iatom,jatom,katom,icoeff)%ABGD,   (3,3,3,3)); Const4%AsrRot4(iatom,jatom,katom,icoeff)%ABGD(:,:,:,:)   =zero
!FB            ABI_MALLOC(Const4%AsrRot4(iatom,jatom,katom,icoeff)%ABGDE,(3,3,3,3,3)); Const4%AsrRot4(iatom,jatom,katom,icoeff)%ABGDE(:,:,:,:,:)=zero
          end do  
        end do  
      end do  
    end do  
!FB    do ishell=1,nshell2at
!FB      do iatom=1,natom_unitcell
!FB        if (Shell2at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
!FB        do iatshell=1,Shell2at%neighbours(iatom,ishell)%n_interactions
!FB          jatom=Shell2at%neighbours(iatom,ishell)%atomj_in_shell(iatshell) 
!FB          if (iatom==jatom) cycle
!FB          isym=Shell2at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
!FB          itrans=Shell2at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
!FB          ncoeff     =Shell2at%ncoeff(ishell)
!FB          ncoeff_prev=Shell2at%ncoeff_prev(ishell)+CoeffMoore%ncoeff1st
!FB          iat_mod=mod(jatom+natom_unitcell-1,natom_unitcell)+1
!FB!         1/ Rotational invariances (coming from the 2nd order). Number of constraints = natom_unitcell*natom*3**3
!FB          if (InVar%Order.ge.3) then
!FB            do alpha=1,3
!FB              do beta=1,3
!FB                do gama=1,3
!FB                  do lambda=1,3
!FB                    do icoeff=1,ncoeff
!FB                      terme1=zero ; terme2=zero ; terme3=zero ; terme4=zero ;
!FB                      if (alpha.eq.lambda) terme1=sum(SS_ref(gama  ,:,beta  ,isym,itrans)*proj2nd(:,icoeff,ishell))
!FB                      if (beta.eq.lambda)  terme2=sum(SS_ref(alpha ,:,gama  ,isym,itrans)*proj2nd(:,icoeff,ishell))
!FB                      if (alpha.eq.gama)   terme3=sum(SS_ref(lambda,:,beta  ,isym,itrans)*proj2nd(:,icoeff,ishell))
!FB                      if (beta.eq.gama)    terme4=sum(SS_ref(alpha ,:,lambda,isym,itrans)*proj2nd(:,icoeff,ishell))
!FB                      if (distance(iatom,jatom,1).lt.InVar%Rcut3) then
!FB                        Const4%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)=&
!FB&                       Const4%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)+terme1+terme2-terme3-terme4
!FB                      end if
!FB                      Const4%AsrRot3(iatom,iatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)=&
!FB&                     Const4%AsrRot3(iatom,iatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)-terme1-terme2+terme3+terme4
!FB                    end do
!FB                  end do
!FB                end do
!FB              end do    
!FB            end do  
!FB          end if !proj3rd
!FB        end do !iatshell
!FB      end do !iatom
!FB    end do !ishell
    do ishell=1,nshell4at
      do iatom=1,natom_unitcell
        if (Shell4at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell4at%neighbours(iatom,ishell)%n_interactions
          jatom=Shell4at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
          katom=Shell4at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
          latom=Shell4at%neighbours(iatom,ishell)%atoml_in_shell(iatshell)
          isym =Shell4at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
          itrans=Shell4at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
          ncoeff     =Shell4at%ncoeff(ishell)
          ncoeff_prev=Shell4at%ncoeff_prev(ishell)+CoeffMoore%ncoeff3rd+CoeffMoore%ncoeff2nd+CoeffMoore%ncoeff1st
!         2/ Acoustic sum rules (4th order). Number of constraints = natom_unitcell*natom**2*3**4
          do alpha=1,3
            do beta=1,3
              do gama=1,3
                do delta=1,3
                  do icoeff=1,ncoeff
                    terme =sum(Const4%Sprod(isym,itrans)%SSSS(alpha,:,beta,gama,delta)*proj4th(:,icoeff,ishell))
                      Const4%AsrRot4(iatom,katom,latom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,delta)=&
&                     Const4%AsrRot4(iatom,katom,latom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,delta)+terme
                  end do
                end do
              end do
            end do
          end do
!FB!         2/ Rotational invariances (coming from the 3rd order). Number of constraints = natom_unitcell*natom*3**4
!FB          do alpha=1,3
!FB            do beta=1,3
!FB              do gama=1,3
!FB                do lambda=1,3
!FB                  do icoeff=1,ncoeff
!FB                    terme1=sum(Const4%Sprod(isym,itrans)%SSS(alpha,:,beta,gama  )*proj3rd(:,icoeff,ishell))*distance(iatom,katom,lambda+1)
!FB                    terme2=sum(Const4%Sprod(isym,itrans)%SSS(alpha,:,beta,lambda)*proj3rd(:,icoeff,ishell))*distance(iatom,katom,gama+1)
!FB                    Const4%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)=&
!FB&                   Const4%AsrRot3(iatom,jatom,icoeff+ncoeff_prev)%ABGD(alpha,beta,gama,lambda)+terme1-terme2
!FB                  end do
!FB                end do
!FB              end do
!FB            end do    
!FB          end do  
        end do !iatshell
      end do !iatom   
    end do !ishell   
  end if !Order=4

  if (order2.or.order3) ABI_FREE(SS_ref)
  if (order3.or.order4) then
    do isym=1,Sym%nsym
      do itrans=1,6
        ABI_FREE(Const3%Sprod(isym,itrans)%SSS)
      end do
    end do
    ABI_FREE(Const3%Sprod)
  end if  
  if (order4) then
    do isym=1,Sym%nsym
      do itrans=1,24
        ABI_FREE(Const4%Sprod(isym,itrans)%SSSS)
      end do
    end do
    ABI_FREE(Const4%Sprod)
  end if  
  ABI_MALLOC(CoeffMoore%const ,(CoeffMoore%ntotconst,ntotcoeff)); CoeffMoore%const (:,:)=0.d0 

! Reduce the number of constraints by selecting the non-zero equations
  write(InVar%stdout,*) '################## Reduce the number of constraints #########################'
  iconst_new=0
  if (order2) then
!   1/ For Rotational invariances (1st order)
    iconst=0
    ABI_MALLOC(vect,(ntotcoeff,CoeffMoore%nconst_1st)) ; vect(:,:)=zero
    do alpha=1,3
      do beta=1,3
        iconst=iconst+1 
        vect(:,iconst)=const_rot1st(alpha,beta,:)
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,CoeffMoore%nconst_1st,nconst_loc,vect)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        CoeffMoore%const(iconst_new,:)=vect(:,iconst_loc)
      end do
    end if 
    ABI_FREE(vect)
    ABI_FREE(const_rot1st)
    CoeffMoore%nconst_1st=nconst_loc

!   2/ For Rotational invariances (2nd order)
    iconst=0
    ABI_MALLOC(vect,(ntotcoeff,CoeffMoore%nconst_rot2nd)) ; vect(:,:)=zero
    do iatom=1,natom_unitcell
      do alpha=1,3
        do beta=1,3
          do gama=1,3
            iconst=iconst+1  
            vect(:,iconst)=const_rot2nd(alpha,beta,gama,iatom,:)
          end do
        end do
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,CoeffMoore%nconst_rot2nd,nconst_loc,vect)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        CoeffMoore%const(iconst_new,:)=vect(:,iconst_loc)
      end do  
    end if
    ABI_FREE(vect)
    ABI_FREE(const_rot2nd)
    CoeffMoore%nconst_rot2nd=nconst_loc

!   3/ For symetry of the dynamical matrix
    iconst=0
    ABI_MALLOC(vect,(ntotcoeff,CoeffMoore%nconst_dynmat)) ; vect(:,:)=zero
    do iatom=1,natom_unitcell
      do jatom=1,natom_unitcell
        do alpha=1,3
          do beta=1,3
            iconst=iconst+1
            vect(:,iconst)=const_dynmat(alpha,beta,iatom,jatom,:)
          end do
        end do
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,CoeffMoore%nconst_dynmat,nconst_loc,vect)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        CoeffMoore%const(iconst_new,:)=vect(:,iconst_loc)
      end do  
    end if
    ABI_FREE(vect)
    ABI_FREE(const_dynmat)
    CoeffMoore%nconst_dynmat=nconst_loc

!   4/ For Huang invariances 
    iconst=0
    ABI_MALLOC(vect,(ntotcoeff,CoeffMoore%nconst_huang)) ; vect(:,:)=zero
    do alpha=1,3
      do beta=1,3
        do gama=1,3
          do lambda=1,3
            iconst=iconst+1
            vect(:,iconst)=const_huang(alpha,beta,gama,lambda,:)
          end do
        end do
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,CoeffMoore%nconst_huang,nconst_loc,vect)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        CoeffMoore%const(iconst_new,:)=vect(:,iconst_loc)
      end do  
    end if
    ABI_FREE(vect)
    ABI_FREE(const_huang)
    CoeffMoore%nconst_huang=nconst_loc
    CoeffMoore%nconst_2nd=CoeffMoore%nconst_rot2nd+CoeffMoore%nconst_dynmat+CoeffMoore%nconst_huang
  end if  

  if (order3) then
!   1/ For acoustic sum rules (3rd order)
    iconst=0
    ABI_MALLOC(vect ,(ntotcoeff,CoeffMoore%nconst_asr3rd)) ; vect (:,:)=zero
    do iatom=1,natom_unitcell
      do jatom=1,natom
        do alpha=1,3
          do beta=1,3
            do gama=1,3 
              iconst=iconst+1
              do ii=1,ntotcoeff
                vect(ii,iconst)=Const3%AsrRot3(iatom,jatom,ii)%ABG(alpha,beta,gama)
              end do  
            end do
          end do
        end do
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,CoeffMoore%nconst_asr3rd,nconst_loc,vect)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        CoeffMoore%const(iconst_new,:)=vect(:,iconst_loc)
      end do  
    end if
    ABI_FREE(vect)
    CoeffMoore%nconst_asr3rd=nconst_loc

!   2/ For Rotational invariances (3rd order)
    iconst=0
    ABI_MALLOC(vect ,(ntotcoeff,CoeffMoore%nconst_rot3rd)) ; vect(:,:)=zero
    do iatom=1,natom_unitcell
      do jatom=1,natom
        do alpha=1,3
          do beta=1,3
            do gama=1,3
              do lambda=1,3
                iconst=iconst+1
                do ii=1,ntotcoeff
                  vect(ii,iconst)=Const3%AsrRot3(iatom,jatom,ii)%ABGD(alpha,beta,gama,lambda)
                  end do  
              end do
            end do
          end do
        end do
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,CoeffMoore%nconst_rot3rd,nconst_loc,vect)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        CoeffMoore%const(iconst_new,:)=vect(:,iconst_loc)
      end do  
    end if
    ABI_FREE(vect)
    do iatom=1,natom_unitcell
      do jatom=1,natom
        do icoeff=1,ntotcoeff
          ABI_FREE(Const3%AsrRot3(iatom,jatom,icoeff)%ABG)
          ABI_FREE(Const3%AsrRot3(iatom,jatom,icoeff)%ABGD)
        end do  
      end do  
    end do  
    ABI_FREE(Const3%AsrRot3)
    CoeffMoore%nconst_rot3rd=nconst_loc
    CoeffMoore%nconst_3rd=CoeffMoore%nconst_asr3rd+CoeffMoore%nconst_rot3rd
  end if

  if (order4) then
!   1/ For acoustic sum rules (4th order)
    iconst=0
    ABI_MALLOC(vect ,(ntotcoeff,CoeffMoore%nconst_asr4th)) ; vect (:,:)=zero
    do iatom=1,natom_unitcell
      do jatom=1,natom
        do katom=1,natom
          do alpha=1,3
            do beta=1,3
              do gama=1,3 
                do delta=1,3 
                  iconst=iconst+1
                  do ii=1,ntotcoeff
                    vect(ii,iconst)=Const4%AsrRot4(iatom,jatom,katom,ii)%ABGD(alpha,beta,gama,delta)
                  end do  
                end do
              end do
            end do
          end do
        end do
      end do
    end do
    call tdep_calc_orthonorm(ntotcoeff,CoeffMoore%nconst_asr4th,nconst_loc,vect)
    if (nconst_loc.ne.0) then
      do iconst_loc=1,nconst_loc
        iconst_new=iconst_new+1
        CoeffMoore%const(iconst_new,:)=vect(:,iconst_loc)
      end do  
    end if
    ABI_FREE(vect)
    CoeffMoore%nconst_asr4th=nconst_loc

!FB!   2/ For Rotational invariances (3rd order)
!FB    iconst=0
!FB    ABI_MALLOC(vect ,(ntotcoeff,CoeffMoore%nconst_rot3rd)) ; vect(:,:)=zero
!FB    do iatom=1,natom_unitcell
!FB      do jatom=1,natom
!FB        do alpha=1,3
!FB          do beta=1,3
!FB            do gama=1,3
!FB              do lambda=1,3
!FB                iconst=iconst+1
!FB                vect(:,iconst)=Const4%AsrRot3(iatom,jatom,:)%ABGD(alpha,beta,gama,lambda)
!FB              end do
!FB            end do
!FB          end do
!FB        end do
!FB      end do
!FB    end do
!FB    call tdep_calc_orthonorm(ntotcoeff,CoeffMoore%nconst_rot3rd,nconst_loc,vect)
!FB    if (nconst_loc.ne.0) then
!FB      do iconst_loc=1,nconst_loc
!FB        iconst_new=iconst_new+1
!FB        CoeffMoore%const(iconst_new,:)=vect(:,iconst_loc)
!FB      end do  
!FB    end if
!FB    ABI_FREE(vect)
!FB    CoeffMoore%nconst_rot3rd=nconst_loc
    do iatom=1,natom_unitcell
      do jatom=1,natom
        do katom=1,natom
          do icoeff=1,ntotcoeff
            ABI_FREE(Const4%AsrRot4(iatom,jatom,katom,icoeff)%ABGD)
!FB            ABI_FREE(Const4%AsrRot4(iatom,jatom,katom,icoeff)%ABGDE)
          end do  
        end do  
      end do  
    end do  
    ABI_FREE(Const4%AsrRot4)
    CoeffMoore%nconst_rot4th=0
    CoeffMoore%nconst_4th=CoeffMoore%nconst_asr4th+CoeffMoore%nconst_rot4th
  end if

! Finalize the orthonormalization   
!FB  nconst=iconst_new
!FB  ABI_MALLOC(vect,(ntotcoeff,nconst)) ; vect(:,:)=zero
!FB  do iconst=1,nconst
!FB    vect(:,iconst)=CoeffMoore%const(iconst,:)
!FB  end do
!FB  ABI_FREE(CoeffMoore%const)
!FB  call tdep_calc_orthonorm(ntotcoeff,nconst,nconst_loc,vect)
!FB  CoeffMoore%ntotconst=nconst_loc
!FB  if (nconst_loc.ne.0) then
!FB    ABI_MALLOC(CoeffMoore%const ,(CoeffMoore%ntotconst,ntotcoeff)); CoeffMoore%const (:,:)=0.d0 
!FB    do iconst_loc=1,nconst_loc
!FB      CoeffMoore%const(iconst_loc,:)=vect(:,iconst_loc)
!FB    end do  
!FB  end if
!FB  ABI_FREE(vect)
  CoeffMoore%ntotconst=iconst_new
  if (MPIdata%iam_master) then
    open(unit=16,file=trim(InVar%output_prefix)//'constraints.dat')
    write(16,*) ' ======== Constraints at the 1st order (Rotational Invariances) ========'
    write(16,*) ' Number of constraints =',CoeffMoore%nconst_1st
    write(16,*) ' ======== Constraints at the 2nd order (Rotational Invariances) ========'
    write(16,*) ' Number of constraints =',CoeffMoore%nconst_rot2nd
    write(16,*) ' ======== Constraints at the 2nd order (Dynamical Matrix) =============='
    write(16,*) ' Number of constraints =',CoeffMoore%nconst_dynmat
    write(16,*) ' ======== Constraints at the 2nd order (Huang) ========================='
    write(16,*) ' Number of constraints =',CoeffMoore%nconst_huang
    if (InVar%Order.ge.3) then
      write(16,*) ' ======== Constraints at the 3rd order (Acoustic sum rules) ============'
      write(16,*) ' Number of constraints =',CoeffMoore%nconst_asr3rd
      write(16,*) ' ======== Constraints at the 3rd order (Rotational Invariances) ========'
      write(16,*) ' Number of constraints =',CoeffMoore%nconst_rot3rd
    end if
    if (InVar%Order.ge.4) then
      write(16,*) ' ======== Constraints at the 4th order (Acoustic sum rules) ============'
      write(16,*) ' Number of constraints =',CoeffMoore%nconst_asr4th
      write(16,*) ' ======== Constraints at the 4th order (Rotational Invariances) ========'
      write(16,*) ' Number of constraints =',CoeffMoore%nconst_rot4th
    end if
    write(16,*) ' ======================================================================='
    write(16,*) ' Total number of constraints =',CoeffMoore%ntotconst
    close(16)
  end if  

end subroutine tdep_calc_constraints

!====================================================================================================
 subroutine tdep_check_constraints(distance,InVar,Phi2,Phi1,nshell3at,nshell4at,&
&                 Phi3_ref,Phi4_ref,Shell3at,Shell4at,Sym,ucart)

  type(Input_Variables_type),intent(in) :: InVar
  integer, intent(in)  :: nshell3at,nshell4at
  double precision, intent(in)  :: Phi2(3*InVar%natom,3*InVar%natom)
  double precision, intent(in)  :: Phi1(3*InVar%natom)
  double precision, intent(in) :: distance(InVar%natom,InVar%natom,4)
  type(Shell_Variables_type),intent(in) :: Shell3at
  type(Shell_Variables_type),intent(in) :: Shell4at
  type(Symetries_Variables_type),intent(in) :: Sym
  double precision, intent(in) :: Phi3_ref(3,3,3,nshell3at)
  double precision, intent(in) :: Phi4_ref(3,3,3,3,nshell4at)
  double precision, intent(in) :: ucart(3,InVar%natom,InVar%my_nstep)
  
  integer :: ii,jj,kk,ll,iatom,jatom,katom,latom,isym,itrans
  integer :: alpha,beta,gama,lambda
  integer :: ishell,iatshell,istep
  double precision :: norm1
  double precision :: Kroenecker(3,3),Phi3_333(3,3,3),Phi4_3333(3,3,3,3)
  double precision, allocatable :: asr3(:,:,:,:,:),rot3(:,:,:,:,:)
  double precision, allocatable :: asr4(:,:,:,:,:,:)
  integer :: ierr
  logical :: order2,order3,order4

  ierr = 0;
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Compute the acoustic sum rule and the !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!! invariance under an arbitrary rotation of the system !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Kroenecker(:,:)=zero
  Kroenecker(1,1)=1 ; Kroenecker(2,2)=1 ; Kroenecker(3,3)=1 

  order2 = .false.
  order3 = .false.
  order4 = .false.
  if (InVar%Order.ge.2) order2=.true.
  if (InVar%Order.ge.3) order3=.true.
!FB4th  if (InVar%Order.ge.4) order4=.true.
  if (InVar%Order.ge.4) order4=.false.

  if (order2) then
!   FIRST ORDER
!   Acoustic sum rule (first order)'
    do alpha=1,3
      norm1=zero
      do iatom=1,InVar%natom_unitcell
        norm1=norm1+Phi1((iatom-1)*3+alpha)
      enddo
      if (abs(norm1).gt.tol8) then
        write(std_out,'(a,1(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ASR1) : alpha,Sum =',alpha,norm1
        if (abs(norm1).gt.tol6) then
          MSG_ERROR('The acoustic sum rule is not fulfilled (order 1)')
        end if
      end if
    enddo
!   Invariance under an arbitrary rotation (first order)'
    do alpha=1,3
      do beta=1,3
        norm1=zero
        do iatom=1,InVar%natom_unitcell
          norm1=norm1+Phi1(3*(iatom-1)+alpha)*distance(1,iatom,beta +1)-&
&                     Phi1(3*(iatom-1)+beta )*distance(1,iatom,alpha+1)
        end do !iatom
        if (abs(norm1).gt.tol8) then
          write(std_out,'(a,2(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ROT1) : alpha,beta,Sum=',alpha,beta,norm1
          if (abs(norm1).gt.tol6) then
            MSG_ERROR('The invariance under arbitrary rotation is not fulfilled (order 1)')
          end if 
        end if 
      end do !beta
    end do !alpha

!   SECOND ORDER
!   Acoustic sum rule (second order)
    do iatom=1,InVar%natom
      do alpha=1,3
        do beta=1,3
          norm1=zero
          do jatom=1,InVar%natom
            norm1=norm1+Phi2((iatom-1)*3+alpha,3*(jatom-1)+beta)
          enddo
          if (abs(norm1).gt.tol8) then
            write(std_out,'(a,3(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ASR2) : iatom,alpha,beta,Sum =',iatom,alpha,beta,norm1
            if (abs(norm1).gt.tol6) then
              MSG_ERROR('The acoustic sum rule is not fulfilled (order 2)')
            end if
          end if
        enddo
      enddo
    enddo
!   Invariance under an arbitrary rotation (first and second order)'
    do iatom=1,InVar%natom
      do alpha=1,3
        do beta=1,3
          do gama=1,3
            norm1=zero
            do jatom=1,InVar%natom
              norm1=norm1+Phi2(3*(iatom-1)+alpha,3*(jatom-1)+beta)*distance(iatom,jatom,gama+1)-&
&                         Phi2(3*(iatom-1)+alpha,3*(jatom-1)+gama)*distance(iatom,jatom,beta+1)              
            end do !jatom
            norm1=norm1+Phi1(3*(iatom-1)+beta)*Kroenecker(alpha,gama)&
&                      -Phi1(3*(iatom-1)+gama)*Kroenecker(alpha,beta)
            if (abs(norm1).gt.tol8) then
              write(std_out,'(a,4(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ROT2) : iatom,alpha,beta,gama,Sum =',iatom,alpha,beta,gama,norm1
              if (abs(norm1).gt.tol4) then
                MSG_ERROR('The invariance under arbitrary rotation is not fulfilled (order 2)')
              end if 
            end if 
          end do !gama
        end do !beta
      end do !alpha
    end do !iatom  
  end if !Order=1,2  
  
  if (order3) then
!   THIRD ORDER
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
&                    Phi2(3*(iatom-1)+gama  ,3*(jatom-1)+beta  )*Kroenecker(alpha,lambda)+&
&                    Phi2(3*(iatom-1)+alpha ,3*(jatom-1)+gama  )*Kroenecker(beta,lambda)-&
&                    Phi2(3*(iatom-1)+lambda,3*(jatom-1)+beta  )*Kroenecker(alpha,gama)-&
&                    Phi2(3*(iatom-1)+alpha ,3*(jatom-1)+lambda)*Kroenecker(beta,gama)
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
          itrans=Shell3at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
          call tdep_build_phi3_333(isym,InVar,Phi3_ref(:,:,:,ishell),Phi3_333,Sym,itrans) 
!         Compute the first ASR : sum_k Psi_ijk=0 
!              --> Psi_iji+sum_{k.ne.i} Psi_ijk=0
!              --> if i.eq.j Psi_iii+sum_{k.ne.i} Psi_iik
          asr3(1,jatom,:,:,:)=asr3(1,jatom,:,:,:)+Phi3_333(:,:,:)
!         Compute the second ASR : sum_j Psi_ijk=0 
          asr3(2,katom,:,:,:)=asr3(2,katom,:,:,:)+Phi3_333(:,:,:)
!         Compute the rotational invariance (third order)
          do alpha=1,3
            do beta=1,3
              do gama=1,3
                do lambda=1,3
                  rot3(jatom,alpha,beta,gama,lambda)=rot3(jatom,alpha,beta,gama,lambda)+&
&                      Phi3_333(alpha,beta,gama  )*distance(iatom,katom,lambda+1)-&
&                      Phi3_333(alpha,beta,lambda)*distance(iatom,katom,gama  +1)
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
                write(std_out,'(a,1x,5(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ASR3) --->',ii,jj,kk,iatom,jatom,asr3(1,jatom,ii,jj,kk)
                if (abs(asr3(1,jatom,ii,jj,kk)).gt.tol6)&
&                  MSG_ERROR('The acoustic sum rule is not fulfilled (order 3, 3rd dim)')
              end if
            end do !jatom  
!           Check the second acoustic sum rule
            do katom=1,InVar%natom
              if (abs(asr3(2,katom,ii,jj,kk)).gt.tol8) then
                write(std_out,'(a,1x,5(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ASR3) --->',ii,jj,kk,iatom,katom,asr3(2,katom,ii,jj,kk)
                if (abs(asr3(2,katom,ii,jj,kk)).gt.tol6)&                
&                  MSG_ERROR('The acoustic sum rule is not fulfilled (order 3, 2nd dim)')
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
                  write(std_out,'(a,6(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ROT3) ---> iatom,jatom,alpha,beta,gama,lambda,norm =',&
&                                iatom,jatom,alpha,beta,gama,lambda,rot3(jatom,alpha,beta,gama,lambda)
                  if (abs(rot3(jatom,alpha,beta,gama,lambda)).gt.tol6)& 
&                     MSG_ERROR('The invariance under arbitrary rotation is not fulfilled (order 3)')
                end if  
              end do !lambda 
            end do !gama 
          end do !beta 
        end do !alpha 
      end do !jatom
    end do !iatom  
    ABI_FREE(asr3)
    ABI_FREE(rot3)
  end if !Order=3

! FOURTH ORDER
  if (order4) then
!FB    ABI_MALLOC(rot3,(Invar%natom,3,3,3,3)) ; rot3(:,:,:,:,:)=0.d0
    ABI_MALLOC(asr4,(Invar%natom,InVar%natom,3,3,3,3)) ; asr4(:,:,:,:,:,:)=0.d0
    do iatom=1,InVar%natom
      asr4(:,:,:,:,:,:)=0.d0
!FB      rot3(:,:,:,:,:)=0.d0
!FB      do jatom=1,InVar%natom
!FB        if (distance(iatom,jatom,1).gt.InVar%Rcut3) cycle
!FB!       Compute the rotational invariance (third order)
!FB        do alpha=1,3
!FB          do beta=1,3
!FB            do gama=1,3
!FB              do lambda=1,3
!FB                rot3(jatom,alpha,beta,gama,lambda)=rot3(jatom,alpha,beta,gama,lambda)+&
!FB&                    Phi2(3*(iatom-1)+gama  ,3*(jatom-1)+beta  )*Kroenecker(alpha,lambda)+&
!FB&                    Phi2(3*(iatom-1)+alpha ,3*(jatom-1)+gama  )*Kroenecker(beta,lambda)-&
!FB&                    Phi2(3*(iatom-1)+lambda,3*(jatom-1)+beta  )*Kroenecker(alpha,gama)-&
!FB&                    Phi2(3*(iatom-1)+alpha ,3*(jatom-1)+lambda)*Kroenecker(beta,gama)
!FB              end do !lambda  
!FB            end do !gama
!FB          end do !beta
!FB        end do !alpha
!FB      end do !jatom
      do ishell=1,nshell4at
!       Build the 3x3x3x3 IFC of an atom in this shell    
        if (Shell4at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell4at%neighbours(iatom,ishell)%n_interactions
          jatom=Shell4at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
          katom=Shell4at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
          latom=Shell4at%neighbours(iatom,ishell)%atoml_in_shell(iatshell)
          isym =Shell4at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
          itrans=Shell4at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
          call tdep_build_phi4_3333(isym,InVar,Phi4_ref(:,:,:,:,ishell),Phi4_3333,Sym,itrans) 
!         Compute the first ASR : sum_l Psi_ijkl=0 
          asr4(jatom,katom,:,:,:,:)=asr4(jatom,katom,:,:,:,:)+Phi4_3333(:,:,:,:)
!FB!         Compute the rotational invariance (third order)
!FB          do alpha=1,3
!FB            do beta=1,3
!FB              do gama=1,3
!FB                do lambda=1,3
!FB                  rot3(jatom,alpha,beta,gama,lambda)=rot3(jatom,alpha,beta,gama,lambda)+&
!FB&                      Phi3_333(alpha,beta,gama  )*distance(iatom,katom,lambda+1)-&
!FB&                      Phi3_333(alpha,beta,lambda)*distance(iatom,katom,gama  +1)
!FB                end do  
!FB              end do  
!FB            end do  
!FB          end do  
        end do !iatshell
      end do !ishell  
!     Check the acoustic sum rules (fourth order)
      do ii=1,3
        do jj=1,3
          do kk=1,3
            do ll=1,3
              do jatom=1,InVar%natom
                do katom=1,InVar%natom
                  if (abs(asr4(jatom,katom,ii,jj,kk,ll)).gt.tol8) then
                    write(std_out,'(a,1x,7(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ASR4) --->',ii,jj,kk,ll,iatom,jatom,katom,asr4(jatom,katom,ii,jj,kk,ll)
                    if (abs(asr4(jatom,katom,ii,jj,kk,ll)).gt.tol6)&
&                      MSG_ERROR('The acoustic sum rule is not fulfilled (order 4)')
                  end if
                end do !katom  
              end do !jatom  
            end do !ll
          end do !kk
        end do !jj
      end do !ii  
!FB!     Check the rotational invariance (third order)
!FB      do jatom=1,InVar%natom
!FB        do alpha=1,3
!FB          do beta=1,3
!FB            do gama=1,3
!FB              do lambda=1,3
!FB                if (abs(rot3(jatom,alpha,beta,gama,lambda)).gt.tol8) then
!FB                  write(std_out,'(a,6(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ROT3) ---> iatom,jatom,alpha,beta,gama,lambda,norm =',&
!FB&                                iatom,jatom,alpha,beta,gama,lambda,rot3(jatom,alpha,beta,gama,lambda)
!FB                  if (abs(rot3(jatom,alpha,beta,gama,lambda)).gt.tol6)& 
!FB&                     MSG_ERROR('The invariance under arbitrary rotation is not fulfilled (order 3)')
!FB                end if  
!FB              end do !lambda 
!FB            end do !gama 
!FB          end do !beta 
!FB        end do !alpha 
!FB      end do !jatom
    end do !iatom  
    ABI_FREE(asr4)
!FB    ABI_FREE(rot3)
  end if !Order=4

 end subroutine tdep_check_constraints
!====================================================================================================

end module m_tdep_constraints
