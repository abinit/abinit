
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_phi4

  use defs_basis
  use m_errors
  use m_abicore
  use m_numeric_tools
  use m_linalg_interfaces
  use m_io_tools
  use m_crystal,          only : crystal_t
  use m_tdep_readwrite,   only : Input_type, MPI_enreg_type
  use m_tdep_latt,        only : Lattice_type
  use m_tdep_shell,       only : Shell_type
  use m_tdep_sym,         only : Symetries_type
  use m_tdep_utils,       only : Coeff_Moore_type, Constraints_type

  implicit none

  public :: tdep_calc_ftot4
  public :: tdep_calc_phi4fcoeff
  public :: tdep_calc_phi4ref
  public :: tdep_write_phi4
  public :: tdep_build_phi4_3333

contains

!====================================================================================================
 subroutine tdep_calc_ftot4(Forces_TDEP,Invar,Phi4_ref,Phi4UiUjUkUl,Shell4at,ucart,Sym) 

  implicit none 

  type(Input_type),intent(in) :: Invar
  type(Shell_type),intent(in) :: Shell4at
  type(Symetries_type),intent(in) :: Sym
  double precision, intent(in)  :: ucart(3,Invar%natom,Invar%my_nstep)
  double precision, intent(in)  :: Phi4_ref(3,3,3,3,Shell4at%nshell)
  double precision, intent(out) :: Phi4UiUjUkUl(Invar%my_nstep)
  double precision, intent(inout) :: Forces_TDEP(3*Invar%natom*Invar%my_nstep)

  integer :: iatom,jatom,katom,latom,isym,itrans,ishell,iatshell
  integer :: ii,jj,kk,ll,istep
  double precision, allocatable :: Phi4_3333(:,:,:,:)
  double precision, allocatable :: ucart_blas(:)
  double precision, allocatable :: ftot4(:,:)

  ABI_MALLOC(Phi4_3333,(3,3,3,3)) ; Phi4_3333(:,:,:,:)=0.d0
  ABI_MALLOC(ftot4,(3*Invar%natom,Invar%my_nstep)); ftot4(:,:)=0.d0
  do iatom=1,Invar%natom
    do ishell=1,Shell4at%nshell
!     Build the 3x3x3x3 IFC of an atom in this shell    
      if (Shell4at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
      do iatshell=1,Shell4at%neighbours(iatom,ishell)%n_interactions
        jatom=Shell4at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
        katom=Shell4at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
        latom=Shell4at%neighbours(iatom,ishell)%atoml_in_shell(iatshell)
        isym =Shell4at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
        itrans=Shell4at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
        call tdep_build_phi4_3333(isym,Phi4_ref(:,:,:,:,ishell),Phi4_3333,Sym,itrans) 
!       Calculation of the force components (third order)
        do istep=1,Invar%my_nstep
          do ii=1,3
            do jj=1,3
              do kk=1,3
                do ll=1,3
                  ftot4(3*(iatom-1)+ii,istep)=ftot4(3*(iatom-1)+ii,istep)+&
&                      Phi4_3333(ii,jj,kk,ll)*ucart(jj,jatom,istep)*ucart(kk,katom,istep)*ucart(ll,latom,istep)
                end do !ll 
              end do !kk 
            end do !jj 
          end do !ii 
        end do !istep
      end do !iatshell
    end do !ishell  
  end do !iatom  
  ABI_FREE(Phi4_3333)
  ftot4(:,:)=ftot4(:,:)/6.d0

  ABI_MALLOC(ucart_blas  ,(3*Invar%natom)) ; ucart_blas  (:)=0.d0
  do istep=1,Invar%my_nstep
    ucart_blas(:)=0.d0 
    do jatom=1,Invar%natom
      do jj=1,3
        ucart_blas(3*(jatom-1)+jj)=ucart(jj,jatom,istep)
      end do
    end do
    call DGEMM('T','N',1,1,3*Invar%natom,1./4.d0,ftot4(:,istep),3*Invar%natom,ucart_blas,&
&              3*Invar%natom,0.d0,Phi4UiUjUkUl(istep),3*Invar%natom)
    Forces_TDEP(3*Invar%natom*(istep-1)+1:3*Invar%natom*istep)=&
&   Forces_TDEP(3*Invar%natom*(istep-1)+1:3*Invar%natom*istep)-ftot4(:,istep)
  end do
  ABI_FREE(ucart_blas)
  ABI_FREE(ftot4)

 end subroutine tdep_calc_ftot4

!====================================================================================================
subroutine tdep_calc_phi4fcoeff(CoeffMoore,Invar,proj,Shell4at,Sym,ucart)

  implicit none

  type(Input_type),intent(in) :: Invar
  type(Symetries_type),intent(in) :: Sym
  type(Shell_type),intent(in) :: Shell4at
  type(Coeff_Moore_type), intent(inout) :: CoeffMoore
  double precision, intent(in) :: ucart(3,Invar%natom,Invar%my_nstep)
  double precision, intent(in) :: proj(81,81,Shell4at%nshell)

  integer :: ishell,ncoeff,ncoeff_prev,istep,iatom,jatom,katom,latom
  integer :: icoeff,isym,iatshell,itrans,counter
  integer :: mu,nu,xi,zeta,alpha,beta,gama,delta,iindex_l,iindex_h
  integer :: ncoeff_prev_l,ncoeff_prev_h
  double precision :: temp
  double precision, allocatable :: SSSS_proj(:,:,:,:,:)
  type(Constraints_type) :: Const

  ABI_MALLOC(Const%Sprod,(Sym%nsym,24))
  do isym=1,Sym%nsym
    do itrans=1,24
      ABI_MALLOC(Const%Sprod(isym,itrans)%SSSS,(3,81,3,3,3)); Const%Sprod(isym,itrans)%SSSS(:,:,:,:,:)=zero
    end do  
  end do  

! For each couple of atoms, transform the Phi4 (3x3x3) ifc matrix using the symetry operation (S)
! Note: iatom=1 is excluded in order to take into account the atomic sum rule (see below)
  do isym=1,Sym%nsym
    do mu=1,3
      do alpha=1,3
        do nu=1,3
          do beta=1,3
            do xi=1,3
              do gama=1,3
                do zeta=1,3
                  do delta=1,3
                    counter=delta+(gama-1)*3+(beta-1)*9+(alpha-1)*27
                    temp=Sym%S_ref(mu,alpha,isym,1)*Sym%S_ref(nu  ,beta ,isym,1)*&
&                        Sym%S_ref(xi,gama ,isym,1)*Sym%S_ref(zeta,delta,isym,1)
                    Const%Sprod(isym,1 )%SSSS(mu,counter,nu,xi,zeta)=temp !\Phi4_efgh
                    Const%Sprod(isym,2 )%SSSS(mu,counter,xi,nu,zeta)=temp !\Phi4_egfh
                    Const%Sprod(isym,3 )%SSSS(nu,counter,mu,xi,zeta)=temp !\Phi4_fegh
                    Const%Sprod(isym,4 )%SSSS(nu,counter,xi,mu,zeta)=temp !\Phi4_fgeh
                    Const%Sprod(isym,5 )%SSSS(xi,counter,mu,nu,zeta)=temp !\Phi4_gefh
                    Const%Sprod(isym,6 )%SSSS(xi,counter,nu,mu,zeta)=temp !\Phi4_gfeh

                    Const%Sprod(isym,7 )%SSSS(mu,counter,nu,zeta,xi)=temp !\Phi4_efhg
                    Const%Sprod(isym,8 )%SSSS(mu,counter,xi,zeta,nu)=temp !\Phi4_eghf
                    Const%Sprod(isym,9 )%SSSS(nu,counter,mu,zeta,xi)=temp !\Phi4_fehg
                    Const%Sprod(isym,10)%SSSS(nu,counter,xi,zeta,mu)=temp !\Phi4_fghe
                    Const%Sprod(isym,11)%SSSS(xi,counter,mu,zeta,nu)=temp !\Phi4_gehf
                    Const%Sprod(isym,12)%SSSS(xi,counter,nu,zeta,mu)=temp !\Phi4_gfhe

                    Const%Sprod(isym,13)%SSSS(mu,counter,zeta,nu,xi)=temp !\Phi4_ehfg
                    Const%Sprod(isym,14)%SSSS(mu,counter,zeta,xi,nu)=temp !\Phi4_ehgf
                    Const%Sprod(isym,15)%SSSS(nu,counter,zeta,mu,xi)=temp !\Phi4_fheg
                    Const%Sprod(isym,16)%SSSS(nu,counter,zeta,xi,mu)=temp !\Phi4_fhge
                    Const%Sprod(isym,17)%SSSS(xi,counter,zeta,mu,nu)=temp !\Phi4_ghef
                    Const%Sprod(isym,18)%SSSS(xi,counter,zeta,nu,mu)=temp !\Phi4_ghfe

                    Const%Sprod(isym,19)%SSSS(zeta,counter,mu,nu,xi)=temp !\Phi4_hefg
                    Const%Sprod(isym,20)%SSSS(zeta,counter,mu,xi,nu)=temp !\Phi4_hegf
                    Const%Sprod(isym,21)%SSSS(zeta,counter,nu,mu,xi)=temp !\Phi4_hfeg
                    Const%Sprod(isym,22)%SSSS(zeta,counter,nu,xi,mu)=temp !\Phi4_hfge
                    Const%Sprod(isym,23)%SSSS(zeta,counter,xi,mu,nu)=temp !\Phi4_hgef
                    Const%Sprod(isym,24)%SSSS(zeta,counter,xi,nu,mu)=temp !\Phi4_hgfe

                  end do
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do  

  write(Invar%stdout,*) ' Compute the coefficients (at the 4th order) used in the Moore-Penrose...'
  do ishell=1,Shell4at%nshell
    do iatom=1,Invar%natom
      if (Shell4at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
      do iatshell=1,Shell4at%neighbours(iatom,ishell)%n_interactions
        jatom=Shell4at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
        katom=Shell4at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
        latom=Shell4at%neighbours(iatom,ishell)%atoml_in_shell(iatshell)
        isym =Shell4at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
        itrans=Shell4at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
        ncoeff     =Shell4at%ncoeff(ishell)
        ncoeff_prev=Shell4at%ncoeff_prev(ishell)+CoeffMoore%ncoeff3rd+CoeffMoore%ncoeff2nd+CoeffMoore%ncoeff1st
        ncoeff_prev_l=ncoeff_prev+1
        ncoeff_prev_h=ncoeff_prev+ncoeff
  
        ABI_MALLOC(SSSS_proj,(3,3,3,3,ncoeff)) ; SSSS_proj(:,:,:,:,:)=zero
        do mu=1,3
          do nu=1,3
            do xi=1,3
              do zeta=1,3
                do icoeff=1,ncoeff
                  SSSS_proj(mu,nu,xi,zeta,icoeff)=DDOT(81,Const%Sprod(isym,itrans)%SSSS(mu,:,nu,xi,zeta),1,proj(:,icoeff,ishell),1)
                end do
              end do
            end do
          end do
        end do
        do istep=1,Invar%my_nstep
          iindex_l=3*(iatom-1)+3*Invar%natom*(istep-1)+1
          iindex_h=3*(iatom-1)+3*Invar%natom*(istep-1)+3
!         F_i^{\mu}(t)=\sum_{\alpha\beta\gamma\delta,jkl,\nu\xi\zeta} S^{\mu\alpha}.S^{\nu\beta}.S^{\xi\gamma}.S^{\zeta\delta}.
!                      \Phi4_{ijkl}^{\alpha\beta\gamma\delta}.u_l^\zeta(t).u_k^\xi(t).u_j^\nu(t)
          do nu=1,3
            do xi=1,3
              do zeta=1,3
                CoeffMoore%fcoeff(iindex_l:iindex_h,ncoeff_prev_l:ncoeff_prev_h)= &
&               CoeffMoore%fcoeff(iindex_l:iindex_h,ncoeff_prev_l:ncoeff_prev_h)+&
&               SSSS_proj(1:3,nu,xi,zeta,1:ncoeff)*ucart(nu,jatom,istep)*ucart(xi,katom,istep)*ucart(zeta,latom,istep)/6.d0 *&
&               Invar%weights(istep)
              end do  
            end do  
          end do  
        end do !istep
        ABI_FREE(SSSS_proj)
      end do !iatshell
    end do !iatom
  end do !ishell
  write(Invar%stdout,*) ' ------- achieved'
  do isym=1,Sym%nsym
    do itrans=1,24
      ABI_FREE(Const%Sprod(isym,itrans)%SSSS)
    end do
  end do
  ABI_FREE(Const%Sprod)

end subroutine tdep_calc_phi4fcoeff

!=====================================================================================================
subroutine tdep_calc_phi4ref(ntotcoeff,proj,Phi4_coeff,Phi4_ref,Shell4at)

  implicit none

  type(Shell_type),intent(in) :: Shell4at
  integer,intent(in) :: ntotcoeff
  double precision, intent(in) :: proj(81,81,Shell4at%nshell)
  double precision, intent(in) :: Phi4_coeff(ntotcoeff,1)
  double precision, intent(inout) :: Phi4_ref(3,3,3,3,Shell4at%nshell)

  integer :: ishell,ncoeff,ncoeff_prev
  integer :: ii,jj,kk,ll,kappa

  do ishell=1,Shell4at%nshell
!   Build the 3x3x3x3 IFC per shell    
    ncoeff     =Shell4at%ncoeff(ishell)
    ncoeff_prev=Shell4at%ncoeff_prev(ishell)
    kappa=0  
    do ii=1,3
      do jj=1,3
        do kk=1,3
          do ll=1,3
            kappa=kappa+1
            Phi4_ref(ii,jj,kk,ll,ishell)=sum(proj(kappa,1:ncoeff,ishell)*Phi4_coeff(ncoeff_prev+1:ncoeff_prev+ncoeff,1))
          end do  
        end do  
      end do  
    end do  
!  Remove the rounding errors before writing (for non regression testing purposes)
    do ii=1,3
      do jj=1,3
        do kk=1,3
          do ll=1,3
            if (abs(Phi4_ref(ii,jj,kk,ll,ishell)).lt.tol8) Phi4_ref(ii,jj,kk,ll,ishell)=zero
          end do
        end do
      end do
    end do  
  end do  

end subroutine tdep_calc_phi4ref

!=====================================================================================================
subroutine tdep_write_phi4(distance,Invar,Phi4_ref,Shell4at,Sym)

  implicit none

  type(Input_type),intent(in) :: Invar
  type(Symetries_type),intent(in) :: Sym
  type(Shell_type),intent(in) :: Shell4at
  double precision, intent(in) :: distance(Invar%natom,Invar%natom,4)
  double precision, intent(in) :: Phi4_ref(3,3,3,3,Shell4at%nshell)

  integer :: ishell,isym,jatom,katom,latom
  integer :: iatref,jatref,katref,latref,iatshell,itrans
  integer :: ii,jj,kk
  double precision :: tmp1,tmp2,tmp3
  double precision, allocatable :: Phi4_3333(:,:,:,:)

  write(Invar%stdout,*) ' '
  write(Invar%stdout,*) '#############################################################################'
  write(Invar%stdout,*) '#### For each shell, list of coefficients (IFC), number of neighbours... ####'
  write(Invar%stdout,*) '#############################################################################'

! Write the IFCs in the data.out file (with others specifications: 
! number of atoms in a shell, Trace...)
  ABI_MALLOC(Phi4_3333,(3,3,3,3)) ; Phi4_3333(:,:,:,:)=0.d0
  do ishell=1,Shell4at%nshell
    iatref=Shell4at%iatref(ishell)
    if (Shell4at%neighbours(iatref,ishell)%n_interactions.ne.0) then
      jatref=Shell4at%jatref(ishell)
      katref=Shell4at%katref(ishell)
      latref=Shell4at%latref(ishell)
      write(Invar%stdout,'(a,i4,a,i4,a)') ' ======== NEW SHELL (ishell=',ishell,&
&           '): There are',Shell4at%neighbours(iatref,ishell)%n_interactions,' atoms on this shell'
      do iatshell=1,Shell4at%neighbours(iatref,ishell)%n_interactions
        jatom =Shell4at%neighbours(iatref,ishell)%atomj_in_shell(iatshell)
        katom =Shell4at%neighbours(iatref,ishell)%atomk_in_shell(iatshell)
        latom =Shell4at%neighbours(iatref,ishell)%atoml_in_shell(iatshell)
        isym  =Shell4at%neighbours(iatref,ishell)%sym_in_shell(iatshell)
        itrans=Shell4at%neighbours(iatref,ishell)%transpose_in_shell(iatshell)
        call tdep_build_phi4_3333(isym,Phi4_ref(:,:,:,:,ishell),Phi4_3333,Sym,itrans) 
        write(Invar%stdout,'(a,i4,a,i4)') '  For iatcell=',iatref,' ,with type=',mod(iatref-1,Invar%natom_unitcell)+1
        write(Invar%stdout,'(a,i4,a,i4)') '  For jatom  =',jatom ,' ,with type=',mod(jatom -1,Invar%natom_unitcell)+1
        write(Invar%stdout,'(a,i4,a,i4)') '  For katom  =',katom ,' ,with type=',mod(katom -1,Invar%natom_unitcell)+1
        write(Invar%stdout,'(a,i4,a,i4)') '  For latom  =',latom ,' ,with type=',mod(latom -1,Invar%natom_unitcell)+1
        do ii=1,3
          do jj=1,3
            write(Invar%stdout,'(a,i2,i2,a)') '  \Phi4^{',ii,jj,'kl}='
            do kk=1,3
              if (abs(Phi4_3333(ii,jj,1,kk)).lt.5.d-7) then
                tmp1=0.d0
              else
                tmp1=Phi4_3333(ii,jj,1,kk)
              end if
              if (abs(Phi4_3333(ii,jj,2,kk)).lt.5.d-7) then
                tmp2=0.d0
              else
                tmp2=Phi4_3333(ii,jj,2,kk)
              end if
              if (abs(Phi4_3333(ii,jj,3,kk)).lt.5.d-7) then
                tmp3=0.d0
              else
                tmp3=Phi4_3333(ii,jj,3,kk)
              end if
              write(Invar%stdout,'(2x,3(f9.6,1x))') tmp1,tmp2,tmp3
            end do
          end do
        end do
        write(Invar%stdout,'(a,3(f9.6,1x))') '  (i,j) vector components:', (distance(iatref,jatom,jj+1),jj=1,3)
        write(Invar%stdout,'(a,3(f9.6,1x))') '  (j,k) vector components:', (distance(jatom ,katom,jj+1),jj=1,3)
        write(Invar%stdout,'(a,3(f9.6,1x))') '  (j,k) vector components:', (distance(katom ,latom,jj+1),jj=1,3)
        write(Invar%stdout,'(a,3(f9.6,1x))') '  (k,i) vector components:', (distance(latom,iatref,jj+1),jj=1,3)
        write(Invar%stdout,*) ' '
      end do !iatshell
    end if !n_interactions 
  end do !ishell 
  ABI_FREE(Phi4_3333)

end subroutine tdep_write_phi4

!=====================================================================================================
subroutine tdep_build_phi4_3333(isym,Phi4_ref,Phi4_3333,Sym,itrans) 

  implicit none

  type(Symetries_type),intent(in) :: Sym
  double precision, intent(in) :: Phi4_ref(3,3,3,3)
  double precision, intent(out) :: Phi4_3333(3,3,3,3)
  integer,intent(in) :: isym,itrans

  integer :: alpha,beta,gama,delta
  integer :: ii,jj,kk,ll,ee,ff,gg,hh,mu,nu,xi,zeta
  double precision :: Phi4_tmp(3,3,3,3)


! Transform in the new basis wrt S_ref
  Phi4_3333(:,:,:,:)=zero
  do mu=1,3
    do alpha=1,3
      do nu=1,3
        do beta=1,3
          do xi=1,3
            do gama=1,3
              do zeta=1,3
                do delta=1,3
                  Phi4_3333(mu,nu,xi,zeta)=Phi4_3333(mu,nu,xi,zeta)+&
&                 Sym%S_ref(mu,alpha,isym,1)*Sym%S_ref(nu  ,beta ,isym,1)*&
&                 Sym%S_ref(xi,gama ,isym,1)*Sym%S_ref(zeta,delta,isym,1)*Phi4_ref(alpha,beta,gama,delta)
                end do
              end do
            end do
          end do
        end do
      end do
    end do
  end do
  
! Take into account the 6 allowed permutations
  Phi4_tmp(:,:,:,:)=Phi4_3333(:,:,:,:)
  if ((itrans.lt.1).or.(itrans.gt.24)) then  
    ABI_BUG('This value of the symmetry index is not permitted')
  end if
  do ii=1,3
    do jj=1,3
      do kk=1,3
        do ll=1,3
          if (itrans==1) then ; ee=ii ; ff=jj ; gg=kk ; hh=ll ; endif !\Phi4_ijkl
          if (itrans==2) then ; ee=ii ; ff=kk ; gg=jj ; hh=ll ; endif !\Phi4_ikjl
          if (itrans==3) then ; ee=jj ; ff=ii ; gg=kk ; hh=ll ; endif !\Phi4_jikl
          if (itrans==4) then ; ee=jj ; ff=kk ; gg=ii ; hh=ll ; endif !\Phi4_jkil
          if (itrans==5) then ; ee=kk ; ff=ii ; gg=jj ; hh=ll ; endif !\Phi4_kijl
          if (itrans==6) then ; ee=kk ; ff=jj ; gg=ii ; hh=ll ; endif !\Phi4_kjil

          if (itrans==7 ) then ; ee=ii ; ff=jj ; gg=ll ; hh=kk ; endif !\Phi4_ijlk
          if (itrans==8 ) then ; ee=ii ; ff=kk ; gg=ll ; hh=jj ; endif !\Phi4_iklj
          if (itrans==9 ) then ; ee=jj ; ff=ii ; gg=ll ; hh=kk ; endif !\Phi4_jilk
          if (itrans==10) then ; ee=jj ; ff=kk ; gg=ll ; hh=ii ; endif !\Phi4_jkli
          if (itrans==11) then ; ee=kk ; ff=ii ; gg=ll ; hh=jj ; endif !\Phi4_kilj
          if (itrans==12) then ; ee=kk ; ff=jj ; gg=ll ; hh=ii ; endif !\Phi4_kjli

          if (itrans==13) then ; ee=ii ; ff=ll ; gg=jj ; hh=kk ; endif !\Phi4_iljk
          if (itrans==14) then ; ee=ii ; ff=ll ; gg=kk ; hh=jj ; endif !\Phi4_ilkj
          if (itrans==15) then ; ee=jj ; ff=ll ; gg=ii ; hh=kk ; endif !\Phi4_jlik
          if (itrans==16) then ; ee=jj ; ff=ll ; gg=kk ; hh=ii ; endif !\Phi4_jlki
          if (itrans==17) then ; ee=kk ; ff=ll ; gg=ii ; hh=jj ; endif !\Phi4_klij
          if (itrans==18) then ; ee=kk ; ff=ll ; gg=jj ; hh=ii ; endif !\Phi4_klji

          if (itrans==19) then ; ee=ll ; ff=ii ; gg=jj ; hh=kk ; endif !\Phi4_lijk
          if (itrans==20) then ; ee=ll ; ff=ii ; gg=kk ; hh=jj ; endif !\Phi4_likj
          if (itrans==21) then ; ee=ll ; ff=jj ; gg=ii ; hh=kk ; endif !\Phi4_ljik
          if (itrans==22) then ; ee=ll ; ff=jj ; gg=kk ; hh=ii ; endif !\Phi4_ljki
          if (itrans==23) then ; ee=ll ; ff=kk ; gg=ii ; hh=jj ; endif !\Phi4_lkij
          if (itrans==24) then ; ee=ll ; ff=kk ; gg=jj ; hh=ii ; endif !\Phi4_lkji

          Phi4_3333(ee,ff,gg,hh)=Phi4_tmp(ii,jj,kk,ll)
        end do
      end do
    end do  
  end do          

end subroutine tdep_build_phi4_3333

!=====================================================================================================
end module m_tdep_phi4
