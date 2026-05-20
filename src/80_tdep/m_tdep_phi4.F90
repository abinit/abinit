
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
  use m_tdep_dataset,     only : atdep_dataset_type, MPI_enreg_type
  use m_tdep_latt,        only : Lattice_type
  use m_tdep_shell,       only : Shell_type
  use m_tdep_sym,         only : Symmetries_type
  use m_tdep_sampling,    only : tdep_Sampling_type
  use m_tdep_solver,      only : tdep_Solver_type
  use m_tdep_model,       only : tdep_Model_type
  use m_tdep_constraints, only : Constraints_type

  implicit none

  public :: tdep_calc_ftot4
  public :: tdep_calc_phi4ref
  public :: tdep_write_phi4
  public :: tdep_build_phi4_3333

contains

!====================================================================================================

 subroutine tdep_calc_ftot4(Model,Invar,Shell4at,ucart,Sym)

  type(tdep_Model_type),intent(inout) :: Model
  type(atdep_dataset_type),intent(in) :: Invar
  type(Shell_type),intent(in) :: Shell4at
  type(Symmetries_type),intent(in) :: Sym
  double precision, intent(in)  :: ucart(3,Invar%natom,Invar%my_nstep)

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
        call tdep_build_phi4_3333(isym,Model%Phi4(:,:,:,:,ishell),Phi4_3333,Sym,itrans)
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
&              3*Invar%natom,0.d0,Model%Phi4UiUjUkUl(istep),3*Invar%natom)
    Model%Forces(3*Invar%natom*(istep-1)+1:3*Invar%natom*istep)=&
&   Model%Forces(3*Invar%natom*(istep-1)+1:3*Invar%natom*istep)-ftot4(:,istep)
  end do
  ABI_FREE(ucart_blas)
  ABI_FREE(ftot4)

 end subroutine tdep_calc_ftot4

!=====================================================================================================
subroutine tdep_calc_phi4ref(Solver,Shell4at,Phi4_ref)

  type(tdep_Solver_type),intent(in) :: Solver
  type(Shell_type),intent(in) :: Shell4at
  double precision, intent(inout) :: Phi4_ref(3,3,3,3,Shell4at%nshell)

  integer :: ishell,ncoeff,ncoeff_prev
  integer :: ii,jj,kk,ll,kappa
  double precision, allocatable :: Phi4_coeff(:)

  ABI_CALLOC(Phi4_coeff, (Solver%ncoeff4th))
  Phi4_coeff(:) = Solver%theta(Solver%ncoeff1st+Solver%ncoeff2nd+Solver%ncoeff3rd+1:Solver%ntotcoeff)

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
            Phi4_ref(ii,jj,kk,ll,ishell)=sum(Shell4at%proj(kappa,1:ncoeff,ishell)*Phi4_coeff(ncoeff_prev+1:ncoeff_prev+ncoeff))
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

  ABI_FREE(Phi4_coeff)

end subroutine tdep_calc_phi4ref

!=====================================================================================================
subroutine tdep_write_phi4(distance,Invar,Phi4_ref,Shell4at,Sym)

  type(atdep_dataset_type),intent(in) :: Invar
  type(Symmetries_type),intent(in) :: Sym
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
#if defined FC_NVHPC
            if (itrans == -1) write(std_out, *)"NVHPC freezes here that is fixed by this print statement."
#endif

            write(Invar%stdout,'(a,i2,i2,a)') '  Phi4^{',ii,jj,'kl}='
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

  type(Symmetries_type),intent(in) :: Sym
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
#if defined FC_NVHPC
                 if (itrans == -1) write(std_out, *)"NVHPC freezes here that is fixed by this print statement."
#endif
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
#if defined FC_NVHPC
          if (itrans == -1) write(std_out, *)"NVHPC freezes here that is fixed by this print statement."
#endif

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
