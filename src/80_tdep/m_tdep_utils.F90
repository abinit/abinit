
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_tdep_utils

 use defs_basis
 use m_errors
 use m_abicore
 use m_xmpi
 use m_io_tools
 use m_tdep_dataset,     only : atdep_dataset_type, MPI_enreg_type
 use m_tdep_latt,        only : Lattice_type
 use m_tdep_sym,         only : Symmetries_type
 use m_tdep_shell,       only : Shell_type
 use m_tdep_model,       only : tdep_Model_type
 use m_tdep_phi3,        only : tdep_build_phi3_333
 use m_tdep_phi4,        only : tdep_build_phi4_3333

 implicit none

 public :: tdep_check_constraints
 public :: tdep_print_Aknowledgments

contains

!====================================================================================================

 subroutine tdep_check_constraints(Model,distance,Invar,Sym,Shell3at,Shell4at)

  type(tdep_Model_type),intent(in) :: Model
  type(atdep_dataset_type),intent(in) :: Invar
  double precision, intent(in) :: distance(Invar%natom,Invar%natom,4)
  type(Symmetries_type),intent(in) :: Sym
  type(Shell_type), intent(in) :: Shell3at
  type(Shell_type), intent(in) :: Shell4at

  integer :: ii,jj,kk,ll,iatom,jatom,katom,latom,isym,itrans
  integer :: alpha,beta,gama,lambda
  integer :: ishell,iatshell
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
  if (Invar%order.ge.2) order2=.true.
  if (Invar%order.ge.3) order3=.true.
!FB4th  if (Invar%order.ge.4) order4=.true.
  if (Invar%order.ge.4) order4=.false.

  if (order2) then
!   FIRST ORDER
!   Acoustic sum rule (first order)'
    do alpha=1,3
      norm1=zero
      do iatom=1,Invar%natom_unitcell
        norm1=norm1+Model%Phi1((iatom-1)*3+alpha)
      enddo
      if (abs(norm1).gt.tol8) then
        write(std_out,'(a,1(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ASR1) : alpha,Sum =',alpha,norm1
        if (abs(norm1).gt.tol6) then
          ABI_WARNING('The acoustic sum rule is not fulfilled (order 1)')
        end if
      end if
    enddo
!   Invariance under an arbitrary rotation (first order)'
    do alpha=1,3
      do beta=1,3
        norm1=zero
        do iatom=1,Invar%natom_unitcell
          norm1=norm1+Model%Phi1(3*(iatom-1)+alpha)*distance(1,iatom,beta +1)-&
&                     Model%Phi1(3*(iatom-1)+beta )*distance(1,iatom,alpha+1)
        end do !iatom
        if (abs(norm1).gt.tol8) then
          write(std_out,'(a,2(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ROT1) : alpha,beta,Sum=',alpha,beta,norm1
          if (abs(norm1).gt.tol6) then
            ABI_WARNING('The invariance under arbitrary rotation is not fulfilled (order 1)')
          end if 
        end if 
      end do !beta
    end do !alpha

!   SECOND ORDER
!   Acoustic sum rule (second order)
    do iatom=1,Invar%natom
      do alpha=1,3
        do beta=1,3
          norm1=zero
          do jatom=1,Invar%natom
            norm1=norm1+Model%Phi2%SR((iatom-1)*3+alpha,3*(jatom-1)+beta)
          enddo
          if (abs(norm1).gt.tol8) then
            write(std_out,'(a,3(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ASR2) : iatom,alpha,beta,Sum =',iatom,alpha,beta,norm1
            if (abs(norm1).gt.tol6) then
              ABI_WARNING('The acoustic sum rule is not fulfilled (order 2)')
            end if
          end if
        enddo
      enddo
    enddo
!   Invariance under an arbitrary rotation (first and second order)'
    do iatom=1,Invar%natom
      do alpha=1,3
        do beta=1,3
          do gama=1,3
            norm1=zero
            do jatom=1,Invar%natom
              norm1=norm1+Model%Phi2%SR(3*(iatom-1)+alpha,3*(jatom-1)+beta)*distance(iatom,jatom,gama+1)-&
&                         Model%Phi2%SR(3*(iatom-1)+alpha,3*(jatom-1)+gama)*distance(iatom,jatom,beta+1)              
            end do !jatom
            norm1=norm1+Model%Phi1(3*(iatom-1)+beta)*Kroenecker(alpha,gama)&
&                      -Model%Phi1(3*(iatom-1)+gama)*Kroenecker(alpha,beta)
            if (abs(norm1).gt.tol8) then
              write(std_out,'(a,4(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ROT2) : iatom,alpha,beta,gama,Sum =',&
&                                                        iatom,alpha,beta,gama,norm1
              if (abs(norm1).gt.tol4) then
                ABI_WARNING('The invariance under arbitrary rotation is not fulfilled (order 2)')
              end if 
            end if 
          end do !gama
        end do !beta
      end do !alpha
    end do !iatom  
  end if !order=1,2  
  
  if (order3) then
!   THIRD ORDER
    ABI_MALLOC(rot3,(Invar%natom,3,3,3,3)) ; rot3(:,:,:,:,:)=0.d0
    ABI_MALLOC(asr3,(2,Invar%natom,3,3,3)) ; asr3(:,:,:,:,:)=0.d0
    do iatom=1,Invar%natom
      rot3(:,:,:,:,:)=0.d0
      asr3(:,:,:,:,:)=0.d0
      do jatom=1,Invar%natom
        if (distance(iatom,jatom,1).gt.Invar%rcut3) cycle
!       Compute the rotational invariance (third order)
        do alpha=1,3
          do beta=1,3
            do gama=1,3
              do lambda=1,3
                rot3(jatom,alpha,beta,gama,lambda)=rot3(jatom,alpha,beta,gama,lambda)+&
&                    Model%Phi2%SR(3*(iatom-1)+gama  ,3*(jatom-1)+beta  )*Kroenecker(alpha,lambda)+&
&                    Model%Phi2%SR(3*(iatom-1)+alpha ,3*(jatom-1)+gama  )*Kroenecker(beta,lambda)-&
&                    Model%Phi2%SR(3*(iatom-1)+lambda,3*(jatom-1)+beta  )*Kroenecker(alpha,gama)-&
&                    Model%Phi2%SR(3*(iatom-1)+alpha ,3*(jatom-1)+lambda)*Kroenecker(beta,gama)
              end do !lambda
            end do !gama
          end do !beta
        end do !alpha
      end do !jatom
      do ishell=1,Shell3at%nshell
!       Build the 3x3x3 IFC of an atom in this shell
        if (Shell3at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell3at%neighbours(iatom,ishell)%n_interactions
          jatom=Shell3at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
          katom=Shell3at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
          isym =Shell3at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
          itrans=Shell3at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
          call tdep_build_phi3_333(isym,Model%Phi3(:,:,:,ishell),Phi3_333,Sym,itrans) 
!         Compute the first ASR : sum_k Phi3_ijk=0 
!              --> Phi3_iji+sum_{k.ne.i} Phi3_ijk=0
!              --> if i.eq.j Phi3_iii+sum_{k.ne.i} Phi3_iik
          asr3(1,jatom,:,:,:)=asr3(1,jatom,:,:,:)+Phi3_333(:,:,:)
!         Compute the second ASR : sum_j Phi3_ijk=0 
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
            do jatom=1,Invar%natom
              if (abs(asr3(1,jatom,ii,jj,kk)).gt.tol8) then
                write(std_out,'(a,1x,5(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ASR3) --->',&
&                                                             ii,jj,kk,iatom,jatom,asr3(1,jatom,ii,jj,kk)
                if (abs(asr3(1,jatom,ii,jj,kk)).gt.tol6)&
&                  ABI_WARNING('The acoustic sum rule is not fulfilled (order 3, 3rd dim)')
              end if
            end do !jatom
!           Check the second acoustic sum rule
            do katom=1,Invar%natom
              if (abs(asr3(2,katom,ii,jj,kk)).gt.tol8) then
                write(std_out,'(a,1x,5(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ASR3) --->',&
&                                                             ii,jj,kk,iatom,katom,asr3(2,katom,ii,jj,kk)
                if (abs(asr3(2,katom,ii,jj,kk)).gt.tol6)&                
&                  ABI_WARNING('The acoustic sum rule is not fulfilled (order 3, 2nd dim)')
              end if
            end do !jatom
          end do !kk
        end do !jj
      end do !ii
!     Check the rotational invariance (third order)
      do jatom=1,Invar%natom
        do alpha=1,3
          do beta=1,3
            do gama=1,3
              do lambda=1,3
                if (abs(rot3(jatom,alpha,beta,gama,lambda)).gt.tol8) then
                  write(std_out,'(a,6(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ROT3) ---> iatom,jatom,alpha,beta,gama,lambda,norm =',&
&                                iatom,jatom,alpha,beta,gama,lambda,rot3(jatom,alpha,beta,gama,lambda)
                  if (abs(rot3(jatom,alpha,beta,gama,lambda)).gt.tol6)& 
&                     ABI_WARNING('The invariance under arbitrary rotation is not fulfilled (order 3)')
                end if
              end do !lambda
            end do !gama
          end do !beta
        end do !alpha
      end do !jatom
    end do !iatom
    ABI_FREE(asr3)
    ABI_FREE(rot3)
  end if !order=3

! FOURTH ORDER
  if (order4) then
!FB    ABI_MALLOC(rot3,(Invar%natom,3,3,3,3)) ; rot3(:,:,:,:,:)=0.d0
    ABI_MALLOC(asr4,(Invar%natom,Invar%natom,3,3,3,3)) ; asr4(:,:,:,:,:,:)=0.d0
    do iatom=1,Invar%natom
      asr4(:,:,:,:,:,:)=0.d0
!FB      rot3(:,:,:,:,:)=0.d0
!FB      do jatom=1,Invar%natom
!FB        if (distance(iatom,jatom,1).gt.Invar%rcut3) cycle
!FB!       Compute the rotational invariance (third order)
!FB        do alpha=1,3
!FB          do beta=1,3
!FB            do gama=1,3
!FB              do lambda=1,3
!FB                rot3(jatom,alpha,beta,gama,lambda)=rot3(jatom,alpha,beta,gama,lambda)+&
!FB&                    Model%Phi2%SR(3*(iatom-1)+gama  ,3*(jatom-1)+beta  )*Kroenecker(alpha,lambda)+&
!FB&                    Model%Phi2%SR(3*(iatom-1)+alpha ,3*(jatom-1)+gama  )*Kroenecker(beta,lambda)-&
!FB&                    Model%Phi2%SR(3*(iatom-1)+lambda,3*(jatom-1)+beta  )*Kroenecker(alpha,gama)-&
!FB&                    Model%Phi2%SR(3*(iatom-1)+alpha ,3*(jatom-1)+lambda)*Kroenecker(beta,gama)
!FB              end do !lambda  
!FB            end do !gama
!FB          end do !beta
!FB        end do !alpha
!FB      end do !jatom
      do ishell=1,Shell4at%nshell
!       Build the 3x3x3x3 IFC of an atom in this shell    
        if (Shell4at%neighbours(iatom,ishell)%n_interactions.eq.0) cycle
        do iatshell=1,Shell4at%neighbours(iatom,ishell)%n_interactions
          jatom=Shell4at%neighbours(iatom,ishell)%atomj_in_shell(iatshell)
          katom=Shell4at%neighbours(iatom,ishell)%atomk_in_shell(iatshell)
          latom=Shell4at%neighbours(iatom,ishell)%atoml_in_shell(iatshell)
          isym =Shell4at%neighbours(iatom,ishell)%sym_in_shell(iatshell)
          itrans=Shell4at%neighbours(iatom,ishell)%transpose_in_shell(iatshell)
          call tdep_build_phi4_3333(isym,Model%Phi4(:,:,:,:,ishell),Phi4_3333,Sym,itrans) 
!         Compute the first ASR : sum_l Phi3_ijkl=0 
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
              do jatom=1,Invar%natom
                do katom=1,Invar%natom
                  if (abs(asr4(jatom,katom,ii,jj,kk,ll)).gt.tol8) then
                    write(std_out,'(a,1x,7(i3,1x),1(e17.10,1x))') '>>>>> WARNING (ASR4) --->',&
&                      ii,jj,kk,ll,iatom,jatom,katom,asr4(jatom,katom,ii,jj,kk,ll)
                    if (abs(asr4(jatom,katom,ii,jj,kk,ll)).gt.tol6)&
&                      ABI_WARNING('The acoustic sum rule is not fulfilled (order 4)')
                  end if
                end do !katom  
              end do !jatom  
            end do !ll
          end do !kk
        end do !jj
      end do !ii  
!FB!     Check the rotational invariance (third order)
!FB      do jatom=1,Invar%natom
!FB        do alpha=1,3
!FB          do beta=1,3
!FB            do gama=1,3
!FB              do lambda=1,3
!FB                if (abs(rot3(jatom,alpha,beta,gama,lambda)).gt.tol8) then
!FB                  write(std_out,'(a,6(i3,1x),1(e17.10,1x))') &
!FB                        &'>>>>> WARNING (ROT3) ---> iatom,jatom,alpha,beta,gama,lambda,norm =',&
!FB&                                iatom,jatom,alpha,beta,gama,lambda,rot3(jatom,alpha,beta,gama,lambda)
!FB                  if (abs(rot3(jatom,alpha,beta,gama,lambda)).gt.tol6)& 
!FB&                     ABI_WARNING('The invariance under arbitrary rotation is not fulfilled (order 3)')
!FB                end if  
!FB              end do !lambda 
!FB            end do !gama 
!FB          end do !beta 
!FB        end do !alpha 
!FB      end do !jatom
    end do !iatom  
    ABI_FREE(asr4)
!FB    ABI_FREE(rot3)
  end if !order=4

 end subroutine tdep_check_constraints

!=====================================================================================================

 subroutine tdep_print_Aknowledgments(unt)

  integer, intent(in) :: unt

  write(unt,*) ' '
  write(unt,'(a)') ' #############################################################################'
  write(unt,'(a)') ' ######################### CALCULATION COMPLETED #############################'
  write(unt,'(a)') ' #############################################################################'
  write(unt,'(a)') ' Suggested references for the acknowledgment of ABINIT usage.'
  write(unt,'(a)') ' '
  write(unt,'(a)') ' The users of ABINIT have little formal obligations with respect to the ABINIT group'
  write(unt,'(a)') ' (those specified in the GNU General Public License, http://www.gnu.org/copyleft/gpl.txt).'
  write(unt,'(a)') ' However, it is common practice in the scientific literature,'
  write(unt,'(a)') ' to acknowledge the efforts of people that have made the research possible.'
  write(unt,'(a)') ' In this spirit, please find below suggested citations of work written by ABINIT developers,'
  write(unt,'(a)') ' corresponding to implementations inside of ABINIT that you have used in the present run.'
  write(unt,'(a)') ' Note also that it will be of great value to readers of publications presenting these results,'
  write(unt,'(a)') ' to read papers enabling them to understand the theoretical formalism and details'
  write(unt,'(a)') ' of the ABINIT implementation.'
  write(unt,'(a)') ' For information on why they are suggested, see also https://docs.abinit.org/theory/acknowledgments.'
  write(unt,'(a)') ' '
  write(unt,'(a)') ' [1] a-TDEP: Temperature Dependent Effective Potential for Abinit '
  write(unt,'(a)') ' -- Lattice dynamic properties including anharmonicity'
  write(unt,'(a)') ' F. Bottin, J. Bieder and J. Bouchet, Comput. Phys. Comm. 254, 107301 (2020).' ! [[cite:Bottin2020]]
  write(unt,'(a)') ' Strong suggestion to cite this paper in your publications.'
  write(unt,'(a)') ' '
  write(unt,'(a)') ' [2] Thermal evolution of vibrational properties of alpha-U'
  write(unt,'(a)') ' J. Bouchet and F. Bottin, Phys. Rev. B 92, 174108 (2015).' ! [[cite:Bouchet2015]]
  write(unt,'(a)') ' Strong suggestion to cite this paper in your publications.'
  write(unt,'(a)') ' '
  write(unt,'(a)') ' [3] Lattice dynamics of anharmonic solids from first principles'
  write(unt,'(a)') ' O. Hellman, I.A. Abrikosov and S.I. Simak, Phys. Rev. B 84, 180301(R) (2011).' ! [[cite:Hellman2011]]
  write(unt,'(a)') ' '
  write(unt,'(a)') ' [4] Temperature dependent effective potential method for accurate free energy calculations of solids'
  write(unt,'(a)') ' O. Hellman, P. Steneteg, I.A. Abrikosov and S.I. Simak, Phys. Rev. B 87, 104111 (2013).' ! [[cite:Hellman2013]]

 end subroutine tdep_print_Aknowledgments

!=====================================================================================================

end module m_tdep_utils
