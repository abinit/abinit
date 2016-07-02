!{\src2tex{textfont=tt}}
!!****f* ABINIT/wfd_vnlpsi
!! NAME
!! wfd_vnlpsi
!!
!! FUNCTION
!!  Evaluates Vnl |psi> in reciprocal space.
!!
!! COPYRIGHT
!! Copyright (C) 2010-2016 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! Wfd<wfd_t>=Datatype gathering info on the wavefunctions
!! band=Band index.
!! ik_ibz=K-point index.
!! spin=Spin index
!! npw_k=Number of PW for this k-point (used to dimension output arrays)
!! Cryst<crystal_t>= data type gathering info on symmetries and unit cell
!! GS_hamk<gs_hamiltonian_type)=Information about the Hamiltonian for this k-point.
!!  Note that no check is done for the consistency btw ik_ibz and the content of GS_hamk.
!! [Kext]<Kdata_t>=Datatype storing form factors and tables that depend of the K.
!!   In the standard mode the routine uses the tables stored in Wfd.
!!   Kext is used to calculate Vnl e^{ik'r}|u_k> where k is defined by ik_ibz whereas k' is
!!   the k-point that has been used to initialize GS_hamk and Kext. This option is used for the
!!   Bloch-state-based interpolation.
!!
!! OUTPUT
!!  vnl_psi(2,npw_k*Wfd%nspinor)= <G+k|V_nl e^{ik'r}|u_k>
!!  opaw_psi(2,npw_k*Wfd%nspinor*Wfd%usepaw)  <G+k|1+S|Cnk> e^{ik'r}|u_k>
!!
!! PARENTS
!!      m_shirley
!!
!! CHILDREN
!!      load_k_hamiltonian,mkkpg,nonlop,pawcprj_alloc,pawcprj_free,wfd_get_cprj
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wfd_vnlpsi(Wfd,band,ik_ibz,spin,npw_k,Cryst,GS_hamk,vnl_psi,opaw_psi,Kext)

 use defs_basis
 use defs_datatypes
 use m_errors
 use m_profiling_abi

 use m_crystal,        only : crystal_t
 use m_wfd,            only : wfd_t, kdata_t, wfd_get_cprj
 use m_hamiltonian,    only : gs_hamiltonian_type, load_k_hamiltonian
 use m_pawcprj,        only : pawcprj_type, pawcprj_alloc, pawcprj_free

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wfd_vnlpsi'
 use interfaces_66_nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,ik_ibz,spin,npw_k
 type(crystal_t),intent(in) :: Cryst
 type(gs_hamiltonian_type),intent(inout) :: GS_hamk
 type(Kdata_t),optional,target,intent(in) :: Kext
 type(wfd_t),target,intent(inout) :: Wfd
!arrays
 real(dp),intent(out) :: vnl_psi(2,npw_k*Wfd%nspinor)  ! <G|Vnl|Cnk>
 real(dp),intent(out) :: opaw_psi(2,npw_k*Wfd%nspinor*Wfd%usepaw) ! <G|1+S|Cnk>

!Local variables ------------------------------
!scalars
 integer,parameter :: idir0=0,ider0=0,nnlout0=0,tim_nonlop0=0
 integer :: choice,cpopt,cp_dim,istwf_k,natom,nkpg,nspinor
 integer :: paw_opt,signs
 logical :: use_kptext,ltest
 character(len=500) :: msg
!arrays
 integer,ABI_CONTIGUOUS pointer :: kg_k(:,:)
 real(dp) :: kpoint(3),dum_enlout(0),dummy_lambda(1)
 real(dp),ABI_CONTIGUOUS pointer :: ffnl(:,:,:,:),ph3d(:,:,:)
 real(dp),allocatable :: kpg_k(:,:),vectin(:,:)
 type(pawcprj_type),allocatable :: Cprj(:,:)

!************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(Wfd%nspinor==1,"nspinor==2 not coded")

 use_kptext = PRESENT(Kext)

 natom   = Cryst%natom
 nspinor = Wfd%nspinor

!npw_k   =  Wfd%Kdata(ik_ibz)%npw
 kg_k => Wfd%Kdata(ik_ibz)%kg_k

 if (.not.use_kptext) then
   istwf_k =  Wfd%Kdata(ik_ibz)%istwfk
   kpoint  =  Wfd%kibz(:,ik_ibz)
   ffnl => Wfd%Kdata(ik_ibz)%fnl_dir0der0
   ph3d => Wfd%Kdata(ik_ibz)%ph3d
 else
   istwf_k =  Kext%istwfk
   kpoint  =  gs_hamk%kpt_k
   ffnl => Kext%fnl_dir0der0
   ph3d => Kext%ph3d
   ltest = (istwf_k==1 .and. Wfd%istwfk(ik_ibz)==1)
   if (.not.ltest) then
     write(msg,'(2a,2(a,i0))')&
&      "istwfk and Wfd%istwfk(ik_ibz) must be 1 when Kext is present, however, ",ch10,&
&      "istwfk= ",istwf_k," and Wfd%istwfk= ",Wfd%istwfk(ik_ibz)
     MSG_BUG(msg)
   end if
   ltest = (ik_ibz==1 .and. ALL(ABS(Wfd%kibz(:,ik_ibz))<tol6))
   if (.not.ltest) then
     write(msg,'(3a,i0,a,3(f8.3,1x))')&
&      "ik_ibz must be 1 and kibz should be 0 when Kext is present, however, ",ch10,&
&      "ik_ibz= ",ik_ibz," and kpoint= ",Wfd%kibz(:,ik_ibz)
     MSG_BUG(msg)
   end if
 end if
 !
 ! Input wavefunction coefficients <G|Cnk>
 ABI_MALLOC(vectin,(2,npw_k*nspinor))
 vectin(1,:) = DBLE (Wfd%Wave(band,ik_ibz,spin)%ug)
 vectin(2,:) = AIMAG(Wfd%Wave(band,ik_ibz,spin)%ug)

 signs  = 2  ! => apply the non-local operator to a function in G-space.
 choice = 1  ! => <G|V_nonlocal|vectin>.
 cpopt  =-1
 paw_opt= 0
 if (Wfd%usepaw==1) then
   paw_opt=4 ! both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
   cpopt=3   ! <p_lmn|in> are already in memory
   if (use_kptext) cpopt=-1 ! Since we have to change the k-point <p_lmn|in> are recomputed in nonlocal and not saved
 end if

 cp_dim = ((cpopt+5)/5)
 ABI_DT_MALLOC(Cprj,(natom,nspinor*cp_dim))

 if (cp_dim>0) then
   call pawcprj_alloc(Cprj,0,Wfd%nlmn_sort)
   call wfd_get_cprj(Wfd,band,ik_ibz,spin,Cryst,Cprj,sorted=.TRUE.)
 end if
 !
 ! Compute (k+G) vectors
 ! Not needed here, however they might be stored in Kdata_t to avoid the calculation inside the loop over bands
 !nkpg=3*GS_hamk%nloalg(3)
 nkpg=0
 ABI_MALLOC(kpg_k,(npw_k,nkpg))
 if (nkpg>0) then
   call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
 endif

 call load_k_hamiltonian(GS_hamk,kpt_k=kpoint,npw_k=npw_k,istwf_k=istwf_k,&
&                        kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl,ph3d_k=ph3d)

 call nonlop(choice,cpopt,Cprj,dum_enlout,GS_hamk,idir0,dummy_lambda,Wfd%MPI_enreg,1,nnlout0,&
&            paw_opt,signs,opaw_psi,tim_nonlop0,vectin,vnl_psi)

 ABI_FREE(kpg_k)
 ABI_FREE(vectin)

 if (cp_dim>0) then
   call pawcprj_free(Cprj)
 end if
 ABI_DT_FREE(Cprj)

 DBG_EXIT("COLL")

end subroutine wfd_vnlpsi
!!***
