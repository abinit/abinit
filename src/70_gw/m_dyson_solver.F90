!!****m* ABINIT/m_dyson_solver
!! NAME
!!  m_dyson_solver
!!
!! FUNCTION
!!  This module contains procedures to solve the Dyson equation to find QP energies.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_dyson_solver

 use defs_basis
 use m_xmpi
 use m_errors
 use m_abicore
 use m_dtfil

 use m_time,          only : timab
 use m_gwdefs,        only : sigparams_t
 use m_numeric_tools, only : linfit, pade, dpade, newrap_step
 use m_io_tools,      only : open_file
 use m_fstrings,      only : int2char10
 use m_hide_lapack,   only : xheev
 use m_bz_mesh,       only : kmesh_t, get_BZ_item
 use m_sigma,         only : sigma_t

 implicit none

 private

 public :: solve_dyson     ! Solve the Dyson equation for the QP energies.

 integer,private,parameter :: NR_MAX_NITER=1000
  ! Max no of iterations in the Newton-Raphson method.

 real(dp),private,parameter :: NR_ABS_ROOT_ERR=0.0001/Ha_eV
  ! Tolerance on the absolute error on the Newton-Raphson root.

CONTAINS  !====================================================================
!!***

!!****f* m_dyson_solver/solve_dyson
!! NAME
!! solve_dyson
!!
!! FUNCTION
!!  Solve the Dyson equation for the QP energies. Two different methods are coded:
!!  The first one is based on the standard perturbative approach in which the self-energy
!!  is linearly expanded around the previous single-particle energy (KS energy if one-shot)
!!  and the derivative is evaluated by finite differences.
!!  In the second method (AC), the values of the self-energy operator on the real axis are obtained
!!  by means of an analitic continuation based on the Pade extrapolation.
!!
!! INPUTS
!!  ikcalc=Index of the considered k-point in the Sigp%kptgw2bz array.
!!  nomega_sigc=Number of frequencies used to evaluate the correlation part of Sigma.
!!  Sigp<sigparams_t>=Structure gathering parameters on the calculation of Sigma.
!!     %minbnd and %maxbnd= min and Max band index for GW correction (for this k-point)
!!     %gwcalctyp=Type of the GW calculation.
!!     %mbpt_sciss=Scissor energy
!!  Sr<sigma_t>=Structure containing the matrix elements of the self-energy INOUT
!!     %nbnds=Number of bands in G0.
!!     %nsppol=Number of independent spin polarizations.
!!     %nsig_ab=Numner of components in the self-energy operator.
!!     %nomega_r=Number of real frequencies for spectral function.
!!     %nomega4sd=Number of real frequencies used to evalute the derivative of Sigma.
!!     %nomega_i=Number of imaginary frequencies for AC.
!!     %omega_i=Purely imaginary frequencies for AC.
!!  Kmesh<kmesh_t>=Info on the K-mesh for the wavefunctions.
!!     %nkibz=Number of points in the IBZ
!!  sigcme_tmp=(nomega_sigc,ib1:ib2,ib1:ib2,nsppol)=Matrix elements of Sigma_c.
!!  qp_ene(nbnds,nkibz,nsppol)= KS or QP energies, only used in case of calculation with scissor operator.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  Sr<sigma_t>=Structure containing the matrix elements of the self-energy:
!!     %sigxme(ib1:ib2,jkibz,nsspol)=Diagonal elements of Sigma_x
!!     %sigcmee0(ib1:ib2,jkibz,nsppol)=Matrix elements of Sigma_c at the initial energy E0.
!!     %dsigmee0(jb,ib1:ib2,nsppol)=Derivate of sigma at the energy E0.
!!     %ze0(ib1:ib2,jkibz,is)=Renormalization factor at the energy E0.
!!     %degw(ib1:ib2,jkibz,is)= QP correction  i.e DeltaE_GW=E-E0
!!     %egw(ib1:ib2,jkibz,is)=QP energy
!!     %sigmee(ib1:ib2,jkibz,is)=Self-energy evaluated at the QP energy.
!!     %sigcme (ib1:ib2,jkibz,io,is)= Sigma_c as a function of frequency.
!!     %sigxcme(ib1:ib2,jkibz,io,is)= Sigma_xc as a function of frequency.
!!     %sigcme4sd (ib1:ib2,jkibz,io,is)= Diagonal matrix elements of \Sigma_c  at frequencies around the KS eigenvalue
!!     %sigxcme4sd(ib1:ib2,jkibz,io,is)= Diagonal matrix elements of \Sigma_xc at frequencies around the KS eigenvalue
!!    where ib1 and ib2 are the band indeces included in the GW calculation for this k-point.
!!
!! PARENTS
!!      m_sigma_driver
!!
!! CHILDREN
!!      int2char10,wrtout
!!
!! SOURCE

subroutine solve_dyson(ikcalc,minbnd,maxbnd,nomega_sigc,Sigp,Kmesh,sigcme_tmp,qp_ene,Sr,prtvol,Dtfil,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikcalc,nomega_sigc,prtvol,minbnd,maxbnd,comm
 type(kmesh_t),intent(in) :: Kmesh
 type(Datafiles_type),intent(in) :: Dtfil
 type(sigparams_t),intent(in) :: Sigp
 type(sigma_t),intent(inout) :: Sr
!arrays
 real(dp),intent(in) :: qp_ene(Sr%nbnds,Sr%nkibz,Sr%nsppol)
 complex(dpc),intent(in) :: sigcme_tmp(nomega_sigc,minbnd:maxbnd,minbnd:maxbnd,Sigp%nsppol*Sigp%nsig_ab)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: iab,ib1,ib2,ikbz_gw,io,ioe0j,spin,is_idx,isym,iter,itim,jb
 integer :: sk_ibz,kb,ld_matrix,mod10,nsploop,my_rank
 real(dp) :: alpha,beta,smrt
 complex(dpc) :: ctdpc,dct,dsigc,sigc,zz,phase
 logical :: converged,ltest
 character(len=500) :: msg
!arrays
 real(dp) :: kbz_gw(3),tsec(2)
 real(dp),allocatable :: e0pde(:),eig(:),scme(:)
 complex(dpc),allocatable :: hdp(:,:),tmpcdp(:),hhartree(:,:,:),htotal(:,:,:),h_tmp1(:,:),h_tmp2(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(490,1,tsec) ! csigme(Dyson)

 my_rank = xmpi_comm_rank(comm)

 mod10=MOD(Sigp%gwcalctyp,10)

 ltest=(nomega_sigc==Sr%nomega_r+Sr%nomega4sd)
 if (mod10==1) ltest=(nomega_sigc==Sr%nomega_i)
 ABI_CHECK(ltest,'Wrong number of frequencies')

 ! Index of the KS or QP energy.
 ioe0j=Sr%nomega4sd/2+1

 ! min and Max band index for GW corrections (for this k-point).
 ib1=MINVAL(Sigp%minbnd(ikcalc,:))
 ib2=MAXVAL(Sigp%maxbnd(ikcalc,:))

 ! Find the index of the k-point for sigma in the IBZ array.
 ikbz_gw=Sigp%kptgw2bz(ikcalc)
 call get_BZ_item(Kmesh,ikbz_gw,kbz_gw,sk_ibz,isym,itim,phase)

 sigc=czero; dsigc=czero

 ! ===========================================================
 ! ==== Solve the Dyson Equation and store results in Sr% ====
 ! ===========================================================

 if (mod10/=1) then
   ! ===============================
   ! ==== Perturbative approach ====
   ! ===============================
   do spin=1,Sr%nsppol
     do jb=ib1,ib2
       ! === Get matrix elements of Sigma_c at energy E0 ===
       ! * SigC(w) is linearly interpolated and the slope alpha is assumed as dSigC/dE
       do iab=1,Sr%nsig_ab
         is_idx=spin; if (Sr%nsig_ab>1) is_idx=iab

         Sr%sigcmee0(jb,sk_ibz,is_idx) = sigcme_tmp(Sr%nomega_r+ioe0j,jb,jb,is_idx)

         ABI_MALLOC(scme,(Sr%nomega4sd))
         ABI_MALLOC(e0pde,(Sr%nomega4sd))
         e0pde(:) = Sr%omega4sd(jb,sk_ibz,:,spin)
         scme(:)  = REAL(sigcme_tmp(Sr%nomega_r+1:Sr%nomega_r+Sr%nomega4sd,jb,jb,is_idx))

         if (Sr%nomega4sd==1) then
           smrt = zero; alpha = zero
         else
           smrt = linfit(Sr%nomega4sd,e0pde(:),scme(:),alpha,beta)
         end if

         if (smrt>0.1/Ha_eV) then
           write(msg,'(3a,i0,a,i0,2a,2(f22.15,2a))')&
             'Values of Re Sig_c are not linear ',ch10,&
             'band index = ',jb,' spin|component = ',is_idx,ch10,&
             'root mean square= ',smrt,ch10,&
             'estimated slope = ',alpha,ch10,&
             'Omega [eV] SigC [eV]'
           MSG_WARNING(msg)
           do io=1,Sr%nomega4sd
             write(msg,'(2f8.4)')e0pde(io)*Ha_eV,scme(io)*Ha_eV
             call wrtout(std_out,msg,"COLL")
           end do
         end if

         ABI_FREE(scme)
         ABI_FREE(e0pde)
         !
         ! === Evaluate renormalization factor and QP correction ===
         ! * Z=(1-dSigma/domega(E0))^-1
         ! * DeltaE_GW=E-E0= (Sigma(E0)-V_xc)/(1-dSigma/domega)
         ! * If nspinor==2, this part is done at the end.
         !
         Sr%dsigmee0(jb,sk_ibz,is_idx)=CMPLX(alpha,zero)

         if (Sr%nsig_ab==1) then
           Sr%ze0(jb,sk_ibz,spin)=one/(one-Sr%dsigmee0(jb,sk_ibz,spin))

           if (ABS(Sigp%mbpt_sciss) < tol6) then
             Sr%degw(jb,sk_ibz,spin) = Sr%ze0(jb,sk_ibz,spin) * &
&              (Sr%sigxme(jb,sk_ibz,spin) + Sr%sigcmee0(jb,sk_ibz,spin) - Sr%e0(jb,sk_ibz,spin) + &
&               Sr%hhartree(jb,jb,sk_ibz,spin))

             Sr%egw(jb,sk_ibz,spin) = Sr%e0(jb,sk_ibz,spin) + Sr%degw(jb,sk_ibz,spin)

             ! Estimate Sigma at the QP-energy: Sigma(E_qp)=Sigma(E0)+(E_qp-E0)*dSigma/dE
             Sr%sigmee(jb,sk_ibz,spin)= &
&              Sr%sigxme(jb,sk_ibz,spin)+Sr%sigcmee0(jb,sk_ibz,spin)+Sr%degw(jb,sk_ibz,spin)*Sr%dsigmee0(jb,sk_ibz,spin)

           else
             ! If GW+scissor: e0 is replaced by qp_ene which contains the updated energy eigenvalue
             Sr%degw(jb,sk_ibz,spin)= Sr%ze0(jb,sk_ibz,spin) * &
&              (Sr%sigxme(jb,sk_ibz,spin) + Sr%sigcmee0(jb,sk_ibz,spin) - qp_ene(jb,sk_ibz,spin) + &
&               Sr%hhartree(jb,jb,sk_ibz,spin))

             Sr%egw(jb,sk_ibz,spin) = qp_ene(jb,sk_ibz,spin) + Sr%degw(jb,sk_ibz,spin)

             ! Estimate Sigma at the QP-energy: Sigma(E_qp)=Sigma(E0)+(E_qp-E0)*dSigma/dE
             Sr%sigmee(jb,sk_ibz,spin)= &
&              Sr%sigxme(jb,sk_ibz,spin) + Sr%sigcmee0(jb,sk_ibz,spin) + &
&              Sr%degw(jb,sk_ibz,spin) * Sr%dsigmee0(jb,sk_ibz,spin)

             ! RS: In the output, the gw corr with respect to e0 without mbpt_sciss is reported.
             Sr%degw(jb,sk_ibz,spin) = Sr%egw(jb,sk_ibz,spin) - Sr%e0(jb,sk_ibz,spin)
           end if
         end if !Sigp%nsig_ab==1

         ! Spectrum of Sigma
         do io=1,Sr%nomega_r
           Sr%sigcme (jb,sk_ibz,io,is_idx)= sigcme_tmp(io,jb,jb,is_idx)
           Sr%sigxcme(jb,sk_ibz,io,is_idx)= Sr%sigxme(jb,sk_ibz,is_idx)+Sr%sigcme(jb,sk_ibz,io,is_idx)
         end do
         do io=1,Sr%nomega4sd
           Sr%sigcme4sd (jb,sk_ibz,io,is_idx)= sigcme_tmp(Sr%nomega_r+io,jb,jb,is_idx)
           Sr%sigxcme4sd(jb,sk_ibz,io,is_idx)= Sr%sigxme(jb,sk_ibz,is_idx)+Sr%sigcme4sd(jb,sk_ibz,io,is_idx)
         end do

       end do !iab

       if (Sr%nsig_ab>1) then
         ABI_CHECK(ABS(Sigp%mbpt_sciss)<0.1d-4,'Scissor with spinor not coded')
         !TODO this should be allocated with nsppol, recheck this part

         ! Evaluate renormalization factor and QP correction.
         ! Z=(1-dSigma/domega(E0))^-1
         ! DeltaE_GW=E-E0= (Sigma(E0)-V_xc)/(1-dSigma/domega)
         !write(std_out,'(a,i2,10f8.3)')' Correlation',jb,Sr%sigcmee0(jb,sk_ibz,:)*Ha_eV,SUM(Sr%sigcmee0(jb,sk_ibz,:))*Ha_eV

         Sr%ze0 (jb,sk_ibz,1) = one/(one-SUM(Sr%dsigmee0(jb,sk_ibz,:)))

         Sr%degw(jb,sk_ibz,1) = Sr%ze0(jb,sk_ibz,1) * &
&          (SUM(Sr%sigxme(jb,sk_ibz,:)+Sr%sigcmee0(jb,sk_ibz,:)+Sr%hhartree(jb,jb,sk_ibz,:))-Sr%e0(jb,sk_ibz,1))

         Sr%egw(jb,sk_ibz,1)=Sr%e0(jb,sk_ibz,1)+Sr%degw(jb,sk_ibz,1)

         ! Estimate Sigma at the QP-energy.
         do iab=1,Sr%nsig_ab
          Sr%sigmee(jb,sk_ibz,iab)= &
&           Sr%sigxme(jb,sk_ibz,iab)+Sr%sigcmee0(jb,sk_ibz,iab)+Sr%degw(jb,sk_ibz,1)*Sr%dsigmee0(jb,sk_ibz,iab)
         end do
       end if

     end do !jb
   end do !spin

 else
   ! =============================
   ! === Analytic Continuation ===
   ! =============================
   ABI_CHECK(Sr%nsig_ab==1,"AC with spinor not implemented")
   do spin=1,Sr%nsppol
     do jb=ib1,ib2

      ABI_MALLOC(tmpcdp,(Sr%nomega_i))
      ! * Calculate Sigc(E0), dSigc(E0)
      zz=CMPLX(Sr%e0(jb,sk_ibz,spin),zero)

      if (Sigp%mbpt_sciss>0.1d-4) then
       ! RS: e0 is replaced by qp_ene which contains the updated energy eigenvalue
       zz=CMPLX(qp_ene(jb,sk_ibz,spin),zero)
      end if

      ! === Diagonal elements of sigcme_tmp ===
      ! * if zz in 2 or 3 quadrant, avoid poles in the complex plane using Sigma(-iw)=Sigma(iw)*.
      do iab=1,Sr%nsig_ab
        is_idx=spin; if (Sr%nsig_ab>1) is_idx=iab
        if (REAL(zz)>zero) then
          tmpcdp(:)=sigcme_tmp(:,jb,jb,is_idx)
          Sr%sigcmee0(jb,sk_ibz,is_idx)=  pade(Sr%nomega_i,Sr%omega_i,tmpcdp,zz)
          Sr%dsigmee0(jb,sk_ibz,is_idx)= dpade(Sr%nomega_i,Sr%omega_i,tmpcdp,zz)
        else
          tmpcdp(:)=CONJG(sigcme_tmp(:,jb,jb,is_idx))
          Sr%sigcmee0(jb,sk_ibz,is_idx)=  pade(Sr%nomega_i,CONJG(Sr%omega_i),tmpcdp,zz)
          Sr%dsigmee0(jb,sk_ibz,is_idx)= dpade(Sr%nomega_i,CONJG(Sr%omega_i),tmpcdp,zz)
        end if
      end do !iab

      ! Z=(1-dSigma/domega(E0))^-1
      if (Sr%nsig_ab==1) then
        Sr%ze0(jb,sk_ibz,spin) = one/(one-Sr%dsigmee0(jb,sk_ibz,spin))
      else
        Sr%ze0(jb,sk_ibz,1)=one/(one-SUM(Sr%dsigmee0(jb,sk_ibz,:)))
      end if

      ! Find roots of E^0-V_xc-V_U+Sig_x+Sig_c(z)-z, i.e E^qp.
      ! using Newton-Raphson method and starting point E^0
      zz=CMPLX(Sr%e0(jb,sk_ibz,spin),zero)

      if (Sigp%mbpt_sciss>0.1d-4) then ! e0 is replaced by qp_ene which contains the updated energy eigenvalue.
        zz=CMPLX(qp_ene(jb,sk_ibz,spin),0.0)
      end if

      iter=0; converged=.FALSE.; ctdpc=cone
      do while (ABS(ctdpc)>NR_ABS_ROOT_ERR.or.iter<NR_MAX_NITER)
        iter=iter+1
        sigc=czero ; dsigc=czero
        if (REAL(zz)>tol12) then
          tmpcdp(:)=sigcme_tmp(:,jb,jb,spin)
          sigc =  pade(Sr%nomega_i,Sr%omega_i,tmpcdp,zz)
          dsigc= dpade(Sr%nomega_i,Sr%omega_i,tmpcdp,zz)
        else
          tmpcdp(:)=CONJG(sigcme_tmp(:,jb,jb,spin))
          sigc =  pade(Sr%nomega_i,CONJG(Sr%omega_i),tmpcdp,zz)
          dsigc= dpade(Sr%nomega_i,CONJG(Sr%omega_i),tmpcdp,zz)
        end if
        ctdpc = Sr%e0(jb,sk_ibz,spin)-Sr%vxcme(jb,sk_ibz,spin)-Sr%vUme(jb,sk_ibz,spin)+Sr%sigxme(jb,sk_ibz,spin)+sigc-zz
        if (ABS(ctdpc)<NR_ABS_ROOT_ERR) then
         converged=.TRUE.; EXIT
        end if
        dct=dsigc-one
        zz=newrap_step(zz,ctdpc,dct)
      end do

      if (.not.converged) then
        write(msg,'(a,i0,3a,f8.4,a,f8.4)')&
          'Newton-Raphson method not converged after ',NR_MAX_NITER,' iterations. ',ch10,&
          'Absolute Error = ',ABS(ctdpc),' > ',NR_ABS_ROOT_ERR
        MSG_WARNING(msg)
      end if
      !
      ! Store the final result TODO re-shift everything according to efermi
      Sr%egw(jb,sk_ibz,spin)=zz
      Sr%degw(jb,sk_ibz,spin)=Sr%egw(jb,sk_ibz,spin) - Sr%e0(jb,sk_ibz,spin)
      Sr%sigmee(jb,sk_ibz,spin)=Sr%sigxme(jb,sk_ibz,spin) + sigc
      !
      ! Spectra of Sigma, remember that Sr%nomega_r does not contains the frequencies used to evaluate the derivative
      ! each frequency is obtained using the pade_expression
      do io=1,Sr%nomega_r
        zz=Sr%omega_r(io)
        if (REAL(zz)>zero) then
          tmpcdp(:)=sigcme_tmp(:,jb,jb,spin)
          Sr%sigcme(jb,sk_ibz,io,spin) = pade(Sr%nomega_i,Sr%omega_i,tmpcdp,zz)
        else
          tmpcdp(:)=CONJG(sigcme_tmp(:,jb,jb,spin))
          Sr%sigcme(jb,sk_ibz,io,spin) = pade(Sr%nomega_i,CONJG(Sr%omega_i),tmpcdp,zz)
        end if
        Sr%sigxcme(jb,sk_ibz,io,spin)= Sr%sigxme(jb,sk_ibz,spin)+Sr%sigcme(jb,sk_ibz,io,spin)
      end do
      !
      ! === Save sigma values along the imaginary axis ===
      do iab=1,Sr%nsig_ab
        is_idx=spin ; if (Sr%nsig_ab>1) is_idx=iab
        do io=1,Sr%nomega_i
          Sr%sigcmesi (jb,sk_ibz,io,is_idx)= sigcme_tmp(io,jb,jb,is_idx)
          Sr%sigxcmesi(jb,sk_ibz,io,is_idx)= Sr%sigxme(jb,sk_ibz,is_idx)+Sr%sigcmesi(jb,sk_ibz,io,is_idx)
        end do
      end do

      ABI_FREE(tmpcdp)

     end do !jb
   end do !is
 end if ! Analytic continuation.
 !
 ! === Diagonalize the QP Hamiltonian (forced to be Hermitian) ===
 ! * Calculate Sr%en_qp_diago and Sr%eigvec_qp to be written in the QPS file.
 ! TODO in case of AC results are wrong.

 ABI_MALLOC(hhartree,(ib1:ib2,ib1:ib2,Sr%nsppol*Sr%nsig_ab))
 hhartree=Sr%hhartree(ib1:ib2,ib1:ib2,sk_ibz,:)

 ! If non self-consistent erase all off-diagonal elements
 if (Sigp%gwcalctyp<20) then
   do jb=ib1,ib2
     do kb=ib1,ib2
      if (jb==kb) CYCLE
      hhartree(jb,kb,:)=czero
     end do
   end do
 end if

 ABI_MALLOC(htotal,(ib1:ib2,ib1:ib2,Sr%nsppol*Sr%nsig_ab))
 do spin=1,Sr%nsppol*Sr%nsig_ab
   do jb=ib1,ib2
     do kb=ib1,ib2
      htotal(kb,jb,spin) = hhartree(kb,jb,spin) + Sr%x_mat(kb,jb,sk_ibz,spin) + sigcme_tmp(Sr%nomega_r+ioe0j,kb,jb,spin)
     end do
   end do
 end do
 !
 ! === Get the Hermitian part of htotal ===
 ! * In the noncollinear case A_{12}^{ab} = A_{21}^{ba}^* if A is Hermitian.
 ABI_MALLOC(h_tmp1,(ib1:ib2,ib1:ib2))
 ABI_MALLOC(h_tmp2,(ib1:ib2,ib1:ib2))

 nsploop=Sr%nsppol; if (Sr%nsig_ab/=1) nsploop=2
 do spin=1,nsploop
   h_tmp1 = CONJG(htotal(:,:,spin))
   h_tmp2 = TRANSPOSE(h_tmp1)
   h_tmp1 = htotal(:,:,spin)
   htotal(:,:,spin)= half * (h_tmp1 + h_tmp2)
 end do

 ! Print the different matrix elements of sigma if QPSC and prtvol>9
 if (Sigp%gwcalctyp>=20.and.prtvol>9.and.my_rank==master) then
   call print_sigma_melems(ikcalc,ib1,ib2,Sr%nsppol*Sr%nsig_ab,htotal,hhartree,&
&               Sr%x_mat(ib1:ib2,ib1:ib2,sk_ibz,:),sigcme_tmp(Sr%nomega_r+ioe0j,:,:,:),Dtfil%filnam_ds(4))
 end if

 if (Sr%nsig_ab==4) then
   h_tmp1 = CONJG(htotal(:,:,4))
   h_tmp2 = TRANSPOSE(h_tmp1)
   h_tmp1 = htotal(:,:,3)
   htotal(:,:,3)= half * (h_tmp1 + h_tmp2)

   h_tmp1 = CONJG(htotal(:,:,3))
   h_tmp2 = TRANSPOSE(h_tmp1)
   htotal(:,:,4) = h_tmp2
 end if

 ! Solve Herm(htotal)*U = E*U
 ld_matrix=ib2-ib1+1
 ABI_MALLOC(hdp,(ld_matrix,ld_matrix))
 ABI_MALLOC(eig,(ld_matrix))

 do spin=1,Sr%nsppol
   if (Sr%nsig_ab==1) then
     hdp=htotal(ib1:ib2,ib1:ib2,spin)
   else
     hdp=SUM(htotal(ib1:ib2,ib1:ib2,:),DIM=3)
   end if
   if (spin == 3) write(std_out,*) hdp  ! This to work around a compiler bug on tikal_gnu_5.4_mpich

   call xheev("Vectors","Upper",ld_matrix,hdp,eig)

   Sr%eigvec_qp(ib1:ib2,ib1:ib2,sk_ibz,spin)=hdp(:,:)
   Sr%en_qp_diago(ib1:ib2,sk_ibz,spin)=eig(:)
 end do

 ABI_FREE(hdp)
 ABI_FREE(eig)
 ABI_FREE(htotal)
 ABI_FREE(hhartree)
 ABI_FREE(h_tmp1)
 ABI_FREE(h_tmp2)

 call timab(490,2,tsec)

 DBG_EXIT("COLL")

end subroutine solve_dyson
!!***

!----------------------------------------------------------------------

!!****f* m_dyson_solver/print_sigma_melems
!! NAME
!!  print_sigma_melems
!!
!! FUNCTION
!!  This routine prints the Hermitian and the non-hermitian part of the matrix
!!  elements of Sigma, as well as the individual contributions.
!!  The first 14x14 are printed to screen, and the full matrices are printed
!!  to files: sigma_melems_, sigma_nonH_melems_, sigma_Hart_melems_,
!!            sigma_x_melems, and sigma_c_melems
!!
!! INPUTS
!!  ikcalc  : index of k-point
!!  ib1,ib2 : starting and ending band indices
!!  nsp     : no. of spin elements
!!  htotal  : Hermitianised matrix elements of Sigma
!!  hhartree : Hartree contribution to matrix elements
!!  sigxme  : Sigma_x contribution to matrix elements
!!  sigcme  : Sigma_c contribution to matrix elements
!!  prefil : prefix for output files.
!!
!! OUTPUT
!!
!! PARENTS
!!      m_dyson_solver
!!
!! CHILDREN
!!      int2char10,wrtout
!!
!! SOURCE

subroutine print_sigma_melems(ikcalc,ib1,ib2,nsp,htotal,hhartree,sigxme,sigcme,prefil)

! Arguments ------------------------------------
 !scalars
 integer,intent(in) :: ikcalc,ib1,ib2,nsp
 character(len=*),intent(in) :: prefil
 !arrays
 complex(dpc),intent(in) :: htotal(ib1:ib2,ib1:ib2,nsp),hhartree(ib1:ib2,ib1:ib2,nsp)
 complex(dpc),intent(in) :: sigxme(ib1:ib2,ib1:ib2,nsp),sigcme(ib1:ib2,ib1:ib2,nsp)

! Local variables ------------------------------
 integer,parameter :: MAX_NCOLS=14
 integer :: isp,mc,mr,jj,ii,temp_unit,ount
 character(len=10) :: sidx
 character(len=500) :: msg
 character(len=100) :: fmth,fmt1,fmt2,fmthh,kpt_index,fmtfile
 character(len=fnlen) :: filename
! *************************************************************************

 if (nsp==3.or.nsp>4) then
   MSG_ERROR('nsp has wrong value in print_sigma_melems')
 end if

 ount = std_out

 mc = ib2-ib1+1; if (mc>MAX_NCOLS) mc = MAX_NCOLS
 mr = mc

 write(fmthh,*)'(2(a),2(I2,a))'
 write(fmth,*)'(7x,',mc,'(i2,8x))'
 write(fmt1,*)'(3x,i2,',mc,'f10.5)'
 write(fmt2,*)'(5x   ,',mc,'f10.5,a)'

! First print to screen
 do isp=1,nsp
   write(msg,'(a)') ''
   call wrtout(ount,msg,'COLL')
   write(msg,fmthh) ch10,' Hermitianised matrix elements of Sigma (spin ',isp,' of ',nsp,'):'
   call wrtout(ount,msg,'COLL')
   write(msg,fmth)(jj,jj=1,mc)
   call wrtout(ount,msg,'COLL') !header
   do ii=ib1,ib1+mr-1
     write(msg,fmt1)ii-ib1+1,DBLE(htotal(ii,ib1:(ib1+mc-1),isp))
     call wrtout(ount,msg,'COLL') !real part
     write(msg,fmt2)  AIMAG(htotal(ii,ib1:(ib1+mc-1),isp)),ch10
     call wrtout(ount,msg,'COLL') !imag part
   end do
 end do !nsp

 write(msg,'(a,i2,a)')" Max. ",MAX_NCOLS," elements printed. Full matrix output in _HTOTAL files"
 call wrtout(ount,msg,'COLL')

 do isp=1,nsp
   write(msg,fmthh) ch10,' H_Hartree matrix elements (spin ',isp,' of ',nsp,'):'
   call wrtout(ount,msg,'COLL')
   write(msg,fmth)(jj,jj=1,mc)
   call wrtout(ount,msg,'COLL') !header
   do ii=ib1,ib1+mr-1
     write(msg,fmt1)ii-ib1+1,DBLE(hhartree(ii,ib1:(ib1+mc-1),isp))
     call wrtout(ount,msg,'COLL') !real part
     write(msg,fmt2)  AIMAG(hhartree(ii,ib1:(ib1+mc-1),isp)),ch10
     call wrtout(ount,msg,'COLL') !imag part
   end do
 end do !nsp

 write(msg,'(a,i2,a)')" Max. ",MAX_NCOLS," elements printed. Full matrix output in _HHARTREE files"
 call wrtout(ount,msg,'COLL')

 do isp=1,nsp
   write(msg,fmthh) ch10,' Sigma_x matrix elements (spin ',isp,' of ',nsp,'):'
   call wrtout(ount,msg,'COLL')
   write(msg,fmth)(jj,jj=1,mc)
   call wrtout(ount,msg,'COLL') !header
   do ii=ib1,ib1+mr-1
     write(msg,fmt1)ii-ib1+1,DBLE(sigxme(ii,ib1:(ib1+mc-1),isp))
     call wrtout(ount,msg,'COLL') !real part
     write(msg,fmt2)  AIMAG(sigxme(ii,ib1:(ib1+mc-1),isp)),ch10
     call wrtout(ount,msg,'COLL') !imag part
   end do
 end do !nsp

 write(msg,'(a,i2,a)')" Max. ",MAX_NCOLS," elements printed. Full matrix output _SIGX files"
 call wrtout(ount,msg,'COLL')

 do isp=1,nsp
   write(msg,fmthh) ch10,' Sigma_c matrix elements (spin ',isp,' of ',nsp,'):'
   call wrtout(ount,msg,'COLL')
   write(msg,fmth)(jj,jj=1,mc)
   call wrtout(ount,msg,'COLL') !header
   do ii=ib1,ib1+mr-1
     write(msg,fmt1)ii-ib1+1,DBLE(sigcme(ii,ib1:(ib1+mc-1),isp))
     call wrtout(ount,msg,'COLL') !real part
     write(msg,fmt2)  AIMAG(sigcme(ii,ib1:(ib1+mc-1),isp)),ch10
     call wrtout(ount,msg,'COLL') !imag part
   end do
 end do !nsp

 write(msg,'(a,i2,a)')" Max ",MAX_NCOLS," elements printed. Full matrix output _SIGC files"
 call wrtout(ount,msg,'COLL')

 ! Then print to file
 ! Format is: row, column, value; with a blank space for each full
 ! set of columns for easy plotting with the gnuplot splot command
 write(fmtfile,*)'(3X,I6,2X,I6,',nsp,'(2(ES28.16E3,3x)))'

 call int2char10(ikcalc,sidx)
 kpt_index = "_KPT"//TRIM(sidx)

 filename = TRIM(prefil)//'_HTOTAL'//TRIM(kpt_index)

 if (open_file(filename,msg,newunit=temp_unit,form="formatted",status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 msg = '#   row    col.      Re(htotal(r,c)) Im(htotal(r,c))  for spin11   ... spin22 ... spin12 ... spin13'
 call wrtout(temp_unit,msg,'COLL')
 do ii=ib1,ib2
   do jj=ib1,ib2
     write(msg,fmtfile) ii,jj,(htotal(jj,ii,isp),isp=1,nsp)
     call wrtout(temp_unit,msg,'COLL')
   end do
   call wrtout(temp_unit,"",'COLL')
 end do
 close(temp_unit)

 filename = TRIM(prefil)//'_HHARTREE'//TRIM(kpt_index)
 if (open_file(filename,msg,newunit=temp_unit,form="formatted",status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 msg = '#   row    col.      Re(hhartree(r,c))  Im(hhartree(r,c)  for spin11   ... spin22 ... spin12 ... spin13'
 call wrtout(temp_unit,msg,'COLL')
 do ii=ib1,ib2
   do jj=ib1,ib2
     write(msg,fmtfile) ii,jj,(hhartree(jj,ii,isp),isp=1,nsp)
     call wrtout(temp_unit,msg,'COLL')
   end do
   call wrtout(temp_unit,"",'COLL')
 end do
 close(temp_unit)

 filename = TRIM(prefil)//'_SIGX'//TRIM(kpt_index)
 if (open_file(filename,msg,newunit=temp_unit,form="formatted",status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 write(msg,'(a)')'#   row    col.      Re(Sigx(r,c)) Im(Sigx(r,c) for spin11   ... spin22 ... spin12 ... spin13'
 call wrtout(temp_unit,msg,'COLL')
 do ii=ib1,ib2
   do jj=ib1,ib2
     write(msg,fmtfile) ii,jj,(sigxme(jj,ii,isp),isp=1,nsp)
     call wrtout(temp_unit,msg,'COLL')
   end do
   call wrtout(temp_unit,"",'COLL')
 end do
 close(temp_unit)

 filename = TRIM(prefil)//'_SIGC'//TRIM(kpt_index)
 if (open_file(filename,msg,newunit=temp_unit,form="formatted",status="replace",action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 write(msg,'(a)')'#   row    col.      Re(Sigc(r,c)) Im(Sigc(r,c) for spin11   ... spin22 ... spin12 ... spin21'
 call wrtout(temp_unit,msg,'COLL')
 do ii=ib1,ib2
   do jj=ib1,ib2
     write(msg,fmtfile) ii,jj,(sigcme(jj,ii,isp),isp=1,nsp)
     call wrtout(temp_unit,msg,'COLL')
   end do
   call wrtout(temp_unit,"",'COLL')
 end do

 close(temp_unit)

end subroutine print_sigma_melems

!----------------------------------------------------------------------

END MODULE m_dyson_solver
!!***
