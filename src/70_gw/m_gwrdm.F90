!!****m* ABINIT/m_gwrdm
!! NAME
!!  m_gwrdm
!!
!! FUNCTION
!!  Compute density matrix correction Galitskii-Migdal Ecorr, G = Go + Go Sigma Go (imaginary freqs. are used in Sigma_c) 
!!  and associated quantities (natural orbitals, matrix elements, etc.).
!! PARENTS
!! 
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gwrdm

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_hide_blas
 use m_time
 use m_wfd           
 use m_hdr
 use m_bz_mesh,       only : kmesh_t, kmesh_free, littlegroup_t, littlegroup_init, littlegroup_free, &
                             kmesh_init, has_BZ_item, isamek, get_ng0sh, kmesh_print, &
                             get_bz_item, has_IBZ_item, find_qmesh
 use m_dtset

 use defs_datatypes,   only : ebands_t
 use m_sigma,          only : sigma_t
 use m_xctk,           only : xcden  
 implicit none

 private :: no2ks,ks2no,printdm1 
!!***
 
 public :: calc_Ec_GM_k,calc_rdmx,calc_rdmc,natoccs,update_hdr_bst,rotate_ks_no,print_tot_occ 
!!***

contains

!!***
!!****f* ABINIT/Calc_Ec_GM_k
!! NAME
!! calc_Ec_GM_k
!!
!! FUNCTION
!! Calculate Galitskii-Migdal corr. energy integrated in the Imaginary axis Ec = 1/pi \sum_i \int Gii(iv)*Sigma_c,ii(iv) + cc. dv
!!
!! INPUTS
!! ib1=min band for given k
!! ib2=max band for given k.
!! ik_ibz= the label of k-point in the IBZ whose Galitskii-Migdal contribution is accounted.
!! weights=array containing the weights used in the quadrature.
!! sigcme_k=array containing Sigma(iw) as Sigma(iw,ib1:ib2,ib1:ib2,nspin)
!! dm1=density matrix, matrix (i,j), where i and j belong to the k-point k (see m_sigma_driver.F90 for more details).
!! Bst=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!! Sr=sigma_t (see the definition of this structured datatype)
!!
!! OUTPUT
!! Compute the Galitskii-Migdal corr energy contribution of this k-point:
!! Ec ^k = 1/(4*pi) * fact_spin * int _{ -Inf }^{ +Inf } dv Sigma_c ^k (iv) * G0(iv) 
!!       = 1/(4*pi) * fact_spin * int _{   0  }^{ +Inf } dv 2 * Re{ Sigma_c ^k (iv) * G0(iv) } 
!!
!! PARENTS
!!  m_sigma_driver.f90
!! CHILDREN
!! SOURCE

function calc_Ec_GM_k(ib1,ib2,ik_ibz,Sr,weights,sigcme_k,BSt) result(Ec_GM_k)

!Arguments ------------------------------------
!scalars
 real(dp) :: Ec_GM_k
 integer,intent(in) :: ib1,ib2,ik_ibz
 type(ebands_t),target,intent(in) :: BSt
 type(sigma_t),intent(in) :: Sr
!arrays
 real(dp),intent(in) :: weights(:)
 complex(dpc),intent(in) :: sigcme_k(:,:,:,:)
!Local variables ------------------------------
!scalars
 integer :: ibdm!,unitt
 real(dp) :: ec_integrated,spin_fact,fact
!arrays
!************************************************************************

 DBG_ENTER("COLL")

 ec_integrated=0.0_dp
 spin_fact=2.0_dp
 fact=spin_fact*(1.0_dp/(two_pi*2.0_dp))

 if (ib1/=1) then
   MSG_WARNING("Unable to compute the Galitskii-Migdal correlation energy because the first band was &
   &not included in bdgw interval. Restart the calculation starting bdgw from 1.")
 else
   ! WARNING: Sigma_c(iv) produced from a previous integration at the screening stage, is numerically not much stable and introduces bumps.
   ! Unfortunately, the Green's function times Sigma_c(iv) does not decay fast enough with iv to overcome the bumps. These bumps are
   ! not pronouced for the linearized density matrix update, as two Green's functions are multiplied making the decay much faster with iv. 
   ! If a better way to produce more stable Sigma_c(iv) values is found, this subroutine can be use to evaluate GM Ecorr in the future. TODO
   do ibdm=1,ib2
     ! Sigma_pp(iv)/[(iv - e_ibdm,k)] + [Sigma_pp(iv)/[(iv - e_ibdm,k)]]^* = 2 Re [Sigma_pp(iv)/(iv - e_ibdm,k)]
     ec_integrated=ec_integrated+2.0_dp*real( sum(weights(:)*sigcme_k(:,ibdm,ibdm,1)/(Sr%omega_i(:)-BSt%eig(ibdm,ik_ibz,1)) ) )
   end do
 endif
 
 Ec_GM_k=fact*ec_integrated

 DBG_EXIT("COLL")

end function calc_Ec_GM_k
!!***

!!****f* ABINIT/calc_rdmx
!! NAME
!! calc_rdmx
!!
!! FUNCTION
!! Calculate density matrix corrections for G = Go + Go (Sigma_x - alpha*Sigma_x - Vxc) Go
!!
!! INPUTS
!! ib1=min band for given k
!! ib2=max band for given k.
!! ik_ibz= the label of k-point in the IBZ.
!! dm1=density matrix, matrix (i,j), where i and j belong to the k-point k (see m_sigma_driver.F90 for more details). 
!! pot=Self-energy-Potential difference, matrix size (i,j), where i and j belong to k.
!! Bst=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!!
!! OUTPUT
!! Updated dm1 matrix array with Go (Sigma_x - alpha*Sigma_x - Vxc) Go
!! PARENTS
!!  m_sigma_driver.f90
!! CHILDREN
!! SOURCE

subroutine calc_rdmx(ib1,ib2,ik_ibz,pot,dm1,BSt)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2,ik_ibz
 type(ebands_t),target,intent(in) :: BSt
!arrays
 complex(dpc),intent(in) :: pot(:,:)
 complex(dpc),intent(inout) :: dm1(:,:)
!Local variables ------------------------------
!scalars
 character(len=500) :: msg
 integer :: ib1dm,ib2dm
 real(dp) :: spin_fact,tol8
!arrays
!************************************************************************

 DBG_ENTER("COLL")
 tol8=1.0e-8
 spin_fact=2.0_dp

 write(msg,'(a58,3f10.5)')' Computing the 1-RDM correction for  Sx-Vxc  and k-point: ',BSt%kptns(1:,ik_ibz)
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 write(msg,'(a11,i5,a8,i5)')'from band ',ib1,' to band',ib2
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 do ib1dm=ib1,ib2-1  
   do ib2dm=ib1dm+1,ib2
     if ((BSt%occ(ib1dm,ik_ibz,1)>tol8) .and. (BSt%occ(ib2dm,ik_ibz,1)<tol8)) then
       dm1(ib1dm,ib2dm)=spin_fact*pot(ib1dm,ib2dm)/(BSt%eig(ib1dm,ik_ibz,1)-BSt%eig(ib2dm,ik_ibz,1)+tol8)
       ! Dji = Dij^*
       dm1(ib2dm,ib1dm)=conjg(dm1(ib1dm,ib2dm))
     end if
   end do
 end do
 
 DBG_EXIT("COLL")

end subroutine calc_rdmx
!!***

!!****f* ABINIT/calc_rdmc
!! NAME
!! calc_rdmc
!!
!! FUNCTION
!! Calculate density matrix corrections for G = Go + int Go(iw) Sigma_c(iw) Go(iw) dw
!!
!! INPUTS
!! ib1=min band for given k
!! ib2=max band for given k.
!! ik_ibz= the label of k-point in the IBZ.
!! weights=array containing the weights used in the quadrature.
!! sigcme_k=array containing Sigma(iw) as Sigma(iw,ib1:ib2,ib1:ib2,nspin) 
!! dm1=density matrix, matrix (i,j), where i and j belong to the k-point k (see m_sigma_driver.F90 for more details). 
!! Bst=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!! Sr=sigma_t (see the definition of this structured datatype)
!!
!! OUTPUT
!! Updated dm1 matrix array with int Go(iw) Sigma_c(iw) Go(iw) dw
!! PARENTS
!!  m_sigma_driver.f90
!! CHILDREN
!! SOURCE

subroutine calc_rdmc(ib1,ib2,ik_ibz,Sr,weights,sigcme_k,BSt,dm1)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2,ik_ibz
 type(ebands_t),target,intent(in) :: BSt
 type(sigma_t) :: Sr
!arrays
 real(dp),intent(in) :: weights(:)
 complex(dpc),intent(inout) :: dm1(:,:)
 complex(dpc),intent(in) :: sigcme_k(:,:,:,:)
!Local variables ------------------------------
!scalars
 real(dp) :: spin_fact,fact
 character(len=500) :: msg
 integer :: ib1dm,ib2dm 
!arrays
!************************************************************************

 DBG_ENTER("COLL")

 spin_fact=2.0_dp
 fact=spin_fact*(1.0_dp/two_pi)

 write(msg,'(a58,3f10.5)')' Computing the 1-RDM correction for  Sc(iw)  and k-point: ',BSt%kptns(1:,ik_ibz)
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 write(msg,'(a11,i5,a8,i5)')'from band ',ib1,' to band',ib2
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 dm1(:,:)=czero
 do ib1dm=ib1,ib2  
   do ib2dm=ib1dm,ib2 
     ! Sigma_pq/[(denominator)] + [Sigma_qp/[(denominator)]]^*
     dm1(ib1dm,ib2dm)=fact*sum(weights(:)*( sigcme_k(:,ib1dm,ib2dm,1)/&
                 &( (Sr%omega_i(:)-BSt%eig(ib1dm,ik_ibz,1))*(Sr%omega_i(:)-BSt%eig(ib2dm,ik_ibz,1)) )&
                                    +conjg( sigcme_k(:,ib2dm,ib1dm,1)/& 
                 &( (Sr%omega_i(:)-BSt%eig(ib1dm,ik_ibz,1))*(Sr%omega_i(:)-BSt%eig(ib2dm,ik_ibz,1)) ) ) ) ) 
     ! Dji = Dij^*
     dm1(ib2dm,ib1dm)=conjg(dm1(ib1dm,ib2dm))
   end do  
 end do

 DBG_EXIT("COLL")

end subroutine calc_rdmc
!!***

!!****f* ABINIT/natoccs
!! NAME
!! natoccs
!!
!! FUNCTION
!! Calculate natural orbitals and occ. numbers for a given k-point 
!!
!! INPUTS
!! ib1=min band for given k
!! ib2=max band for given k.
!! ik_ibz= the label of k-point in the IBZ.
!! iinfo=use Sigma_x or Sigma_c phaser
!! weights=array containing the weights used in the quadrature.
!! nateigv=array containing the natural eigenvectors in columns (nbands,nband,k-point,nspin)
!! dm1=density matrix, matrix (i,j), where i and j belong to the k-point k (see m_sigma_driver.F90 for more details).
!! occs = array containing the occ numbers for a given k-point occs(nband,k-point).
!! Bst=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!! checksij=check the orthonormality of the nat. orbitals
!!
!! OUTPUT
!! Compute the nat. orbitals and occ. numbers from the dm1 matrix (for exchange and correlations)
!! PARENTS
!!  m_sigma_driver.f90
!! CHILDREN
!! SOURCE

subroutine natoccs(ib1,ib2,dm1,nateigv,occs,BSt,ik_ibz,iinfo,checksij)
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2,ik_ibz,iinfo
 integer,intent(in),optional::checksij
 type(ebands_t),target,intent(in) :: BSt
!arrays
 real(dp),intent(inout) :: occs(:,:)
 complex(dpc),intent(inout) :: dm1(:,:),nateigv(:,:,:,:)
!Local variables ------------------------------
!scalars
 integer::ndim,ib1dm,ib2dm,ib3dm,lwork,info
 logical::check_Sijmat
 character(len=500) :: msg
 real(dp) :: toccs_k,tol10
 complex(dp) :: Sib1k_ib2k
!arrays
 real(dp),allocatable :: occs_tmp(:),occs_tmp2(:),rwork(:)
 complex(dpc),allocatable :: work(:),dm1_tmp(:,:),eigenvect(:,:)
!************************************************************************

 DBG_ENTER("COLL")
 
 check_Sijmat=.false.
 if (present(checksij)) then
  check_Sijmat=.true.
 end if
 tol10=1.0d-10

 ndim=ib2-ib1+1
 lwork=2*ndim-1
 ABI_MALLOC(occs_tmp,(ndim))
 ABI_MALLOC(occs_tmp2,(ndim))
 ABI_MALLOC(work,(lwork))
 ABI_MALLOC(dm1_tmp,(ndim,ndim))
 ABI_MALLOC(eigenvect,(ndim,ndim))
 ABI_MALLOC(rwork,(3*ndim-2))

 dm1_tmp=0.0_dp
 do ib2dm=1,ndim
   do ib1dm=ib2dm,ndim
     dm1_tmp(ib1dm,ib2dm)=dm1(ib1+(ib1dm-1),ib1+(ib2dm-1))
     ! Dji = Dij^*
     dm1_tmp(ib2dm,ib1dm)=conjg(dm1_tmp(ib1dm,ib2dm))
   end do
 end do

 work=0.0_dp
 occs_tmp=0.0_dp
 info=0
 call zheev('v','u',ndim,dm1_tmp,ndim,occs_tmp,work,lwork,rwork,info)
 if (info/=0) then
   MSG_WARNING("Failed the diagonalization of the updated GW 1-RDM")
 end if

 ! Uncomment for debug 
 !write(msg,'(a6)') 'Eigvec'
 !call wrtout(std_out,msg,'COLL')
 !call printdm1(1,10,dm1_tmp) 
 !eigenvect=dm1_tmp
 !occs_tmp2=occs_tmp
 !Order from highest occ to lowest occ
 do ib1dm=1,ndim
  occs_tmp2(ib1dm)=occs_tmp(ndim-(ib1dm-1))
  do ib2dm=1,ndim
   eigenvect(ib2dm,ib1dm)=dm1_tmp(ib2dm,(ndim-(ib1dm-1)))
  end do
  if (abs(occs_tmp2(ib1dm))<tol10) then
    occs_tmp2(ib1dm)=0.0_dp
  end if
 end do

 if (check_Sijmat) then 
   do ib1dm=1,ndim
     do ib2dm=1,ib1dm
       Sib1k_ib2k=czero
       do ib3dm=1,ndim
         Sib1k_ib2k=Sib1k_ib2k+conjg(eigenvect(ib3dm,ib1dm))*eigenvect(ib3dm,ib2dm)
       end do
       if (ib1dm==ib2dm) then
         if(abs(Sib1k_ib2k-cmplx(1.0_dp,0.0_dp))>tol10) then
           write(msg,'(a45,i5,a1,i5,f10.5)') 'Large deviation from identity for bands ',ib1dm,' ',ib2dm,real(Sib1k_ib2k) 
           call wrtout(std_out,msg,'COLL')
         endif
       else
         if (abs(Sib1k_ib2k)>tol10) then
           write(msg,'(a45,i5,a1,i5,f10.5)') 'Large deviation from identity for bands ',ib1dm,' ',ib2dm,real(Sib1k_ib2k) 
           call wrtout(std_out,msg,'COLL')
         end if
       end if
     end do
   end do
 end if

 if (info==0) then
   if (iinfo==0) then       
     write(msg,'(a51,3f10.5)') 'Occs. after updating with Sx-Vxc corr. at k-point:',BSt%kptns(1:,ik_ibz)
   else
     write(msg,'(a51,3f10.5)') 'Occs. after updating with S_c correct. at k-point:',BSt%kptns(1:,ik_ibz)
   endif 
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   ib1dm=ndim-(ndim/10)*10
   do ib2dm=1,(ndim/10)*10,10
     write(msg,'(f11.5,9f10.5)') occs_tmp2(ib2dm:ib2dm+9)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   end do  
   ib1dm=(ndim/10)*10+1
   write(msg,'(f11.5,*(f10.5))') occs_tmp2(ib1dm:)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 else
   write(msg,'(a36,3f10.5)') 'Error computing occs. for k-point: ',BSt%kptns(1:,ik_ibz)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if

 ! Store natural orbital eigenvectors matrix and occs. Also compute total number of electrons for this k-point
 toccs_k=0.0_dp
 do ib1dm=1,ndim
   do ib2dm=1,ndim
     nateigv(ib1+(ib1dm-1),ib1+(ib2dm-1),ik_ibz,1)=eigenvect(ib1dm,ib2dm)
   end do
   occs(ib1+(ib1dm-1),ik_ibz)=occs_tmp2(ib1dm)  ! Overwrite the initial KS-DFT occs from ib1 to ib2
   toccs_k=toccs_k+occs_tmp2(ib1dm)
 end do

 write(msg,'(a22,i5,a3,i5,a21,f10.5)') ' Total occ. from band ',ib1,' to', ib2,' at current k-point: ',toccs_k
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 write(msg,'(a5)') ' '
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 ABI_FREE(rwork)
 ABI_FREE(work)
 ABI_FREE(dm1_tmp)
 ABI_FREE(eigenvect)
 ABI_FREE(occs_tmp)
 ABI_FREE(occs_tmp2)

 DBG_EXIT("COLL")

end subroutine natoccs
!!***

!!****f* ABINIT/update_hdr_bst
!! NAME
!! update_hdr_bst
!!
!! FUNCTION
!! Update the Hdr for the WFK and DEN files and the occ. numbers in the BSt file for a given k-point
!!
!! INPUTS
!! Wfd<wfd_t>=Datatype gathering data on QP amplitudes.
!! ngfft_in(18)=information on the fine FFT grid used for densities and potentials.
!! b1gw=min band for given k in the interval where we update.
!! b2gw=max band for given k in the interval where we update.
!! occs= array containing the occ numbers for a given k-point occs_ks(nband,k-point).
!! Bst=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!!
!! OUTPUT
!! Updated Hdr and BSt information
!! PARENTS
!!  m_sigma_driver.f90
!! CHILDREN
!! SOURCE

subroutine update_hdr_bst(Wfd,occs,b1gw,b2gw,BSt,Hdr,ngfft_in)
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: b1gw,b2gw
 integer,intent(in),dimension(3) :: ngfft_in
 type(ebands_t),target,intent(inout) :: BSt
 type(Hdr_type),intent(inout) :: Hdr
 type(wfd_t),intent(in) :: Wfd
!arrays
 real(dp),intent(in) :: occs(:,:)
!Local variables ------------------------------
!scalars
 integer :: ib1dm,ib2dm,dim_bands,ikpoint
!arrays
!************************************************************************
 DBG_ENTER("COLL")

 ! BSt occ (QP_BSt ones) are changed and never recoverd
 do ikpoint=1,BSt%nkpt
   BSt%occ(b1gw:b2gw,ikpoint,1) = occs(b1gw:b2gw,ikpoint) ! Spins summed, occ in [0:2] 
 enddo
 MSG_COMMENT("QP_BSt: occupancies were updated with nat. orb. ones")
 if ((size(Hdr%occ(:))/BSt%nkpt) < (b2gw-b1gw+1)) then
   !Actually, we should never reach this point because the code should stop during Wfd initialization in m_sigma_driver.F90
   MSG_ERROR("Impossible to use the existing read WFK to build a new one!")
 end if
 
 ! Update occ in Hdr before printing
 ib1dm=1
 do ikpoint=1,BSt%nkpt
   dim_bands=size(BSt%occ(:,ikpoint,1))
   do ib2dm=1,dim_bands   
     Hdr%occ(ib1dm)=BSt%occ(ib2dm,ikpoint,1) ! Because Hdr%occ is a 1-D array
     ib1dm=ib1dm+1
   end do
 end do

 Hdr%npwarr(:)=Wfd%npwarr(:)                                   ! Use the npw and ngfft = ones used in GW calc
 Hdr%ngfft(1:3)=ngfft_in(1:3)
 MSG_COMMENT("Hdr_sigma: occupancies, npw, and ngfft were updated")

end subroutine update_hdr_bst
!!***

!!****f* ABINIT/print_tot_occ
!! NAME
!! print_tot_occ
!!
!! FUNCTION
!! Compute and print the total (averaged) occ. from all k-points
!!
!! INPUTS
!! Kmesh <kmesh_t>=Structure describing the k-point sampling.
!! Sr=sigma_t (see the definition of this structured datatype)
!! BSt=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!!
!! OUTPUT
!! Print the total (averaged) occ. = sum_k weight_k * Nelec_k 
!!
!! PARENTS
!!  m_sigma_driver.f90
!! CHILDREN
!! SOURCE

subroutine print_tot_occ(sigma,kmesh,BSt)

!Arguments ------------------------------------
!scalars
 type(sigma_t),intent(in) :: sigma
 type(kmesh_t),intent(in) :: kmesh
 type(ebands_t),intent(in) :: BSt
!Local variables-------------------------------
!scalars
 character(len=500) :: msg
 integer :: ik,ib,spin
 real(dp) :: wtk,occ_bks,tot_occ

! *************************************************************************

 tot_occ=zero

 do spin=1,sigma%nsppol
   do ik=1,sigma%nkibz
     wtk = kmesh%wt(ik)
     do ib=sigma%b1gw,sigma%b2gw
       occ_bks = BSt%occ(ib,ik,spin)
       if (sigma%nsig_ab==1) then ! Only closed-shell restricted is programed
         tot_occ=tot_occ+occ_bks*wtk
       end if
     end do
   end do
 end do

 write(msg,'(a1)') ' '
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 write(msg,'(a39,f10.5)') ' Total averaged occ. from all k-points: ',tot_occ
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')
 write(msg,'(a1)') ' '
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

end subroutine print_tot_occ
!!***

!!****f* ABINIT/rotate_ks_no
!! NAME
!! rotate_ks_no
!!
!! FUNCTION
!! Rotate a matrix from KS to NO basis and vicerversa.
!!
!! INPUTS
!! ib1=min band for given k
!! ib2=max band for given k.
!! Umat=array containing the eigenvectors in Columns (a unitary matrix)
!! Mat=initially an array containing the matrix elements in KS or NO basis
!! option=0 rotate from NO -> KS | 1 rotate from KS -> NO
!!
!! OUTPUT
!! Rotate a matrix from KS to NO basis and vicerversa and save the new matrix on Mat.
!! Mat=at the end an array containing the matrix elements in NO or KS basis
!!
!! PARENTS
!!  m_sigma_driver.f90
!! CHILDREN
!! SOURCE

subroutine rotate_ks_no(ib1,ib2,Mat,Umat,option)
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2,option
!arrays
 complex(dpc),intent(in) :: Umat(:,:)
 complex(dpc),intent(inout) :: Mat(:,:)
!Local variables ------------------------------
!scalars
 integer:: ndim
!arrays
!************************************************************************
 ndim=ib2-ib1+1
 if (option==0) then
   call no2ks(ndim,Mat,Umat)
 else
   call ks2no(ndim,Mat,Umat)
 end if
end subroutine rotate_ks_no
!!***

!!****f* ABINIT/ks2no
!! NAME
!! ks2no
!!
!! FUNCTION
!! Transform the matrix mat from KS to NO basis
!!
!! INPUTS
!! dim=dimension of the matrices 
!! mat=array in the KS basis
!! rot=unitary matrix containg the eigenvectors in NO basis
!!
!! OUTPUT
!! mat=array in the NO basis 
!!
!! PARENTS
!!  rotate_ks_no
!! CHILDREN
!! SOURCE

subroutine ks2no(ndim,mat,rot) 
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndim 
!arrays
 complex(dpc),dimension(:,:),intent(in) :: rot
 complex(dpc),dimension(:,:),intent(inout) :: mat
!Local variables ------------------------------
!scalars
!arrays
 complex(dpc),allocatable :: res(:,:)
!************************************************************************
 ABI_MALLOC(res,(ndim,ndim))
 res=czero

 ! <NO|Op|NO> =  (U^t)* <KS|Op|KS> U
 res=matmul(conjg(transpose(rot)),mat)
 mat=matmul(res,rot)

 ABI_FREE(res)

end subroutine ks2no
!!***

!!****f* ABINIT/no2ks
!! NAME
!! no2ks
!!
!! FUNCTION
!! Transform the matrix mat from NO to KS basis
!!
!! INPUTS
!! dim=dimension of the matrices
!! mat=array in the KS basis
!! rot=unitary matrix containg the eigenvectors in NO basis
!!
!! OUTPUT
!! mat=array in the KS basis
!!
!! PARENTS
!!  rotate_ks_no
!! CHILDREN
!! SOURCE

subroutine no2ks(ndim,mat,rot) 
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndim 
!arrays
 complex(dpc),dimension(:,:),intent(in) :: rot
 complex(dpc),dimension(:,:),intent(inout) :: mat
!Local variables ------------------------------
!scalars
!arrays
 complex(dpc),allocatable :: res(:,:)
!************************************************************************
 ABI_MALLOC(res,(ndim,ndim))
 res=czero

 ! <KS|Op|KS> = U <NO|Op|NO> (U^t)*
 res=matmul(rot,mat)
 mat=matmul(res,conjg(transpose(rot)))

 ABI_FREE(res)

end subroutine no2ks
!!***

!!****f* ABINIT/printdm1
!! NAME
!! printdm1
!!
!! FUNCTION
!! Print the DM1 matrix
!!
!! INPUTS
!! ib1=min band.
!! ib2=max band.
!! dm1=array containing the 1-RDM matrix
!!
!! OUTPUT
!! Print the 1-RDM matrix
!! PARENTS
!! CHILDREN
!! SOURCE

subroutine printdm1(ib1,ib2,dm1) ! Only used for debug on this file, do not use it with large arrays!
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2
!arrays
 complex(dpc),intent(in) :: dm1(:,:)
!Local variables ------------------------------
!scalars
 integer::ib1dm
 character(len=500) :: msg
!arrays
!************************************************************************
 do ib1dm=ib1,ib2
   write(msg,'(a2,*(f10.5))') '  ',Real(dm1(ib1dm,ib1:ib2))
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end do
end subroutine printdm1
!!***

end module m_gwrdm
!!***
