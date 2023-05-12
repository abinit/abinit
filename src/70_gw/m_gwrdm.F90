!!****m* ABINIT/m_gwrdm
!! NAME
!!  m_gwrdm
!!
!! FUNCTION
!!  Compute density matrix correction Galitskii-Migdal Ecorr, G = Go + Go Sigma Go (imaginary freqs. are used in Sigma_c)
!!  and associated quantities (natural orbitals, matrix elements, etc.).
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gwrdm

 use defs_basis
 use m_gwdefs
 use m_abicore
 use m_xmpi
 use m_errors
 use m_hide_blas
 use m_time
 use m_wfd
 use m_hdr
 use m_dtset

 use m_fstrings,       only : sjoin, itoa
 use m_melemts,        only : melements_t
 use m_bz_mesh,        only : kmesh_t
 use defs_datatypes,   only : ebands_t
 use m_sigma,          only : sigma_t
 use m_xctk,           only : xcden
 use m_gaussian_quadrature, only: cgqf

 implicit none

 private :: no2ks,ks2no,printrdm_k,rotate_ks_no
!!***

 public :: quadrature_sigma_cw,calc_Ec_GM_k,calc_rdmx,calc_rdmc,natoccs,update_hdr_bst,print_tot_occ,change_matrix
 public :: print_chkprdm,print_total_energy,print_band_energies,get_chkprdm
!!***

contains
!!***

!!****f* ABINIT/quadrature_sigma_cw
!! NAME
!! quadrature_sigma_cw
!!
!! FUNCTION
!!  Quadrature frequencies used for Sigma_c(iw) integration
!!
!! INPUTS
!! Sigp<sigparams_t>=Parameters governing the self-energy calculation.
!! Sr=sigma_t (see the definition of this structured datatype)
!! weights=real quadrature weights.
!!
!! OUTPUT
!! Update Sigp and Sr imaginary frequencies with iw, and weights with the quadrature weights
!!
!! SOURCE

subroutine quadrature_sigma_cw(Sigp,Sr,weights)
!Arguments ------------------------------------
!scalars
 type(sigparams_t),intent(inout) :: Sigp
 type(sigma_t),intent(inout) :: Sr
!arrays
 real(dp),intent(inout) :: weights(:)

!Local variables ------------------------------
!scalars
 integer :: ifreqs,order_int,gaussian_kind,units(2)
 real(dp) :: gwalpha,gwbeta,wmin,wmax
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: freqs(:)
!************************************************************************

 units = [std_out, ab_out]

 order_int=Sigp%nomegasi
 write(msg,'(a45,i9)')' number of imaginary frequencies for Sigma_c ',order_int
 call wrtout(units, msg)
 write(msg,'(a1)')' '
 call wrtout(units, msg)
 order_int=Sigp%nomegasi
 ABI_MALLOC(freqs,(order_int))
 gaussian_kind=1
 gwalpha=zero
 gwbeta=zero
 wmin=zero
 wmax=one
 call cgqf(order_int,gaussian_kind,gwalpha,gwbeta,wmin,wmax,freqs,weights)
 ! From  0 to 1 -> 0 to infinity
 weights(:)=weights(:)/(one-freqs(:))**two
 freqs(:)=freqs(:)/(one-freqs(:))
 ! Form complex frequencies from 0 to iInf and print them in the log file
 write(msg,'(a52)')'           Re(iw)           Im(iw)           Weight  '
 call wrtout(std_out,msg)
 write(msg,'(a52)')'          --------         --------         -------- '
 call wrtout(std_out,msg)
 do ifreqs=1,order_int
   Sigp%omegasi(ifreqs)=cmplx(zero,freqs(ifreqs))
   Sr%omega_i(ifreqs)=Sigp%omegasi(ifreqs)
   write(msg,'(3f17.5)') Sr%omega_i(ifreqs),weights(ifreqs)
   call wrtout(std_out,msg)
 enddo
 ABI_FREE(freqs)

end subroutine quadrature_sigma_cw
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
!! rdm_k=density matrix, matrix (i,j), where i and j belong to the k-point k (see m_sigma_driver for more details).
!! ebands=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!! Sr=sigma_t (see the definition of this structured datatype)
!!
!! OUTPUT
!! Compute the Galitskii-Migdal corr energy contribution of this k-point:
!! Ec ^k = 1/(4*pi) * fact_spin * int _{ -Inf }^{ +Inf } dv Sigma_c ^k (iv) * G0(iv)
!!       = 1/(4*pi) * fact_spin * int _{   0  }^{ +Inf } dv 2 * Re{ Sigma_c ^k (iv) * G0(iv) }
!!
!! SOURCE

function calc_Ec_GM_k(ib1,ib2,ik_ibz,Sr,weights,sigcme_k,ebands) result(Ec_GM_k)

!Arguments ------------------------------------
!scalars
 real(dp) :: Ec_GM_k
 integer,intent(in) :: ib1,ib2,ik_ibz
 type(ebands_t),target,intent(in) :: ebands
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

 ec_integrated=zero
 spin_fact=two
 fact=spin_fact*(one/(two_pi*two))

 if (ib1/=1) then
   ABI_WARNING("Unable to compute the Galitskii-Migdal correlation energy because the first band was &
   &not included in bdgw interval. Restart the calculation starting bdgw from 1.")
 else
   ! WARNING: Sigma_c(iv) produced from a previous integration at the screening stage, is numerically not much stable and introduces bumps.
   ! Unfortunately, the Green's function times Sigma_c(iv) does not decay fast enough with iv to overcome the bumps. These bumps are
   ! not pronouced for the linearized density matrix update, as two Green's functions are multiplied making the decay much faster with iv.
   ! If a better way to produce more stable Sigma_c(iv) values is found, this subroutine can be use to evaluate GM Ecorr in the future. TODO
   do ibdm=1,ib2
     ! Sigma_pp(iv)/[(iv - e_ibdm,k)] + [Sigma_pp(iv)/[(iv - e_ibdm,k)]]^* = 2 Re [Sigma_pp(iv)/(iv - e_ibdm,k)]
     ec_integrated=ec_integrated+two*real( sum(weights(:)*sigcme_k(:,ibdm,ibdm,1)/(Sr%omega_i(:)-ebands%eig(ibdm,ik_ibz,1)) ) )
   end do
 endif

 Ec_GM_k=fact*ec_integrated

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
!! rdm_k=density matrix, matrix (i,j), where i and j belong to the k-point k (see m_sigma_driver for more details).
!! pot=Self-energy-Potential difference, matrix size (i,j), where i and j belong to k.
!! ebands=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!!
!! OUTPUT
!! Updated rdm_k matrix array with Go (Sigma_x - alpha*Sigma_x - Vxc) Go
!!
!! SOURCE

subroutine calc_rdmx(ib1,ib2,ik_ibz,pot,rdm_k,ebands)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2,ik_ibz
 type(ebands_t),target,intent(in) :: ebands
!arrays
 complex(dpc),intent(in) :: pot(:,:)
 complex(dpc),intent(inout) :: rdm_k(:,:)

!Local variables ------------------------------
!scalars
 character(len=500) :: msg
 integer :: ib1dm,ib2dm, units(2)
 real(dp) :: spin_fact,tol8
!************************************************************************

 tol8=1.0e-8
 spin_fact=two
 units = [std_out, ab_out]

 write(msg,'(a58,3f10.5)')' Computing the 1-RDM correction for  Sx-Vxc  and k-point: ',ebands%kptns(:,ik_ibz)
 call wrtout(units, msg)
 write(msg,'(a11,i5,a8,i5)')'from band ',ib1,' to band',ib2
 call wrtout(units, msg)

 rdm_k(:,:)=czero
 do ib1dm=ib1,ib2-1
   do ib2dm=ib1dm+1,ib2
     if ((ebands%occ(ib1dm,ik_ibz,1)>tol8) .and. (ebands%occ(ib2dm,ik_ibz,1)<tol8)) then
       rdm_k(1+(ib1dm-ib1),1+(ib2dm-ib1))=spin_fact&
            &*pot(1+(ib1dm-ib1),1+(ib2dm-ib1))/(ebands%eig(ib1dm,ik_ibz,1)-ebands%eig(ib2dm,ik_ibz,1)+tol8)
       ! Dji = Dij^*
       rdm_k(1+(ib2dm-ib1),1+(ib1dm-ib1))=conjg(rdm_k(1+(ib1dm-ib1),1+(ib2dm-ib1)))
     end if
   end do
 end do

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
!! omega_i=Frequencies along the imaginary axis.
!! weights=array containing the weights used in the quadrature.
!! sigcme_k=array containing Sigma(iw) as Sigma(iw,ib1:ib2,ib1:ib2,nspin)
!! rdm_k=density matrix, matrix (i,j), where i and j belong to the k-point k (see m_sigma_driver for more details).
!! ebands=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!
!! OUTPUT
!! Updated rdm_k matrix array with int Go(iw) Sigma_c(iw) Go(iw) dw
!!
!! SOURCE

subroutine calc_rdmc(ib1,ib2,ik_ibz,omega_i,weights,sigcme_k,ebands,rdm_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2,ik_ibz
 type(ebands_t),target,intent(in) :: ebands
 complex(dpc),intent(in) :: omega_i(:)
!arrays
 real(dp),intent(in) :: weights(:)
 complex(dpc),intent(inout) :: rdm_k(:,:)
 complex(dpc),intent(in) :: sigcme_k(:,:,:,:)

!Local variables ------------------------------
!scalars
 real(dp) :: spin_fact,fact
 integer :: ib1dm, ib2dm, units(2)
 character(len=500) :: msg
!************************************************************************

 spin_fact = two
 fact = spin_fact * (one/two_pi)
 units = [std_out, ab_out]

 write(msg,'(a58,3f10.5)')' Computing the 1-RDM correction for  Sc(iw)  and k-point: ',ebands%kptns(1:,ik_ibz)
 call wrtout(units , msg)
 write(msg,'(a11,i5,a8,i5)')'from band ',ib1,' to band',ib2
 call wrtout(units , msg)

 rdm_k(:,:)=czero
 do ib1dm=ib1,ib2
   do ib2dm=ib1dm,ib2
     ! Sigma_pq/[(denominator)] + [Sigma_qp/[(denominator)]]^*
     rdm_k(1+(ib1dm-ib1),1+(ib2dm-ib1))=fact*sum(weights(:)*( sigcme_k(:,1+(ib1dm-ib1),1+(ib2dm-ib1),1)/&
                 &( (omega_i(:)-ebands%eig(ib1dm,ik_ibz,1))*(omega_i(:)-ebands%eig(ib2dm,ik_ibz,1)) )&
                                    +conjg( sigcme_k(:,1+(ib2dm-ib1),1+(ib1dm-ib1),1)/&
                 &( (omega_i(:)-ebands%eig(ib1dm,ik_ibz,1))*(omega_i(:)-ebands%eig(ib2dm,ik_ibz,1)) ) ) ) )
     ! Dji = Dij^*
     rdm_k(1+(ib2dm-ib1),1+(ib1dm-ib1))=conjg(rdm_k(1+(ib1dm-ib1),1+(ib2dm-ib1)))
   end do
 end do

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
!! rdm_k=density matrix, matrix (i,j), where i and j belong to the k-point k (see m_sigma_driver for more details).
!! occs = array containing the occ numbers for a given k-point occs(nband,k-point).
!! ebands=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!! checksij=check the orthonormality of the nat. orbitals
!!
!! OUTPUT
!! Compute the nat. orbitals and occ. numbers from the rdm_k matrix (for exchange and correlations)
!!
!! SOURCE

subroutine natoccs(ib1,ib2,rdm_k,nateigv,occs,ebands,ik_ibz,iinfo,checksij)
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2,ik_ibz,iinfo
 integer,intent(in),optional :: checksij
 type(ebands_t),target,intent(in) :: ebands
!arrays
 real(dp),intent(inout) :: occs(:,:)
 complex(dpc),intent(inout) :: rdm_k(:,:),nateigv(:,:,:,:)
!Local variables ------------------------------
!scalars
 integer:: ndim,ib1dm,ib2dm,ib3dm,lwork,info
 logical:: check_Sijmat
 character(len=500) :: msg
 real(dp) :: toccs_k,tol10
 complex(dp) :: Sib1k_ib2k
!arrays
 integer :: units(2)
 real(dp),allocatable :: occs_tmp(:),occs_tmp2(:),rwork(:)
 complex(dpc),allocatable :: work(:),tmp_mat(:,:),eigenvect(:,:)
!************************************************************************

 check_Sijmat=.false.; if (present(checksij)) check_Sijmat=.true.
 units = [std_out, ab_out]
 tol10=1.0e-10

 ndim=ib2-ib1+1
 lwork=2*ndim-1
 ABI_MALLOC(occs_tmp,(ndim))
 ABI_MALLOC(occs_tmp2,(ndim))
 ABI_MALLOC(work,(lwork))
 ABI_MALLOC(tmp_mat,(ndim,ndim))
 ABI_MALLOC(eigenvect,(ndim,ndim))
 ABI_MALLOC(rwork,(3*ndim-2))

 tmp_mat=zero
 do ib2dm=1,ndim
   do ib1dm=ib2dm,ndim
     tmp_mat(ib1dm,ib2dm)=rdm_k(ib1dm,ib2dm)
     ! Dji = Dij^*
     tmp_mat(ib2dm,ib1dm)=conjg(tmp_mat(ib1dm,ib2dm))
   end do
 end do

 work=zero
 occs_tmp=zero
 info=0
 call zheev('v','u',ndim,tmp_mat,ndim,occs_tmp,work,lwork,rwork,info)
 ABI_CHECK(info == 0, sjoin("Failed the diagonalization of the updated GW 1-RDM with info:", itoa(info)))

 ! Sort in descending order
 do ib1dm=1,ndim
  occs_tmp2(ib1dm)=occs_tmp(ndim-(ib1dm-1))
  do ib2dm=1,ndim
   eigenvect(ib2dm,ib1dm)=tmp_mat(ib2dm,(ndim-(ib1dm-1)))
  end do
  if (abs(occs_tmp2(ib1dm))<tol10) then
    occs_tmp2(ib1dm)=zero
  end if
 end do

 ! Check orthonormality?
 if (check_Sijmat) then
   do ib1dm=1,ndim
     do ib2dm=1,ib1dm
       Sib1k_ib2k=czero
       do ib3dm=1,ndim
         Sib1k_ib2k=Sib1k_ib2k+conjg(eigenvect(ib3dm,ib1dm))*eigenvect(ib3dm,ib2dm)
       end do
       if (ib1dm==ib2dm) then
         if(abs(Sib1k_ib2k-cmplx(one,zero))>tol10) then
           write(msg,'(a45,i5,a1,i5,f10.5)') 'Large deviation from identity for bands ',ib1dm,' ',ib2dm,real(Sib1k_ib2k)
           call wrtout(std_out,msg)
         endif
       else
         if (abs(Sib1k_ib2k)>tol10) then
           write(msg,'(a45,i5,a1,i5,f10.5)') 'Large deviation from identity for bands ',ib1dm,' ',ib2dm,real(Sib1k_ib2k)
           call wrtout(std_out,msg)
         end if
       end if
     end do
   end do
 end if

 ! Print results
 if (info==0) then
   if (iinfo==0) then
     write(msg,'(a51,3f10.5)') 'Occs. after updating with Sx-Vxc corr. at k-point:',ebands%kptns(1:,ik_ibz)
   else
     write(msg,'(a51,3f10.5)') 'Occs. after updating with S_c correct. at k-point:',ebands%kptns(1:,ik_ibz)
   endif
   call wrtout(units, msg)
   ib1dm=ndim-(ndim/10)*10
   do ib2dm=1,(ndim/10)*10,10
     write(msg,'(f11.5,9f10.5)') occs_tmp2(ib2dm:ib2dm+9)
     call wrtout(units, msg)
   end do
   ib1dm=(ndim/10)*10+1
   write(msg,'(f11.5,*(f10.5))') occs_tmp2(ib1dm:)
   call wrtout(units, msg)
 else
   write(msg,'(a36,3f10.5)') 'Error computing occs. for k-point: ',ebands%kptns(1:,ik_ibz)
   call wrtout(units, msg)
 end if

 ! Store natural orbital eigenvectors matrix and occs. Also compute total number of electrons for this k-point
 toccs_k=zero
 do ib1dm=1,ndim
   do ib2dm=1,ndim
     nateigv(ib1+(ib1dm-1),ib1+(ib2dm-1),ik_ibz,1)=eigenvect(ib1dm,ib2dm)
   end do
   occs(ib1+(ib1dm-1),ik_ibz)=occs_tmp2(ib1dm)  ! Overwrite the initial KS-DFT occs from ib1 to ib2
   toccs_k=toccs_k+occs_tmp2(ib1dm)
 end do

 write(msg,'(a22,i5,a3,i5,a21,f10.5)') ' Total occ. from band ',ib1,' to', ib2,' at current k-point: ',toccs_k
 call wrtout(units, msg)
 write(msg,'(a5)') ' '
 call wrtout(units, msg)

 ABI_FREE(rwork)
 ABI_FREE(work)
 ABI_FREE(tmp_mat)
 ABI_FREE(eigenvect)
 ABI_FREE(occs_tmp)
 ABI_FREE(occs_tmp2)

end subroutine natoccs
!!***

!!****f* ABINIT/update_hdr_bst
!! NAME
!! update_hdr_bst
!!
!! FUNCTION
!! Update the Hdr for the WFK and DEN files and the occ. numbers in the ebands file for a given k-point
!!
!! INPUTS
!! Wfd<wfd_t>=Datatype gathering data on QP amplitudes.
!! ngfft_in(18)=information on the fine FFT grid used for densities and potentials.
!! b1gw=min band for given k in the interval where we update.
!! b2gw=max band for given k in the interval where we update.
!! occs= array containing the occ numbers for a given k-point occs_ks(nband,k-point).
!! ebands=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!!
!! OUTPUT
!! Updated Hdr and ebands information
!!
!! SOURCE

subroutine update_hdr_bst(Wfd,occs,b1gw,b2gw,ebands,Hdr,ngfft_in)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: b1gw,b2gw
 integer,intent(in),dimension(3) :: ngfft_in
 type(ebands_t),target,intent(inout) :: ebands
 type(Hdr_type),intent(inout) :: Hdr
 class(wfd_t),intent(in) :: Wfd
!arrays
 real(dp),intent(in) :: occs(:,:)
!Local variables ------------------------------
!scalars
 integer :: ib1dm,ib2dm,dim_bands,ikpoint
!arrays
!************************************************************************

 ! ebands occ (QP_ebands ones) are changed and never recoverd
 do ikpoint=1,ebands%nkpt
   ebands%occ(b1gw:b2gw,ikpoint,1) = occs(b1gw:b2gw,ikpoint) ! Spins summed, occ in [0:2]
 enddo
 ABI_COMMENT("QP_ebands: occupancies were updated with nat. orb. ones")
 if ((size(Hdr%occ(:))/ebands%nkpt) < (b2gw-b1gw+1)) then
   !Actually, we should never reach this point because the code should stop during Wfd initialization in m_sigma_driver
   ABI_ERROR("Impossible to use the existing read WFK to build a new one!")
 end if

 ! Update occ in Hdr before printing
 ib1dm=1
 do ikpoint=1,ebands%nkpt
   dim_bands=size(ebands%occ(:,ikpoint,1))
   do ib2dm=1,dim_bands
     Hdr%occ(ib1dm)=ebands%occ(ib2dm,ikpoint,1) ! Because Hdr%occ is a 1-D array
     ib1dm=ib1dm+1
   end do
 end do

 Hdr%npwarr(:)=Wfd%npwarr(:)                                   ! Use the npw and ngfft = ones used in GW calc
 Hdr%ngfft(1:3)=ngfft_in(1:3)
 ABI_COMMENT("Hdr_sigma: occupancies, npw, and ngfft were updated")

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
!! ebands=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!!
!! OUTPUT
!! Print the total (averaged) occ. = sum_k weight_k * Nelec_k
!!
!! SOURCE

subroutine print_tot_occ(ebands)

!Arguments ------------------------------------
 type(ebands_t),intent(in) :: ebands

!Local variables-------------------------------
!scalars
 character(len=500) :: msg
 integer :: ik,spin, units(2)
 real(dp) :: wtk,occ_bks,tot_occ
! *************************************************************************

 units = [std_out, ab_out]

 tot_occ=zero

 do spin=1,ebands%nsppol
   do ik=1,ebands%nkpt
     wtk = ebands%wtk(ik)
     occ_bks = sum(ebands%occ(:,ik,spin))
     !if (sigma%nsig_ab==1) then ! Only closed-shell restricted is programed
     tot_occ=tot_occ+occ_bks*wtk
     !end if
   end do
 end do

 write(msg,'(a1)') ' '
 call wrtout(units, msg)
 write(msg,'(a39,f10.5)') ' Total averaged occ. from all k-points: ',tot_occ
 call wrtout(units, msg)
 write(msg,'(a1)') ' '
 call wrtout(units, msg)

end subroutine print_tot_occ
!!***

!!****f* ABINIT/get_chkprdm
!! NAME
!! get_chkprdm
!!
!! FUNCTION
!!  Read all checkpoint files built on previous runs
!!
!! INPUTS
!! Wfd<wfd_t>=Wave function descriptor see file 69_wfd/m_wfd.F90
!! Kmesh <kmesh_t>=Structure describing the k-point sampling.
!! Sigp<sigparams_t>=Parameters governing the self-energy calculation.
!! ebands=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!! occs = occ. numbers array occs(Wfd%mband,Wfd%nkibz)
!! nateigv = natural orbital eigenvectors nateigv(Wfd%mband,Wfd%mband,Wfd%nkibz,Sigp%nsppol))
!! sigmak_todo = integer array initialized to 1 and its components are set to 0 if the kpoint
!! is read from the checkpoint sigmak_todo(Wfd%nkibz)
!! my_rank = rank of the mpi process.
!!
!! OUTPUT
!! occ are updated if they are read from any checkpoint file
!! nateigv are stored if they are read from any checkpoint file
!! sigmak_todo components set to 1 if the kpoint is read from any checkpoint file
!!
!! SOURCE

subroutine get_chkprdm(Wfd,Kmesh,Sigp,ebands,occs,nateigv,sigmak_todo,my_rank,gw1rdm_fname_in)
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_rank
 class(wfd_t),intent(in) :: Wfd
 type(kmesh_t),intent(in) :: Kmesh
 type(sigparams_t),intent(in) :: Sigp
 type(ebands_t),intent(in) :: ebands
 character(len=fnlen),intent(in) :: gw1rdm_fname_in
!arrays
 integer,intent(inout) :: sigmak_todo(:)
 real(dp),intent(inout) :: occs(:,:)
 complex(dpc),intent(inout) :: nateigv(:,:,:,:)
!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,iunit=666314
 integer :: ierr,ib1,ib2,ib3,ikcalc,istat,ik_ibz,ik_ibz_read,iread,iread_eigv
 real(dp) :: auxl_read
 character(len=fnlen) :: gw1rdm_fname
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: occ_tmp(:),eigvect_tmp(:)
! *************************************************************************

 if (my_rank==master) then
   iread_eigv=Wfd%mband
   iread_eigv=iread_eigv*(2*iread_eigv)
   ABI_MALLOC(occ_tmp,(Wfd%mband))
   ABI_MALLOC(eigvect_tmp,(iread_eigv))

   do ikcalc=1,Sigp%nkptgw
     ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Irred k-point for GW
     if(ik_ibz<10) then
       write(gw1rdm_fname,"(a,i1)") trim(gw1rdm_fname_in),ik_ibz
     else if(ik_ibz<100 .and. ik_ibz>=10) then
       write(gw1rdm_fname,"(a,i2)") trim(gw1rdm_fname_in),ik_ibz
     else if(ik_ibz<1000 .and. ik_ibz>=100) then
       write(gw1rdm_fname,"(a,i3)") trim(gw1rdm_fname_in),ik_ibz
     else
       ABI_ERROR("The maximum k-point label for the checkpoint file to read is 999.")
     end if
     write(msg,'(a1)')' '
     call wrtout(std_out,msg)
     write(msg,'(a25,a)')' Reading checkpoint file ',gw1rdm_fname
     call wrtout(std_out,msg)
     write(msg,'(a1)')' '
     call wrtout(std_out,msg)
     occ_tmp(:)=zero;eigvect_tmp(:)=zero;
     open(unit=iunit,form='unformatted',file=gw1rdm_fname,iostat=istat,status='old')
     iread=0;ik_ibz_read=0;
     if (istat==0) then
       do
         if (iread<Wfd%mband) then
           iread=iread+1
           read(iunit,iostat=istat) auxl_read
           if (istat==0) then
             occ_tmp(iread)=auxl_read
           end if
         else if (iread<(iread_eigv+Wfd%mband)) then
           iread=iread+1
           read(iunit,iostat=istat) auxl_read
           if (istat==0) then
             eigvect_tmp(iread-Wfd%mband)=auxl_read
           end if
         else
           read(iunit,iostat=istat) ik_ibz_read
           if (istat==0 .and. ik_ibz_read/=0) then
            iread=0
            sigmak_todo(ik_ibz_read)=0
            ib3=1
            do ib1=1,Wfd%mband
              occs(ib1,ik_ibz_read)=occ_tmp(ib1)
              do ib2=1,Wfd%mband
                nateigv(ib1,ib2,ik_ibz_read,1)=cmplx(eigvect_tmp(ib3),eigvect_tmp(ib3+1))
                ib3=ib3+2
              end do
            end do
            ik_ibz_read=0
            occ_tmp=zero;eigvect_tmp=zero;
           end if
         end if
         if(istat/=0) then
           exit
         end if
       end do
     end if
     close(iunit)
   end do
   write(msg,'(a1)')' '
   call wrtout(std_out,msg)
   write(msg,'(a49)')' List of k-points read from all checkpoint files '
   call wrtout(std_out,msg)
   do ikcalc=1,Sigp%nkptgw
     ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Irred k-point for GW
     if (sigmak_todo(ik_ibz)==0) then
       write(msg,'(3f10.5)') ebands%kptns(1:,ik_ibz)
       call wrtout(std_out,msg)
     end if
   enddo
   write(msg,'(a1)')' '
   call wrtout(std_out,msg)
   ABI_FREE(occ_tmp)
   ABI_FREE(eigvect_tmp)
 end if

! Broadcast from master the information stored in occs and nateigv to all processes.
 call xmpi_barrier(Wfd%comm)
 ierr=0
 call xmpi_bcast(sigmak_todo(:),master,Wfd%comm,ierr)
 if(ierr/=0) then
   ABI_ERROR("Error distributing the sigmak_todo table.")
 endif
 call xmpi_bcast(occs(:,:),master,Wfd%comm,ierr)
 if(ierr/=0) then
   ABI_ERROR("Error distributing the occs read from checkpoint file(s).")
 endif
 call xmpi_bcast(nateigv(:,:,:,:),master,Wfd%comm,ierr)
 if(ierr/=0) then
   ABI_ERROR("Error distributing the natural orbital eigenvectors read from checkpoint file(s).")
 endif

end subroutine get_chkprdm
!!***

!!****f* ABINIT/print_chkprdm
!! NAME
!! print_chkprdm
!!
!! FUNCTION
!!  Write the checkpoint file for a given k-point
!!
!! INPUTS
!! Wfd<wfd_t>=Wave function descriptor see file 69_wfd/m_wfd.F90
!! occs = occ. numbers array occs(Wfd%mband,Wfd%nkibz)
!! nateigv = natural orbital eigenvectors nateigv(Wfd%mband,Wfd%mband,Wfd%nkibz,Sigp%nsppol))
!! my_rank = rank of the mpi process.
!! gw1rdm_fname_out = name of the gw1rdm checkpoint out without k-point extension
!!
!! OUTPUT
!!
!! SOURCE

subroutine print_chkprdm(Wfd,occs,nateigv,ik_ibz,my_rank,gw1rdm_fname_out)
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_ibz,my_rank
 class(wfd_t),intent(in) :: Wfd
 character(len=fnlen),intent(in) :: gw1rdm_fname_out
!arrays
 real(dp),intent(in) :: occs(:,:)
 complex(dpc),intent(in) :: nateigv(:,:,:,:)
!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,iunit=666314
 integer :: iwrite,iwrite2
 character(len=fnlen) :: gw1rdm_fname
 character(len=500) :: msg
!arrays
! *************************************************************************

 if (my_rank==master) then
   if(ik_ibz<10) then
     write(gw1rdm_fname,"(a,i1)") trim(gw1rdm_fname_out),ik_ibz
   else if(ik_ibz<100 .and. ik_ibz>=10) then
     write(gw1rdm_fname,"(a,i2)") trim(gw1rdm_fname_out),ik_ibz
   else if(ik_ibz<1000 .and. ik_ibz>=100) then
     write(gw1rdm_fname,"(a,i3)") trim(gw1rdm_fname_out),ik_ibz
   else
     ABI_ERROR("The maximum k-point label for the checkpoint file to write is 999.")
   end if
   write(msg,'(a1)')' '
   call wrtout(std_out,msg)
   write(msg,'(a25,a)')' Writing checkpoint file ',gw1rdm_fname
   call wrtout(std_out,msg)
   write(msg,'(a1)')' '
   call wrtout(std_out,msg)
   open(unit=iunit,form='unformatted',file=gw1rdm_fname)
   do iwrite=1,Wfd%mband
     write(iunit) occs(iwrite,ik_ibz)
   end do
   do iwrite=1,Wfd%mband
     do iwrite2=1,Wfd%mband
       write(iunit) real(nateigv(iwrite,iwrite2,ik_ibz,1))
       write(iunit) aimag(nateigv(iwrite,iwrite2,ik_ibz,1))
     end do
   end do
   write(iunit) ik_ibz
   close(iunit)
 end if

 call xmpi_barrier(Wfd%comm)

end subroutine print_chkprdm
!!***

!!****f* ABINIT/change_matrix
!! NAME
!! change_matrix
!!
!! FUNCTION
!!  Transform integrals from KS -> NO and NO -> KS orbitals
!!
!!   Transform <NO_i|K[NO]|NO_j> -> <KS_i|K[NO]|KS_j>,
!!             <KS_i|J[NO]|KS_j> -> <NO_i|J[NO]|NO_j>,
!!   and         <KS_i|T|KS_j>   ->   <NO_i|T|NO_j>
!!
!!
!! INPUTS
!! Kmesh <kmesh_t>=Structure describing the k-point sampling.
!! Sigp<sigparams_t>=Parameters governing the self-energy calculation.
!! nateigv = natural orbital eigenvectors nateigv(Wfd%mband,Wfd%mband,Wfd%nkibz,Sigp%nsppol))
!!
!! OUTPUT
!!  Mels
!!   %kinetic=matrix elements of $t$.
!!   %vhartr =matrix elements of $v_H$.
!! Sr=sigma_t (see the definition of this structured datatype)
!!
!! SOURCE

subroutine change_matrix(Sigp,Sr,Mels,Kmesh,nateigv)

!Arguments ------------------------------------
!scalars
 type(kmesh_t),intent(in) :: Kmesh
 type(sigparams_t),intent(in) :: Sigp
 type(sigma_t),intent(inout) :: Sr
 type(melements_t),intent(inout) :: Mels
!arrays
 complex(dpc),intent(in) :: nateigv(:,:,:,:)
!Local variables-------------------------------
!scalars
 integer :: ikcalc,ik_ibz,ib1,ib2,ib1dm,ib2dm
!arrays
 complex(dpc),allocatable :: mat2rot(:,:),Umat(:,:)
! *************************************************************************

  do ikcalc=1,Sigp%nkptgw
    ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Index of the irreducible k-point for GW
    ib1=MINVAL(Sigp%minbnd(ikcalc,:))       ! min and max band indices for GW corrections (for this k-point)
    ib2=MAXVAL(Sigp%maxbnd(ikcalc,:))
    ABI_MALLOC(mat2rot,(ib2-ib1+1,ib2-ib1+1))
    ABI_MALLOC(Umat,(ib2-ib1+1,ib2-ib1+1))
    ! <NO_i|K[NO]|NO_j> -> <KS_i|K[NO]|KS_j>
    do ib1dm=1,ib2-ib1+1
      do ib2dm=1,ib2-ib1+1
        Umat(ib1dm,ib2dm)=nateigv(ib1+(ib1dm-1),ib1+(ib2dm-1),ik_ibz,1)
        mat2rot(ib1dm,ib2dm)=Sr%x_mat(ib1+(ib1dm-1),ib1+(ib2dm-1),ik_ibz,1)
      end do
    end do
    call rotate_ks_no(ib1,ib2,mat2rot,Umat,0)
    do ib1dm=1,ib2-ib1+1
      do ib2dm=1,ib2-ib1+1
        Sr%x_mat(ib1+(ib1dm-1),ib1+(ib2dm-1),ik_ibz,1)=mat2rot(ib1dm,ib2dm)
      end do
    end do
    ! <KS_i|J[NO]|KS_j> -> <NO_i|J[NO]|NO_j>
    do ib1dm=1,ib2-ib1+1
      do ib2dm=1,ib2-ib1+1
        mat2rot(ib1dm,ib2dm)=Mels%vhartree(ib1+(ib1dm-1),ib1+(ib2dm-1),ik_ibz,1)
      end do
    end do
    call rotate_ks_no(ib1,ib2,mat2rot,Umat,1)
    do ib1dm=1,ib2-ib1+1
      do ib2dm=1,ib2-ib1+1
        Mels%vhartree(ib1+(ib1dm-1),ib1+(ib2dm-1),ik_ibz,1)=mat2rot(ib1dm,ib2dm)
      end do
    end do
    ! <KS_i|T|KS_j> -> <NO_i|T|NO_j>
    do ib1dm=1,ib2-ib1+1
      do ib2dm=1,ib2-ib1+1
        mat2rot(ib1dm,ib2dm)=Mels%kinetic(ib1+(ib1dm-1),ib1+(ib2dm-1),ik_ibz,1)
      end do
    end do
    call rotate_ks_no(ib1,ib2,mat2rot,Umat,1)
    do ib1dm=1,ib2-ib1+1
      do ib2dm=1,ib2-ib1+1
        Mels%kinetic(ib1+(ib1dm-1),ib1+(ib2dm-1),ik_ibz,1)=mat2rot(ib1dm,ib2dm)
      end do
    end do
    ABI_FREE(Umat)
    ABI_FREE(mat2rot)
  end do
end subroutine change_matrix
!!***

!!****f* ABINIT/print_total_energy
!! NAME
!! print_total_energy
!!
!! FUNCTION
!!  Print total energy and energy components
!!
!! INPUTS
!! all energy terms are self-explanatory
!!
!! OUTPUT
!!
!! SOURCE

subroutine print_total_energy(ekin_energy,evext_energy,evextnl_energy,e_corepsp,eh_energy,ex_energy,&
                              exc_mbb_energy,e_ewald,etot,etot2,den_int)
!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: ekin_energy,evext_energy,evextnl_energy,e_corepsp,eh_energy,ex_energy
 real(dp),intent(in) :: exc_mbb_energy,e_ewald,etot,etot2,den_int

!Local variables-------------------------------
 character(len=500) :: msg
 integer :: units(2)

!************************************************************************

 units = [std_out, ab_out]

 write(msg,'(a1)')' '
 call wrtout(units, msg)
 write(msg,'(a98)')'---------------------------------------------------------------&
         &----------------------------------'
 call wrtout(units, msg)
 write(msg,'(a,f10.5,a,f10.3,a)')' Ekinetic   = : ',ekin_energy,' Ha ,',ekin_energy*Ha_eV,' eV'
 call wrtout(units, msg)
 write(msg,'(a,f10.5,a,f10.3,a)')' Evext_l    = : ',evext_energy,' Ha ,',evext_energy*Ha_eV,' eV'
 call wrtout(units, msg)
 write(msg,'(a,f10.5,a,f10.3,a)')' Evext_nl   = : ',evextnl_energy,' Ha ,',evextnl_energy*Ha_eV,' eV'
 call wrtout(units, msg)
 write(msg,'(a,f10.5,a,f10.3,a)')' Epsp_core  = : ',e_corepsp,' Ha ,',e_corepsp*Ha_eV,' eV'
 call wrtout(units, msg)
 write(msg,'(a,f10.5,a,f10.3,a)')' Ehartree   = : ',eh_energy,' Ha ,',eh_energy*Ha_eV,' eV'
 call wrtout(units, msg)
 write(msg,'(a,f10.5,a,f10.3,a)')' Ex[SD]     = : ',ex_energy,' Ha ,',ex_energy*Ha_eV,' eV'
 call wrtout(units, msg)
 write(msg,'(a,f10.5,a,f10.3,a)')' Exc[MBB]   = : ',exc_mbb_energy,' Ha ,',exc_mbb_energy*Ha_eV,' eV'
 call wrtout(units, msg)
 write(msg,'(a,f10.5,a,f10.3,a)')' Enn        = : ',e_ewald,' Ha ,',e_ewald*Ha_eV,' eV'
 call wrtout(units, msg)
 write(msg,'(a98)')'-----------------------------------------------------------------&
         &--------------------------------'
 call wrtout(units, msg)
 write(msg,'(a,f10.5,a,f10.3,a)')' Etot[SD]   = : ',etot,' Ha ,',etot*Ha_eV,' eV'
 call wrtout(units, msg)
 write(msg,'(a,f10.5,a,f10.3,a)')' Etot[MBB]  = : ',etot2,' Ha ,',etot2*Ha_eV,' eV'
 call wrtout(units, msg)
 write(msg,'(a,f10.5,a,f10.3,a)')' Vee[SD]    = : ',(ex_energy+eh_energy),' Ha ,',(ex_energy+eh_energy)*Ha_eV,' eV'
 call wrtout(units, msg)
 write(msg,'(a,f10.5,a,f10.3,a)')' Vee[MBB]   = : ',(exc_mbb_energy+eh_energy),' Ha ,',&
         &(exc_mbb_energy+eh_energy)*Ha_eV,' eV'
 call wrtout(units, msg)
 write(msg,'(a,f10.5)')  ' Density    = : ',den_int
 call wrtout(units, msg)
 write(msg,'(a)')' Vee[SD] (= Ehartree + Ex[SD]) energy obtained using GW 1-RDM:'
 call wrtout(units, msg)
 write(msg,'(a)')' Vee[MBB] (= Ehartree + Exc[MBB]) energy obtained using GW 1-RDM:'
 call wrtout(units, msg)
 write(msg,'(a98)')'-------------------------------------------------------------------&
         &------------------------------'
 call wrtout(units, msg)

end subroutine print_total_energy
!!***

!!****f* ABINIT/print_band_energies
!! NAME
!! print_band_energies
!!
!! FUNCTION
!!  Print updated band energies
!!
!!
!! INPUTS
!! Kmesh <kmesh_t>=Structure describing the k-point sampling.
!! Sigp<sigparams_t>=Parameters governing the self-energy calculation.
!!  Mels
!!   %kinetic=matrix elements of $t$.
!!   %vhartr =matrix elements of $v_H$.
!! Sr=sigma_t (see the definition of this structured datatype)
!! ebands=<ebands_t>=Datatype gathering info on the QP energies (KS if one shot)
!!  eig(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=KS or QP energies for k-points, bands and spin
!!  occ(Sigp%nbnds,Kmesh%nibz,Wfd%nsppol)=occupation numbers, for each k point in IBZ, each band and spin
!!
!! OUTPUT
!!
!! SOURCE

subroutine print_band_energies(b1gw,b2gw,Sr,Sigp,Mels,Kmesh,ebands,new_hartr,old_purex)
!Arguments ------------------------------------
!scalars
 type(kmesh_t),intent(in) :: Kmesh
 type(sigparams_t),intent(in) :: Sigp
 type(sigma_t),intent(in) :: Sr
 type(ebands_t),intent(in) :: ebands
 type(melements_t),intent(in) :: Mels
 integer,intent(in) :: b1gw,b2gw
!arrays
 complex(dpc),intent(in) :: old_purex(:,:),new_hartr(:,:)
!Local variables-------------------------------
!scalars
 integer :: ib,ikcalc,ik_ibz, units(2)
 real(dp) :: eik_new
 complex(dpc) :: delta_band_ibik
 character(len=500) :: msg
!************************************************************************

 units = [std_out, ab_out]

 write(msg,'(a1)')  ' '
 call wrtout(units, msg)
 write(msg,'(a42)')  ' Computing band corrections Delta eik (eV)'
 call wrtout(units, msg)
 write(msg,'(a42)')  ' -----------------------------------------'
 call wrtout(units, msg)
 write(msg,'(a1)')  ' '
 call wrtout(units, msg)
 write(msg,'(a1)')  ' '
 call wrtout(units, msg)
 write(msg,'(a110)') ' Band corrections Delta eik = <KS_i|K[NO]-a*K[KS]+vH[NO]&
       &-vH[KS]-Vxc[KS]|KS_i> and eik^new = eik^GS + Delta eik'
 call wrtout(units, msg)
 write(msg,'(a1)')  ' '
 call wrtout(units, msg)
 do ikcalc=1,Sigp%nkptgw
   ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Index of the irreducible k-point for GW
   write(msg,'(a127)')'---------------------------------------------------------&
           &--------------------------------------------------------------------'
   call wrtout(units, msg)
   write(msg,'(a)')'  k-point  band      eik^GS        eik^new     Delta eik  &
     &       K[NO]       a*K[KS]         Vxc[KS]       vH[NO]        vH[KS]'
   call wrtout(units, msg)
   do ib=b1gw,b2gw
     delta_band_ibik=(new_hartr(ib,ikcalc)-Mels%vhartree(ib,ib,ik_ibz,1))&
     &+Sr%x_mat(ib,ib,ik_ibz,1)-Mels%vxcval(ib,ib,ik_ibz,1)-old_purex(ib,ikcalc)
     eik_new=real(ebands%eig(ib,ik_ibz,1))+real(delta_band_ibik)
     write(msg,'(i5,4x,i5,8(4x,f10.3))') &
     & ik_ibz,ib,real(ebands%eig(ib,ik_ibz,1))*Ha_eV,eik_new*Ha_eV,real(delta_band_ibik)*Ha_eV,&
     & real(Sr%x_mat(ib,ib,ik_ibz,1))*Ha_eV,real(old_purex(ib,ikcalc))*Ha_eV,&
     & real(Mels%vxcval(ib,ib,ik_ibz,1))*Ha_eV,&
     & real(new_hartr(ib,ikcalc))*Ha_eV,real(Mels%vhartree(ib,ib,ik_ibz,1))*Ha_eV
     call wrtout(units, msg)
   enddo
 enddo
 write(msg,'(a127)')'---------------------------------------------------------&
         &--------------------------------------------------------------------'
 call wrtout(units, msg)

end subroutine print_band_energies
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

!!****f* ABINIT/printrdm_k
!! NAME
!! printrdm_k
!!
!! FUNCTION
!! Print the DM1 matrix
!!
!! INPUTS
!! ib1=min band.
!! ib2=max band.
!! rdm_k=array containing the 1-RDM matrix
!!
!! OUTPUT
!! Print the 1-RDM matrix
!! SOURCE

subroutine printrdm_k(ib1,ib2,rdm_k) ! Only used for debug on this file, do not use it with large arrays!
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ib1,ib2
!arrays
 complex(dpc),intent(in) :: rdm_k(:,:)
!Local variables ------------------------------
!scalars
 integer::ib1dm
 character(len=500) :: msg
!arrays
!************************************************************************

 do ib1dm=ib1,ib2
   write(msg,'(*(f12.5))') real(rdm_k(ib1dm,ib1:ib2))
   call wrtout(std_out, msg)
 end do

end subroutine printrdm_k
!!***

end module m_gwrdm
!!***
