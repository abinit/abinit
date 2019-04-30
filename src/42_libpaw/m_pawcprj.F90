!{\src2tex{textfont=tt}}
!!****m* abinit/m_pawcprj
!! NAME
!!  m_pawcprj
!!
!! FUNCTION
!!  This module contains functions used to manipulate variables of
!!   structured datatype pawcprj_type.
!!   pawcprj_type variables are <p_lmn|Cnk> projected quantities,
!!   where |p_lmn> are non-local projectors
!!         |Cnk> are wave functions
!!
!! COPYRIGHT
!! Copyright (C) 2012-2019 ABINIT group (MT,JWZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

module m_pawcprj

 USE_DEFS
 USE_MSG_HANDLING
 USE_MPI_WRAPPERS
 USE_MEMORY_PROFILING

 use m_pawtab, only : pawtab_type

 implicit none

 private
!!***

!!****t* m_pawcprj/pawcprj_type
!! NAME
!! pawcprj_type
!!
!! FUNCTION
!! This structured datatype contains <p_lmn|Cnk> projected scalars and derivatives
!!             where |p_lmn> are non-local projectors for a given atom
!!                   |Cnk> is a wave function
!! Used only for PAW calculations.
!!
!! SOURCE

 type,public :: pawcprj_type

!Integer scalars

  integer :: ncpgr=0
   ! Number of gradients of cp=<p_lmn|Cnk>

  integer :: nlmn=0
   ! Number of (l,m,n) non-local projectors

!Real (real(dp)) arrays

  real(dp), allocatable :: cp (:,:)
   ! cp(2,nlmn)
   ! <p_lmn|Cnk> projected scalars for a given atom and wave function

  real(dp), allocatable :: dcp (:,:,:)
   ! dcp(2,ncpgr,nlmn)
   ! derivatives of <p_lmn|Cnk> projected scalars for a given atom and wave function

 end type pawcprj_type

!public procedures.
 public :: pawcprj_alloc          ! Allocation
 public :: pawcprj_free           ! Deallocation
 public :: pawcprj_set_zero       ! Set to zero all arrays in a cprj datastructure
 public :: pawcprj_copy           ! Copy a cprj datastructure into another
 public :: pawcprj_axpby          ! cprjy(:,:) <- alpha.cprjx(:,:)+beta.cprjy(:,:)
 public :: pawcprj_zaxpby         ! cprjy(:,:) <- alpha.cprjx(:,:)+beta.cprjy(:,:), alpha and beta are COMPLEX scalars
 public :: pawcprj_conjg          ! cprj(:,:) <- conjugate(cprj(:,:))
 public :: pawcprj_symkn          ! construct cprj from that at a symmetry related k point
 public :: pawcprj_lincom         ! Compute a LINear COMbination of cprj datastructure:
 public :: pawcprj_output         ! Output a cprj. Useful for debugging.
 public :: pawcprj_get            ! Read the cprj for a given k-point from memory or from a temporary file
 public :: pawcprj_put            ! Write the cprj for a given set of (n,k) into memory or into a temporary file
 public :: pawcprj_reorder        ! Change the order of a cprj datastructure
 public :: pawcprj_mpi_allgather  ! Perform MPI_ALLGATHER on a pawcprj_type inside a MPI communicator.
 public :: pawcprj_bcast          ! Broadcast a pawcprj_type from master to all nodes inside a MPI communicator.
 public :: pawcprj_transpose      ! Transpose a cprj datastructure FOR A GIVEN (K,SPIN)
 public :: pawcprj_gather_spin    ! Collect spin distributed cprjs.
 public :: pawcprj_mpi_exch       ! Exchange a pawcprj_type between two processors inside a MPI communicator.
 public :: pawcprj_mpi_send       ! Send a pawcprj_type inside a MPI communicator.
 public :: pawcprj_mpi_recv       ! Receive a pawcprj_type inside a MPI communicator.
 public :: pawcprj_mpi_sum        ! Perform MPI_SUM on a pawcprj_type inside a MPI communicator.
 public :: pawcprj_getdim         ! Returns the number of lmn components in the <p_{lmn}^i|\psi> for the i-th atom.
 public :: paw_overlap            ! Compute the onsite contribution to the overlap between two states.
 public :: pawcprj_pack           ! Copy data from a cprj to a simple real buffer
 public :: pawcprj_unpack         ! Copy data from a simple real buffer to a cprj
!!***

CONTAINS

!===========================================================
!!***

!!****f* m_pawcprj/pawcprj_alloc
!! NAME
!! pawcprj_alloc
!!
!! FUNCTION
!! Allocation of a cprj datastructure
!!
!! INPUTS
!!  ncpgr=number of gradients to be allocated
!!  nlmn(:)=sizes of cprj%cp
!!
!! SIDE EFFECTS
!!  cprj(:,:) <type(pawcprj_type)>= cprj datastructure
!!
!! PARENTS
!!      berryphase_new,calc_optical_mels,calc_sigc_me,calc_sigx_me,calc_vhxc_me
!!      calc_wf_qp,cchi0,cchi0q0,cchi0q0_intraband,cgwf,chebfi,chern_number
!!      classify_bands,cohsex_me,ctocprj,d2frnl,m_datafordmft,debug_tools
!!      dfpt_accrho,dfpt_cgwf,dfpt_looppert,dfpt_nstpaw,dfpt_scfcv,dfpt_vtowfk
!!      dfpt_wfkfermi,dotprod_set_cgcprj,dotprodm_sumdiag_cgcprj,energy
!!      exc_build_block,exc_build_ham,exc_plot,extrapwf,fock2ACE,forstr
!!      forstrnps,getgh1c,getgh2c,getghc,getgsc,initberry,ks_ddiago
!!      lincom_cgcprj,m_electronpositron,m_fock,m_invovl,m_io_kss,m_pawcprj
!!      m_plowannier,m_shirley,m_wfd,make_grad_berry,nonlop,optics_paw
!!      optics_paw_core,outkss,partial_dos_fractions_paw,paw_symcprj,pawmkaewf
!!      pawmkrhoij,posdoppler,prep_calc_ucrpa,rf2_init,scfcv,setup_positron
!!      sigma,smatrix_pawinit,suscep_stat,update_e_field_vars,vtorho,vtowfk
!!      wf_mixing,wfd_pawrhoij,wfd_vnlpsi,wvl_hpsitopsi
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

 subroutine pawcprj_alloc(cprj,ncpgr,nlmn)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncpgr
!arrays
 integer,intent(in) :: nlmn(:)
 type(pawcprj_type),intent(inout) :: cprj(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,n1dim,n2dim,nn
 character(len=500) :: msg

! *************************************************************************

 n1dim=size(cprj,dim=1);n2dim=size(cprj,dim=2);nn=size(nlmn,dim=1)
 if (nn/=n1dim) then
   write(msg,*) 'wrong sizes (pawcprj_alloc)! :',nn,n1dim
   MSG_ERROR(msg)
 end if

 do jj=1,n2dim
   do ii=1,n1dim
     if (allocated(cprj(ii,jj)%cp)) then
       LIBPAW_DEALLOCATE(cprj(ii,jj)%cp)
     end if
     if (allocated(cprj(ii,jj)%dcp)) then
       LIBPAW_DEALLOCATE(cprj(ii,jj)%dcp)
     end if
     nn=nlmn(ii)
     cprj(ii,jj)%nlmn=nn
     LIBPAW_ALLOCATE(cprj(ii,jj)%cp,(2,nn))
     cprj(ii,jj)%cp=zero
     cprj(ii,jj)%ncpgr=ncpgr
     if (ncpgr>0) then
       LIBPAW_ALLOCATE(cprj(ii,jj)%dcp,(2,ncpgr,nn))
       cprj(ii,jj)%dcp=zero
     end if
   end do
 end do

end subroutine pawcprj_alloc
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_free
!! NAME
!! pawcprj_free
!!
!! FUNCTION
!! Deallocation of a cprj datastructure
!!
!! SIDE EFFECTS
!!  cprj(:,:) <type(pawcprj_type)>= cprj datastructure
!!
!! PARENTS
!!      berryphase_new,calc_optical_mels,calc_sigc_me,calc_sigx_me,calc_vhxc_me
!!      calc_wf_qp,cchi0,cchi0q0,cchi0q0_intraband,cgwf,chebfi,chern_number
!!      classify_bands,cohsex_me,ctocprj,d2frnl,datafordmft,debug_tools
!!      dfpt_accrho,dfpt_cgwf,dfpt_looppert,dfpt_nstpaw,dfpt_scfcv,dfpt_vtowfk
!!      dfpt_wfkfermi,dotprod_set_cgcprj,dotprodm_sumdiag_cgcprj,energy
!!      exc_build_block,exc_build_ham,exc_plot,extrapwf,fock2ACE,forstr
!!      forstrnps,getgh1c,getgh2c,getghc,getgsc,ks_ddiago,lincom_cgcprj
!!      m_efield,m_electronpositron,m_fock,m_gkk,m_invovl,m_io_kss,m_pawcprj
!!      m_phgamma,m_phpi,m_plowannier,m_scf_history,m_shirley,m_sigmaph,m_wfd
!!      make_grad_berry,nonlop,optics_paw,optics_paw_core,outkss
!!      partial_dos_fractions_paw,paw_symcprj,pawmkaewf,pawmkrhoij,posdoppler
!!      prep_calc_ucrpa,rf2_init,scfcv,setup_positron,sigma,smatrix_pawinit
!!      suscep_stat,update_e_field_vars,vtorho,vtowfk,wf_mixing,wfd_pawrhoij
!!      wfd_vnlpsi
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

 subroutine pawcprj_free(cprj)

!Arguments ------------------------------------
!scalars
!arrays
 type(pawcprj_type),intent(inout) :: cprj(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,n1dim,n2dim

! *************************************************************************

 n1dim=size(cprj,dim=1);n2dim=size(cprj,dim=2)

 do jj=1,n2dim
   do ii=1,n1dim
     if (allocated(cprj(ii,jj)%cp))  then
       LIBPAW_DEALLOCATE(cprj(ii,jj)%cp)
     end if
     if (allocated(cprj(ii,jj)%dcp))  then
       LIBPAW_DEALLOCATE(cprj(ii,jj)%dcp)
     end if
   end do
 end do

end subroutine pawcprj_free
!!***

!----------------------------------------------------------------------

!!      m_scf_history,suscep_stat
!!****f* m_pawcprj/pawcprj_set_zero
!! NAME
!! pawcprj_set_zero
!!
!! FUNCTION
!! Set to zero all arrays in a cprj datastructure
!!
!! SIDE EFFECTS
!!  cprj(:,:) <type(pawcprj_type)>= cprj datastructure
!!
!! PARENTS
!!      ctocprj,dfpt_cgwf,m_fock
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

 subroutine pawcprj_set_zero(cprj)

!Arguments ------------------------------------
!scalars
!arrays
 type(pawcprj_type),intent(inout) :: cprj(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,n1dim,n2dim

! *************************************************************************

 n1dim=size(cprj,dim=1);n2dim=size(cprj,dim=2)

 do jj=1,n2dim
   do ii=1,n1dim
     if (cprj(ii,jj)%nlmn>0)  cprj(ii,jj)%cp(:,:)=zero
     if (cprj(ii,jj)%ncpgr>0) cprj(ii,jj)%dcp(:,:,:)=zero
   end do
 end do

end subroutine pawcprj_set_zero
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_copy
!! NAME
!! pawcprj_copy
!!
!! FUNCTION
!! Copy a cprj datastructure into another
!!
!! INPUTS
!!  icpgr= (optional argument) if present, only component icpgr of
!!         input cprj gradient is copied into output cprj
!!         Not used if cprj(:,:)%ncpgr<icpgr
!!         -1 only copy cp
!!  cprj_in(:,:) <type(pawcprj_type)>= input cprj datastructure
!!
!! OUTPUT
!!  cprj_out(:,:) <type(pawcprj_type)>= output cprj datastructure
!!
!! NOTES
!!  MG: What about an option to report a pointer to cprj_in?
!!
!! PARENTS
!!      berryphase_new,calc_sigc_me,calc_sigx_me,cchi0q0,cchi0q0_intraband,cgwf
!!      chebfi,classify_bands,cohsex_me,corrmetalwf1,dfpt_looppert,dfpt_nstpaw
!!      dfpt_vtowfk,dfpt_wfkfermi,extrapwf,getgsc,m_electronpositron,m_fock
!!      m_pawcprj,m_wfd,make_grad_berry,nonlop,outkss,paw_symcprj,posdoppler
!!      prep_calc_ucrpa,setup_positron,vtowfk
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

 subroutine pawcprj_copy(cprj_in,cprj_out,&
&                    icpgr) ! optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: icpgr
!arrays
 type(pawcprj_type),intent(in) :: cprj_in(:,:)
 type(pawcprj_type),intent(inout) :: cprj_out(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,kk,n1dim_in,n1dim_out,n2dim_in,n2dim_out,ncpgr_in,ncpgr_out,nlmn
 logical :: has_icpgr,copy_dcp
 character(len=500) :: msg

! *************************************************************************

 n1dim_in=size(cprj_in,dim=1); n1dim_out=size(cprj_out,dim=1)
 n2dim_in=size(cprj_in,dim=2); n2dim_out=size(cprj_out,dim=2)
 ncpgr_in=cprj_in(1,1)%ncpgr;  ncpgr_out=cprj_out(1,1)%ncpgr

 if (n1dim_in/=n1dim_out) then
   write(msg,'(a,2(1x,i0))')" Error in pawcprj_copy: n1 wrong sizes ",n1dim_in,n1dim_out
   MSG_ERROR(msg)
 end if
 if (n2dim_in/=n2dim_out) then
   write(msg,'(a,2(1x,i0))')" Error in pawcprj_copy: n2 wrong sizes ",n2dim_in,n2dim_out
   MSG_ERROR(msg)
 end if
 if (ncpgr_in<ncpgr_out)  then
   write(msg,'(a,2(1x,i0))')" Error in pawcprj_copy: ncpgr wrong sizes ",ncpgr_in,ncpgr_out
   MSG_ERROR(msg)
 end if

!Check if icgr is present and if dcp have to be copy
 has_icpgr=present(icpgr)
 copy_dcp = .TRUE.
 if(has_icpgr)then
   copy_dcp = icpgr>=0
 end if

 do jj=1,n2dim_in
   do ii=1,n1dim_in
     nlmn=cprj_in(ii,jj)%nlmn
     cprj_out(ii,jj)%nlmn =nlmn
     do kk=1,nlmn
       cprj_out(ii,jj)%cp(1:2,kk)=cprj_in(ii,jj)%cp(1:2,kk)
     end do
   end do
 end do

 if (ncpgr_in>0.and.copy_dcp) then
   if (has_icpgr) has_icpgr=(ncpgr_out>0.and.icpgr>0.or.icpgr<=ncpgr_in)

   if (has_icpgr) then
     do jj=1,n2dim_in
       do ii=1,n1dim_in
         nlmn=cprj_in(ii,jj)%nlmn
         do kk=1,nlmn
           cprj_out(ii,jj)%dcp(1:2,1,kk)=cprj_in(ii,jj)%dcp(1:2,icpgr,kk)
         end do
       end do
     end do
   else
     if (ncpgr_out>=ncpgr_in) then
       do jj=1,n2dim_in
         do ii=1,n1dim_in
           nlmn=cprj_in(ii,jj)%nlmn
           do kk=1,nlmn
             cprj_out(ii,jj)%dcp(1:2,1:ncpgr_in,kk)=cprj_in(ii,jj)%dcp(1:2,1:ncpgr_in,kk)
           end do
         end do
       end do
     end if
   end if
 end if

end subroutine pawcprj_copy
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_axpby
!! NAME
!! pawcprj_axpby
!!
!! FUNCTION
!! Apply AXPBY (blas-like) operation with 2 cprj datastructures:
!!  cprjy(:,:) <- alpha.cprjx(:,:)+beta.cprjy(:,:)
!!  alpha and beta are REAL scalars
!!
!! INPUTS
!!  alpha,beta= alpha,beta REAL factors
!!  cprjx(:,:) <type(pawcprj_type)>= input cprjx datastructure
!!
!! SIDE EFFECTS
!!  cprjy(:,:) <type(pawcprj_type)>= input/output cprjy datastructure
!!
!! PARENTS
!!      chebfi,dfpt_cgwf,dfpt_wfkfermi,getdc1,m_invovl,wf_mixing
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

 subroutine pawcprj_axpby(alpha,beta,cprjx,cprjy)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: alpha,beta
!arrays
 type(pawcprj_type),intent(in) :: cprjx(:,:)
 type(pawcprj_type),intent(inout) :: cprjy(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,kk,n1dimx,n1dimy,n2dimx,n2dimy,ncpgrx,ncpgry,nlmn
 character(len=500) :: msg

! *************************************************************************

 n1dimy=size(cprjy,dim=1);n2dimy=size(cprjy,dim=2);ncpgry=cprjy(1,1)%ncpgr
 if (abs(alpha)>tol16) then
   n1dimx=size(cprjx,dim=1);n2dimx=size(cprjx,dim=2);ncpgrx=cprjx(1,1)%ncpgr
   msg = ""
   if (n1dimx/=n1dimy) msg = TRIM(msg)//"Error in pawcprj_axpby: n1 wrong sizes !"//ch10
   if (n2dimx/=n2dimy) msg = TRIM(msg)//"Error in pawcprj_axpby: n2 wrong sizes !"//ch10
   if (ncpgrx/=ncpgry) msg = TRIM(msg)//"Error in pawcprj_axpby: ncpgr wrong sizes !"//ch10
   if (LEN_TRIM(msg) > 0) then
     MSG_ERROR(msg)
   end if
 else
   n1dimx=0;n2dimx=0;ncpgrx=0
 end if

 if (abs(alpha)<=tol16) then
   do jj=1,n2dimy
     do ii=1,n1dimy
       nlmn=cprjy(ii,jj)%nlmn
       do kk=1,nlmn
         cprjy(ii,jj)%cp(1:2,kk)=beta*cprjy(ii,jj)%cp(1:2,kk)
       end do
     end do
   end do
   if (ncpgry>0) then
     do jj=1,n2dimy
       do ii=1,n1dimy
         nlmn=cprjy(ii,jj)%nlmn
         do kk=1,nlmn
           cprjy(ii,jj)%dcp(1:2,1:ncpgry,kk)=beta*cprjy(ii,jj)%dcp(1:2,1:ncpgry,kk)
         end do
       end do
     end do
   end if
 else if (abs(beta)<=tol16) then
   do jj=1,n2dimx
     do ii=1,n1dimx
       nlmn=cprjx(ii,jj)%nlmn
       cprjy(ii,jj)%nlmn=nlmn
       do kk=1,nlmn
         cprjy(ii,jj)%cp(1:2,kk)=alpha*cprjx(ii,jj)%cp(1:2,kk)
       end do
     end do
   end do
   if (ncpgrx>0) then
     do jj=1,n2dimx
       do ii=1,n1dimx
         nlmn=cprjx(ii,jj)%nlmn
         do kk=1,nlmn
           cprjy(ii,jj)%dcp(1:2,1:ncpgrx,kk)=alpha*cprjx(ii,jj)%dcp(1:2,1:ncpgrx,kk)
         end do
       end do
     end do
   end if
 else  ! alpha/=0 and beta/=0
   do jj=1,n2dimx
     do ii=1,n1dimx
       nlmn=cprjx(ii,jj)%nlmn
       cprjy(ii,jj)%nlmn=nlmn
       do kk=1,nlmn
         cprjy(ii,jj)%cp(1:2,kk)=alpha*cprjx(ii,jj)%cp(1:2,kk) &
&         +beta *cprjy(ii,jj)%cp(1:2,kk)
       end do
     end do
   end do
   if (ncpgrx>0) then
     do jj=1,n2dimx
       do ii=1,n1dimx
         nlmn=cprjx(ii,jj)%nlmn
         do kk=1,nlmn
           cprjy(ii,jj)%dcp(1:2,1:ncpgrx,kk)=alpha*cprjx(ii,jj)%dcp(1:2,1:ncpgrx,kk) &
&           +beta *cprjy(ii,jj)%dcp(1:2,1:ncpgrx,kk)
         end do
       end do
     end do
   end if
 end if

end subroutine pawcprj_axpby
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_zaxpby
!! NAME
!! pawcprj_zaxpby
!!
!! FUNCTION
!! Apply ZAXPBY (blas-like) operation with 2 cprj datastructures:
!!  cprjy(:,:) <- alpha.cprjx(:,:)+beta.cprjy(:,:)
!!  alpha and beta are COMPLEX scalars
!!
!! INPUTS
!!  alpha(2),beta(2)= alpha,beta COMPLEX factors
!!  cprjx(:,:) <type(pawcprj_type)>= input cprjx datastructure
!!
!! SIDE EFFECTS
!!  cprjy(:,:) <type(pawcprj_type)>= input/output cprjy datastructure
!!
!! PARENTS
!!      corrmetalwf1,extrapwf
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

 subroutine pawcprj_zaxpby(alpha,beta,cprjx,cprjy)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: alpha(2),beta(2)
!arrays
 type(pawcprj_type),intent(in) :: cprjx(:,:)
 type(pawcprj_type),intent(inout) :: cprjy(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,kk,ll,n1dimx,n1dimy,n2dimx,n2dimy,ncpgrx,ncpgry,nlmn
 real(dp) :: cp1,cp2,norma,normb
 character(len=500) :: msg

! *************************************************************************

 norma=alpha(1)**2+alpha(2)**2
 normb=beta(1) **2+beta(2) **2
 n1dimy=size(cprjy,dim=1);n2dimy=size(cprjy,dim=2);ncpgry=cprjy(1,1)%ncpgr
 if (norma>tol16) then
   n1dimx=size(cprjx,dim=1);n2dimx=size(cprjx,dim=2);ncpgrx=cprjx(1,1)%ncpgr
   msg = ""
   if (n1dimx/=n1dimy) msg = TRIM(msg)//"Error in pawcprj_zaxpby: n1 wrong sizes !"//ch10
   if (n2dimx/=n2dimy) msg = TRIM(msg)//"Error in pawcprj_zaxpby: n2 wrong sizes !"//ch10
   if (ncpgrx/=ncpgry) msg = TRIM(msg)//"Error in pawcprj_zaxpby: ncpgr wrong sizes !"//ch10
   if (LEN_TRIM(msg) > 0) then
     MSG_ERROR(msg)
   end if
 end if

 if (norma<=tol16) then
   do jj=1,n2dimy
     do ii=1,n1dimy
       nlmn=cprjy(ii,jj)%nlmn
       do kk=1,nlmn
         cp1=beta(1)*cprjy(ii,jj)%cp(1,kk)-beta(2)*cprjy(ii,jj)%cp(2,kk)
         cp2=beta(1)*cprjy(ii,jj)%cp(2,kk)+beta(2)*cprjy(ii,jj)%cp(1,kk)
         cprjy(ii,jj)%cp(1,kk)=cp1
         cprjy(ii,jj)%cp(2,kk)=cp2
       end do
     end do
   end do
   if (ncpgry>0) then
     do jj=1,n2dimy
       do ii=1,n1dimy
         nlmn=cprjy(ii,jj)%nlmn
         do kk=1,nlmn
           do ll=1,ncpgry
             cp1=beta(1)*cprjy(ii,jj)%dcp(1,ll,kk)-beta(2)*cprjy(ii,jj)%dcp(2,ll,kk)
             cp2=beta(1)*cprjy(ii,jj)%dcp(2,ll,kk)+beta(2)*cprjy(ii,jj)%dcp(1,ll,kk)
             cprjy(ii,jj)%dcp(1,ll,kk)=cp1
             cprjy(ii,jj)%dcp(2,ll,kk)=cp2
           end do
         end do
       end do
     end do
   end if
 else if (normb<=tol16) then
   do jj=1,n2dimx
     do ii=1,n1dimx
       nlmn=cprjx(ii,jj)%nlmn
       cprjy(ii,jj)%nlmn=nlmn
       do kk=1,nlmn
         cprjy(ii,jj)%cp(1,kk)=alpha(1)*cprjx(ii,jj)%cp(1,kk)-alpha(2)*cprjx(ii,jj)%cp(2,kk)
         cprjy(ii,jj)%cp(2,kk)=alpha(1)*cprjx(ii,jj)%cp(2,kk)+alpha(2)*cprjx(ii,jj)%cp(1,kk)
       end do
     end do
   end do
   if (ncpgrx>0) then
     do jj=1,n2dimx
       do ii=1,n1dimx
         nlmn=cprjx(ii,jj)%nlmn
         do kk=1,nlmn
           cprjy(ii,jj)%dcp(1,1:ncpgrx,kk)=alpha(1)*cprjx(ii,jj)%dcp(1,1:ncpgrx,kk) &
&           -alpha(2)*cprjx(ii,jj)%dcp(2,1:ncpgrx,kk)
           cprjy(ii,jj)%dcp(2,1:ncpgrx,kk)=alpha(1)*cprjx(ii,jj)%dcp(2,1:ncpgrx,kk) &
&           +alpha(2)*cprjx(ii,jj)%dcp(1,1:ncpgrx,kk)
         end do
       end do
     end do
   end if
 else
   do jj=1,n2dimx
     do ii=1,n1dimx
       nlmn=cprjx(ii,jj)%nlmn
       cprjy(ii,jj)%nlmn =nlmn
       do kk=1,nlmn
         cp1=alpha(1)*cprjx(ii,jj)%cp(1,kk)-alpha(2)*cprjx(ii,jj)%cp(2,kk) &
&         +beta(1) *cprjy(ii,jj)%cp(1,kk)-beta(2) *cprjy(ii,jj)%cp(2,kk)
         cp2=alpha(1)*cprjx(ii,jj)%cp(2,kk)+alpha(2)*cprjx(ii,jj)%cp(1,kk) &
&         +beta(1) *cprjy(ii,jj)%cp(2,kk)+beta(2) *cprjy(ii,jj)%cp(1,kk)
         cprjy(ii,jj)%cp(1,kk)=cp1
         cprjy(ii,jj)%cp(2,kk)=cp2
       end do
     end do
   end do
   if (ncpgrx>0) then
     do jj=1,n2dimx
       do ii=1,n1dimx
         nlmn=cprjx(ii,jj)%nlmn
         do kk=1,nlmn
           do ll=1,ncpgrx
             cp1=alpha(1)*cprjx(ii,jj)%dcp(1,ll,kk)-alpha(2)*cprjx(ii,jj)%dcp(2,ll,kk) &
&             +beta(1) *cprjy(ii,jj)%dcp(1,ll,kk)-beta(2) *cprjy(ii,jj)%dcp(2,ll,kk)
             cp2=alpha(1)*cprjx(ii,jj)%dcp(2,ll,kk)+alpha(2)*cprjx(ii,jj)%dcp(1,ll,kk) &
&             +beta(1) *cprjy(ii,jj)%dcp(2,ll,kk)+beta(2) *cprjy(ii,jj)%dcp(1,ll,kk)
             cprjy(ii,jj)%dcp(1,ll,kk)=cp1
             cprjy(ii,jj)%dcp(2,ll,kk)=cp2
           end do
         end do
       end do
     end do
   end if
 end if

end subroutine pawcprj_zaxpby
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_symkn
!! NAME
!! pawcprj_symkn
!!
!! FUNCTION
!! compute cprj for a given band and k point based on cprj at a symmetry-related
!! k point.
!!
!! INPUTS
!!  cprj_ikn (pawcprj_type) :: cprj for a single band and k point, typically a k point in the IBZ
!!  cprj_sym(4,nsym,natom) :: 1:3 shift, and 4 final atom, of symmetry isym operating on iatom
!!                            (S^{-1}(R - t) = r0 + L, see symatm.F90
!!  dimlmn(natom) :: ln dimension of each atom
!!  iband :: number of bands to treat, use -1 to treat all nband bands
!!  indlmn(6,lmnmax,ntypat) :: n,l,m dimensions for each atom type (see psps type)
!!  isym :: symmetry element used in current application
!!  itim :: 1 if time reversal also used, 0 else
!!  kpt(3) :: kpt vector used
!!  lmax :: max l value
!!  lmnmax :: max lmn value
!!  mband :: maximum number of bands
!!  natom :: number of atoms in cell
!!  nband :: number of bands in cprj_ikn
!!  nspinor :: number of spinors
!!  nsym :: total number of symmetry elements
!!  ntypat :: number of types of atoms
!!  typat(natom) :: type of each atom
!!  zarot(2*lmax+1,2*lmax+1,lmax+1,nsym) :: elements of rotation matrix for angular momentum states
!!                                          and symmetry operations. See m_paw_sphharm/setsym_ylm.
!!
!! OUTPUT
!!  cprj_fkn (pawcprj_type) :: cprj for a single band and k point where the k point is related to
!!    the input k point by a symmetry operation
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  This routine is based on M. Giantomassi's doctoral dissertation, formula 7.77. It is not clear
!!  whether it is implemented correctly for nonsymmorphic symmetries.
!!
!! PARENTS
!!      berryphase_new,cgwf,m_fock,make_grad_berry
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

 subroutine pawcprj_symkn(cprj_fkn,cprj_ikn,cprj_sym,dimlmn,iband,indlmn,&
&                       isym,itim,kpt,lmax,lmnmax,mband,natom,nband,nspinor,nsym,ntypat,&
&                       typat,zarot)

!Arguments---------------------------
!scalars
 integer,intent(in) :: iband,isym,itim,lmax,lmnmax,mband
 integer,intent(in) :: natom,nband,nspinor,nsym,ntypat

!arrays
 integer,intent(in) :: cprj_sym(4,nsym,natom),dimlmn(natom)
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),typat(natom)
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(in) :: zarot(2*lmax+1,2*lmax+1,lmax+1,nsym)
 type(pawcprj_type),intent(in) :: cprj_ikn(natom,mband*nspinor)
 type(pawcprj_type),intent(inout) :: cprj_fkn(natom,mband*nspinor) !vz_i

!Local variables---------------------------
!scalars
 integer :: iatm,iatom, ibct, ibnd, ibsp, ibst, icpgr, iin, il, il0, im
 integer :: ilmn, iln, iln0, ilpm, indexi, ispinor, itypat, jatm,jatom, mm, nlmn
 real(dp) :: kdotL, phr, phi
!arrays
 real(dp) :: rl(3), t1(2), t2(2)

! *************************************************************************

 if (iband == -1) then
   ibst = 1
   ibnd = nband
 else
   ibst = iband
   ibnd = iband
 end if

 do iatom = 1, natom
   iatm=iatom
   itypat = typat(iatom)
   nlmn = dimlmn(iatm)
   jatom = cprj_sym(4,isym,iatom)
   jatm=jatom
   rl(:) = cprj_sym(1:3,isym,iatom)
   kdotL = dot_product(rl,kpt)
   phr = cos(two_pi*kdotL)
   phi = sin(two_pi*kdotL)

   il0 = -1; iln0 = -1; indexi = 1
   do ilmn = 1, nlmn

     il = indlmn(1,ilmn,itypat)
     im = indlmn(2,ilmn,itypat)
     iin = indlmn(3,ilmn,itypat)
     iln = indlmn(5,ilmn,itypat)
     ilpm = 1 + il + im
     if (iln /= iln0) indexi = indexi + 2*il0 + 1

     do ibct = ibst, ibnd

       do ispinor = 1, nspinor

         ibsp = nspinor*(ibct-1) + ispinor

         t1(:) = zero
         do mm = 1, 2*il+1
           t1(1) = t1(1) + zarot(mm,ilpm,il+1,isym)*cprj_ikn(jatm,ibsp)%cp(1,indexi+mm)
           t1(2) = t1(2) + zarot(mm,ilpm,il+1,isym)*cprj_ikn(jatm,ibsp)%cp(2,indexi+mm)
         end do
         t2(1) = t1(1)*phr - t1(2)*phi
         t2(2) = t1(2)*phr + t1(1)*phi

         if (itim == 1) t2(2) = -t2(2)

         cprj_fkn(iatm,ibsp)%cp(1,ilmn) = t2(1)
         cprj_fkn(iatm,ibsp)%cp(2,ilmn) = t2(2)

! do same transformations for gradients of cprj_ikn
! note that ncpgr = 0 if no gradients present so this loop will not be executed
! in this case

         do icpgr = 1, cprj_ikn(jatom,ibsp)%ncpgr
           t1(:) = zero

           do mm = 1, 2*il+1
             t1(1) = t1(1) + zarot(mm,ilpm,il+1,isym)*cprj_ikn(jatm,ibsp)%dcp(1,icpgr,indexi+mm)
             t1(2) = t1(2) + zarot(mm,ilpm,il+1,isym)*cprj_ikn(jatm,ibsp)%dcp(2,icpgr,indexi+mm)
           end do

           t2(1) = t1(1)*phr - t1(2)*phi
           t2(2) = t1(2)*phr + t1(1)*phi

           if (itim == 1) t2(2) = -t2(2)

           cprj_fkn(iatm,ibsp)%dcp(1,icpgr,ilmn) = t2(1)
           cprj_fkn(iatm,ibsp)%dcp(2,icpgr,ilmn) = t2(2)

         end do ! end loop over ncpgr

       end do ! end loop over nspinor

     end do ! end loop over bands

     il0 = il; iln0 = iln
   end do ! end loop over ilmn
 end do ! end loop over atoms

 end subroutine pawcprj_symkn
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_conjg
!! NAME
!! pawcprj_conjg
!!
!! FUNCTION
!! conjugate a cprj datastructures:
!!  cprj(:,:) <- conjugate(cprj(:,:))
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  cprj(:,:) <type(pawcprj_type)>= input/output cprj datastructure
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

 subroutine pawcprj_conjg(cprj)

!Arguments ------------------------------------
!scalars
!arrays
 type(pawcprj_type),intent(inout) :: cprj(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,jj,kk,n1dim,n2dim,ncpgr,nlmn
!arrays

! *************************************************************************


 n1dim=size(cprj,dim=1);n2dim=size(cprj,dim=2);ncpgr=cprj(1,1)%ncpgr

 do jj=1,n2dim
   do ii=1,n1dim
     nlmn=cprj(ii,jj)%nlmn
     do kk=1,nlmn
       cprj(ii,jj)%cp(2,kk)=-cprj(ii,jj)%cp(2,kk)
     end do
   end do
 end do
 if (ncpgr>0) then
   do jj=1,n2dim
     do ii=1,n1dim
       nlmn=cprj(ii,jj)%nlmn
       do kk=1,nlmn
         cprj(ii,jj)%dcp(2,1:ncpgr,kk)=-cprj(ii,jj)%dcp(2,1:ncpgr,kk)
       end do
     end do
   end do
 end if

end subroutine pawcprj_conjg
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_lincom
!! NAME
!! pawcprj_lincom
!!
!! FUNCTION
!! Compute a LINear COMbination of cprj datastructure:
!!  cprj_out(:,:) <--- Sum_i [ alpha_i . cprj_i(:,:) ]
!!  alpha_i are COMPLEX scalars
!!
!! INPUTS
!!  alpha(2,nn)= alpha COMPLEX factors
!!  cprj_in(:,:) <type(pawcprj_type)>= input cprj_in datastructure
!!  nn= number of cprj involved in the linear combination
!!
!! OUTPUT
!!  cprj_out(:,:) <type(pawcprj_type)>= output cprj_out datastructure
!!
!! NOTES
!!  cprj_in and cprj_out must be dimensionned as cprj_in(n1,n2*nn) and cprj_in(n1,n2)
!!
!! PARENTS
!!      extrapwf,getdc1,lincom_cgcprj,wf_mixing
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

 subroutine pawcprj_lincom(alpha,cprj_in,cprj_out,nn)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn
 real(dp),intent(in) :: alpha(2,nn)
!arrays
 type(pawcprj_type),intent(in) :: cprj_in(:,:)
 type(pawcprj_type),intent(inout) :: cprj_out(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,in,jj,jn,kk,ll,n1in,n1out,n2in,n2out,ncpgrin,ncpgrout,nlmn
 real(dp) :: cp1,cp2
 character(len=500) :: msg

! *************************************************************************

 n1in=size(cprj_in,dim=1);n1out=size(cprj_out,dim=1)
 n2in=size(cprj_in,dim=2);n2out=size(cprj_out,dim=2)
 ncpgrin=cprj_in(1,1)%ncpgr;ncpgrout=cprj_out(1,1)%ncpgr

 msg = ""
 if (n1in/=n1out) msg = TRIM(msg)//"Bug in pawcprj_lincom: n1 wrong sizes!"//ch10
 if (n2in/=n2out*nn) msg = TRIM(msg)//"Bug in pawcprj_lincom: n2 wrong sizes!"//ch10
 if (ncpgrin/=ncpgrout) msg = TRIM(msg)//"Bug in pawcprj_lincom: ncpgr wrong sizes!"//ch10
 if (LEN_TRIM(msg) > 0) then
   MSG_ERROR(msg)
 end if

 do jj=1,n2out
   do ii=1,n1out
     nlmn=cprj_in(ii,jj)%nlmn
     cprj_out(ii,jj)%nlmn=nlmn
     cprj_out(ii,jj)%cp(1:2,1:nlmn)=zero
     jn=jj
     do in=1,nn
       do kk=1,nlmn
         cp1=cprj_out(ii,jj)%cp(1,kk) &
&         +alpha(1,in)*cprj_in(ii,jn)%cp(1,kk)-alpha(2,in)*cprj_in(ii,jn)%cp(2,kk)
         cp2=cprj_out(ii,jj)%cp(2,kk) &
&         +alpha(1,in)*cprj_in(ii,jn)%cp(2,kk)+alpha(2,in)*cprj_in(ii,jn)%cp(1,kk)
         cprj_out(ii,jj)%cp(1,kk)=cp1
         cprj_out(ii,jj)%cp(2,kk)=cp2
       end do
       jn=jn+n2out
     end do
   end do
 end do

 if (ncpgrin>0) then
   do jj=1,n2out
     do ii=1,n1out
       nlmn=cprj_in(ii,jj)%nlmn
       cprj_out(ii,jj)%dcp(1:2,1:ncpgrin,1:nlmn)=zero
       jn=jj
       do in=1,nn
         do kk=1,nlmn
           do ll=1,ncpgrin
             cp1=cprj_out(ii,jj)%dcp(1,ll,kk) &
&             +alpha(1,in)*cprj_in(ii,jn)%dcp(1,ll,kk) &
&             -alpha(2,in)*cprj_in(ii,jn)%dcp(2,ll,kk)
             cp2=cprj_out(ii,jj)%dcp(2,ll,kk) &
&             +alpha(1,in)*cprj_in(ii,jn)%dcp(2,ll,kk) &
             +alpha(2,in)*cprj_in(ii,jn)%dcp(1,ll,kk)
             cprj_out(ii,jj)%dcp(1,ll,kk)=cp1
             cprj_out(ii,jj)%dcp(2,ll,kk)=cp2
           end do
         end do
         jn=jn+n2out
       end do
     end do
   end do
 end if

end subroutine pawcprj_lincom
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_output
!! NAME
!! pawcprj_output
!!
!! FUNCTION
!! Output a cprj. Useful for debugging.
!!
!! INPUTS
!!  cprj(:,:) <type(pawcprj_type)>= cprj datastructure
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

 subroutine pawcprj_output(cprj)

!Arguments ------------------------------------
!scalar
!arrays
 type(pawcprj_type),intent(in) :: cprj(:,:)

!Local variables-------------------------------
!scalar
 integer :: ii,jj,kk,nlmn,n1dim,n2dim

! *************************************************************************

 n1dim=size(cprj,dim=1)
 n2dim=size(cprj,dim=2)

 write(std_out,'(a)')' pawcprj_output '

 do jj=1,n2dim
   do ii=1,n1dim
     write(std_out,'(a,i4,a,i4)')'atom ',ii,' band*k ',jj
     nlmn=cprj(ii,jj)%nlmn
     do kk=1,nlmn
       write(std_out,'(2f12.8)')cprj(ii,jj)%cp(1,kk),cprj(ii,jj)%cp(2,kk)
     end do
   end do
 end do

end subroutine pawcprj_output
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_get
!! NAME
!! pawcprj_get
!!
!! FUNCTION
!! Read the cprj_k for a given k-point from memory in cprj or from a temporary file
!!
!! INPUTS
!!  atind(natom)=index table for atoms (see iorder below)
!!  cprj(dimcp,nspinor*mband*mkmem*nsppol)=input cprj (used if mkmem/=0)
!!  dimcp=first dimension of cprj_k,cprj arrays (1 or natom)
!!  iband1=index of first band in cprj
!!  ibg=shift in cprj array to locate current k-point
!!  [icpgr]= (optional argument) if present, only component icpgr of
!!           input cprj gradient is copied into output cprj
!!           Not used if cprj(:,:)%ncpgr<icpgr (mkmem>0)
!!                    or ncpgr(optional)<icpgr (mkmem=0)
!!  ikpt=index of current k-point (only needed for the parallel distribution)
!!  iorder=0 if cprj ordering does not change during reading
!!         1 if cprj ordering changes during reading, depending on content of atind array:
!!              - if atind=atindx  (type-sorted=>unsorted)
!!              - if atind=atindx1 (unsorted=>type-sorted)
!!  isppol=index of current spin component
!!  mband=maximum number of bands
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!  [mpi_comm]=(optional argument) MPI communicator over (k-pts,bands,spins)
!!             Must be used in association with proc_distrb argument
!!  natom=number of atoms in cell
!!  nband=number of bands to import (usually 1 or nband_k)
!!  nband_k=total number of bands for this k-point
!!  [ncpgr]=(optional argument) second dimension of cprj%dcp(2,ncpgr,nlmn
!!          stored in memory (mkmem>0) or present on disk (mkmem=0))
!!          needed only when optional argument icpgr is present
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  [proc_distrb(nkpt,nband,nsppol)]=(optional argument) processor distribution
!!          Describe how cprj datastructures are distributed over processors
!!          When present, mpicomm argument must be also present
!!  uncp=unit number for cprj data (used if mkmem=0)
!!
!! OUTPUT
!!  cprj_k(dimcp,nspinor*nband) <type(pawcprj_type)>= output cprj datastructure
!!
!! PARENTS
!!      berryphase_new,cgwf,chern_number,datafordmft,dfpt_nstpaw,dfpt_vtowfk
!!      dfpt_wfkfermi,dotprod_set_cgcprj,dotprodm_sumdiag_cgcprj,extrapwf
!!      fock2ACE,forstrnps,m_plowannier,make_grad_berry,optics_paw
!!      optics_paw_core,pawmkrhoij,posdoppler,rf2_init,smatrix_pawinit
!!      suscep_stat,wf_mixing
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

 subroutine pawcprj_get(atind,cprj_k,cprj,dimcp,iband1,ibg,ikpt,iorder,isppol,mband,&
&                    mkmem,natom,nband,nband_k,nspinor,nsppol,uncp,&
&                    icpgr,ncpgr,mpicomm,proc_distrb) ! optionals arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dimcp,iband1,ibg,ikpt,iorder,isppol,mband,mkmem,natom
 integer,intent(in) :: nband,nband_k,nspinor,nsppol,uncp
 integer,intent(in),optional :: icpgr,mpicomm,ncpgr
!arrays
 integer,intent(in) :: atind(natom)
 integer,intent(in),optional :: proc_distrb(:,:,:)
 type(pawcprj_type),intent(in) :: cprj(dimcp,nspinor*mband*mkmem*nsppol)
 type(pawcprj_type),intent(inout) :: cprj_k(dimcp,nspinor*nband)

!Local variables-------------------------------
!scalars
 integer :: iatm,iatom,ib,ibsp,icpgr_,isp,ispinor,jband,me,nband0,ncpgr_
 logical :: has_distrb,has_icpgr
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: tmp(:,:,:)

! *************************************************************************

 ncpgr_=cprj_k(1,1)%ncpgr;if (present(ncpgr)) ncpgr_=ncpgr
 icpgr_=-1;if(present(icpgr)) icpgr_=icpgr
 has_icpgr=(icpgr_>0.and.icpgr_<=ncpgr_)
 if (present(icpgr).and.(.not.present(ncpgr))) then
   msg='ncpgr must be present when icpgr is present (pawcprj_get)!'
   MSG_BUG(msg)
 end if
 if (has_icpgr.and.cprj_k(1,1)%ncpgr<1) then
   msg='cprj_k%ncpgr not consistent with icpgr (pawcprj_get)!'
   MSG_BUG(msg)
 end if

!MPI data
 has_distrb=present(proc_distrb)
 if (has_distrb) then
   if (.not.present(mpicomm)) then
     msg='mpicomm must be present when proc_distrb is present (pawcprj_get)!'
     MSG_BUG(msg)
   end if
   me=xmpi_comm_rank(mpicomm)
 end if

 if (mkmem==0) then

   if (iband1==1) then
     read(uncp) nband0
     if (nband_k/=nband0) then
       msg='_PAW file was not created with the right options (pawcprj_get)!'
       MSG_BUG(msg)
     end if
   end if

   isp=0;jband=iband1-1
   do ib=1,nband
     jband=jband+1
     if (has_distrb) then
       if (abs(proc_distrb(ikpt,jband,isppol)-me)/=0) then
         isp=isp+nspinor
         cycle
       end if
     end if
     do ispinor=1,nspinor
       isp=isp+1
       if (iorder==0) then
         if (ncpgr_==0) then
           do iatom=1,dimcp
             read(uncp) cprj_k(iatom,isp)%cp(:,:)
           end do
         else
           if (has_icpgr) then
             do iatom=1,dimcp
               LIBPAW_ALLOCATE(tmp,(2,ncpgr_,cprj_k(iatom,1)%nlmn))
               read(uncp) cprj_k(iatom,isp)%cp(:,:),tmp(:,:,:)
               cprj_k(iatom,isp)%dcp(:,1,:)=tmp(:,icpgr_,:)
               LIBPAW_DEALLOCATE(tmp)
             end do
           else
             do iatom=1,dimcp
               read(uncp) cprj_k(iatom,isp)%cp(:,:),cprj_k(iatom,isp)%dcp(:,:,:)
             end do
           end if
         end if
       else
         if (ncpgr_==0) then
           do iatom=1,dimcp
             iatm=min(atind(iatom),dimcp)
             read(uncp) cprj_k(iatm,isp)%cp(:,:)
           end do
         else
           if (has_icpgr) then
             do iatom=1,dimcp
               iatm=min(atind(iatom),dimcp)
               LIBPAW_ALLOCATE(tmp,(2,ncpgr_,cprj_k(iatm,1)%nlmn))
               read(uncp) cprj_k(iatm,isp)%cp(:,:),tmp(:,:,:)
               cprj_k(iatm,isp)%dcp(:,1,:)=tmp(:,icpgr_,:)
               LIBPAW_DEALLOCATE(tmp)
             end do
           else
             do iatom=1,dimcp
               iatm=min(atind(iatom),dimcp)
               read(uncp) cprj_k(iatm,isp)%cp(:,:),cprj_k(iatm,isp)%dcp(:,:,:)
             end do
           end if
         end if
       end if
     end do
   end do

 else

   isp=0;ibsp=ibg+nspinor*(iband1-1);jband=iband1-1
   do ib=1,nband
     jband=jband+1
     if (has_distrb) then
       if (abs(proc_distrb(ikpt,jband,isppol)-me)/=0) then
         isp=isp+nspinor;ibsp=ibsp+nspinor
         cycle
       end if
     end if
     do ispinor=1,nspinor
       isp=isp+1;ibsp=ibsp+1
       if (iorder==0) then
         if (ncpgr_==0) then
           do iatom=1,dimcp
             cprj_k(iatom,isp)%cp(:,:)=cprj(iatom,ibsp)%cp(:,:)
           end do
         else
           if (has_icpgr) then
             do iatom=1,dimcp
               cprj_k(iatom,isp)%cp(:,:)   =cprj(iatom,ibsp)%cp(:,:)
               cprj_k(iatom,isp)%dcp(:,1,:)=cprj(iatom,ibsp)%dcp(:,icpgr_,:)
             end do
           else
             do iatom=1,dimcp
               cprj_k(iatom,isp)%cp(:,:)   =cprj(iatom,ibsp)%cp(:,:)
               cprj_k(iatom,isp)%dcp(:,:,:)=cprj(iatom,ibsp)%dcp(:,:,:)
             end do
           end if
         end if
       else
         if (ncpgr_==0) then
           do iatom=1,dimcp
             iatm=min(atind(iatom),dimcp)
             cprj_k(iatm,isp)%cp(:,:)=cprj(iatom,ibsp)%cp(:,:)
           end do
         else
           if (has_icpgr) then
             do iatom=1,dimcp
               iatm=min(atind(iatom),dimcp)
               cprj_k(iatm,isp)%cp(:,:)   =cprj(iatom,ibsp)%cp(:,:)
               cprj_k(iatm,isp)%dcp(:,1,:)=cprj(iatom,ibsp)%dcp(:,icpgr_,:)
             end do
           else
             do iatom=1,dimcp
               iatm=min(atind(iatom),dimcp)
               cprj_k(iatm,isp)%cp(:,:)   =cprj(iatom,ibsp)%cp(:,:)
               cprj_k(iatm,isp)%dcp(:,:,:)=cprj(iatom,ibsp)%dcp(:,:,:)
             end do
           end if
         end if
       end if
     end do
   end do

 end if

end subroutine pawcprj_get
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_put
!! NAME
!! pawcprj_put
!!
!! FUNCTION
!! Write cprj_k for a given set of (n,k) into memory in cprj, or into a temporary file
!!
!! INPUTS
!!  atind(natom)=index table for atoms (see iorder below)
!!  cprj_k(dimcp,nspinor*nband) <type(pawcprj_type)>= input cprj datastructure
!!  dimcp=first dimension of cprj_k,cprjnk arrays (1 or natom)
!!  iband1=index of first band in cprj
!!  ibg=shift in cprj array to locate current k-point
!!  ikpt=index of current k-point (only needed for the parallel distribution)
!!  iorder=0 if cprj ordering does not change during reading
!!         1 if cprj ordering changes during writing, depending on content of atind array:
!!              - if atind=atindx  (type-sorted->unsorted)
!!              - if atind=atindx1 (unsorted->type-sorted)
!!  isppol=index of current spin component
!!  mband=maximum number of bands
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!  [mpi_comm]=(optional argument) MPI communicator over (k-pts,bands,spins)
!!             Must be used in association with proc_distrb argument
!!  [mpi_comm_band]=(optional argument) MPI communicator over bands
!!             Must be used in association with proc_distrb argument
!!  natom=number of atoms in cell
!!  nband=number of bands to export (usually 1, nband_k or nblockbd)
!!  nband_k=total number of bands for this k-point
!!  nlmn(dimcp)=array of dimensions of cprj_k,cprjnk datastructures
!!  nspinor=number of spinorial components of the wavefunctions (on current proc)
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  [proc_distrb(nkpt,nband,nsppol)]=(optional argument) processor distribution
!!          Describe how cprj datastructures are distributed over processors
!!          When present, mpicomm argument must be also present
!!  [to_be_gathered]=(optional argument) TRUE if cprj_k arrays have to be
!!                   gathered between procs (band-fft parallelism only)
!!  uncp=unit number for cprj data (used if mkmem=0)
!!
!! SIDE EFFECTS
!!  cprj(dimcp,nspinor*mband*mkmem*nsppol)=output cprj (used if mkmem/=0)
!!
!! PARENTS
!!      berryphase_new,cgwf,ctocprj,dfpt_vtowfk,extrapwf,vtowfk,wf_mixing
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

 subroutine pawcprj_put(atind,cprj_k,cprj,dimcp,iband1,ibg,ikpt,iorder,isppol,mband,&
&           mkmem,natom,nband,nband_k,nlmn,nspinor,nsppol,uncp,&
&           mpicomm,mpi_comm_band,proc_distrb,to_be_gathered) ! Optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband1,ibg,ikpt,iorder,isppol,dimcp,mband,mkmem
 integer,intent(in) :: natom,nband,nband_k,nspinor,nsppol,uncp
 integer,intent(in),optional :: mpicomm,mpi_comm_band
 logical,optional,intent(in) :: to_be_gathered
!arrays
 integer,intent(in) :: atind(natom),nlmn(dimcp)
 integer,intent(in),optional :: proc_distrb(:,:,:)
 type(pawcprj_type),intent(inout) :: cprj(dimcp,nspinor*mband*mkmem*nsppol)
 type(pawcprj_type),intent(in) :: cprj_k(dimcp,nspinor*nband)

!Local variables-------------------------------
!scalars
 integer :: iatm,iatom,iband,ibsp,icpgr,ierr,ii,ilmn,isp,ispinor,jband,jj
 integer :: lmndim,me,ncpgr,nproc_band
 logical :: has_distrb,to_be_gathered_
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: buffer1(:),buffer2(:)
! *************************************************************************

 ncpgr=cprj_k(1,1)%ncpgr
 to_be_gathered_=.false.;if (present(to_be_gathered)) to_be_gathered_=to_be_gathered

!MPI data
 nproc_band=1;if (present(mpi_comm_band)) nproc_band=xmpi_comm_size(mpi_comm_band)
 has_distrb=present(proc_distrb)
 if (has_distrb) then
   if (.not.present(mpicomm)) then
     msg='mpicomm must be present when proc_distrb is present (pawcprj_put)!'
     MSG_BUG(msg)
   end if
   me=xmpi_comm_rank(mpicomm)
 end if

 if (nproc_band==1.or.(.not.to_be_gathered_)) then

   if (mkmem==0) then

     if (iband1==1) write(uncp) nband_k

     isp=0;jband=iband1-1
     do iband=1,nband
       jband=jband+1
       if (has_distrb) then
         if (abs(proc_distrb(ikpt,jband,isppol)-me)/=0) then
           isp=isp+nspinor
           cycle
         end if
       end if
       do ispinor=1,nspinor
         isp=isp+1
         if (iorder==0) then
           do iatom=1,dimcp
             if (ncpgr==0) then
               write(uncp) cprj_k(iatom,isp)%cp(:,:)
             else
               write(uncp) cprj_k(iatom,isp)%cp(:,:),cprj_k(iatom,isp)%dcp(:,:,:)
             end if
           end do
         else
           do iatom=1,dimcp
             iatm=min(atind(iatom),dimcp)
             if (ncpgr==0) then
               write(uncp) cprj_k(iatm,isp)%cp(:,:)
             else
               write(uncp) cprj_k(iatm,isp)%cp(:,:),cprj_k(iatm,isp)%dcp(:,:,:)
             end if
           end do
         end if
       end do
     end do

   else

     isp=0;ibsp=ibg+nspinor*(iband1-1);jband=iband1-1
     do iband=1,nband
       jband=jband+1
       if (has_distrb) then
         if (abs(proc_distrb(ikpt,jband,isppol)-me)/=0) then
           isp=isp+nspinor;ibsp=ibsp+nspinor
           cycle
         end if
       end if
       do ispinor=1,nspinor
         isp=isp+1;ibsp=ibsp+1
         if (iorder==0) then
           do iatom=1,dimcp
             cprj(iatom,ibsp)%cp(:,:)=cprj_k(iatom,isp)%cp(:,:)
             if (ncpgr>0) cprj(iatom,ibsp)%dcp(:,:,:)=cprj_k(iatom,isp)%dcp(:,:,:)
           end do
         else
           do iatom=1,dimcp
             iatm=min(atind(iatom),dimcp)
             cprj(iatom,ibsp)%cp(:,:)=cprj_k(iatm,isp)%cp(:,:)
             if (ncpgr>0) cprj(iatom,ibsp)%dcp(:,:,:)=cprj_k(iatm,isp)%dcp(:,:,:)
           end do
         end if
       end do
     end do

   end if

 else ! np_band>1

   lmndim=2*sum(nlmn(1:dimcp))*(1+ncpgr)*nspinor
   LIBPAW_ALLOCATE(buffer1,(lmndim))
   LIBPAW_ALLOCATE(buffer2,(lmndim*nproc_band))
   isp=0;ibsp=ibg+nspinor*(iband1-1)
   do iband=1,nband  ! must be nblockbd for band-fft parallelism
     jj=1
     do ispinor=1,nspinor
       isp=isp+1
       do iatom=1,dimcp
         if (iorder==0) then
           iatm=iatom
         else
           iatm=min(atind(iatom),dimcp)
         end if
         do ilmn=1,nlmn(iatm)
           buffer1(jj:jj+1)=cprj_k(iatm,isp)%cp(1:2,ilmn)
           jj=jj+2
         end do
         if (ncpgr>0) then
           do ilmn=1,nlmn(iatm)
             do icpgr=1,ncpgr
               buffer1(jj:jj+1)=cprj_k(iatm,isp)%dcp(1:2,icpgr,ilmn)
               jj=jj+2
             end do
           end do
         end if
       end do !iatom
     end do !ispinor
     call xmpi_allgather(buffer1,lmndim,buffer2,mpi_comm_band,ierr)
     jj=1
     do ii=1,nproc_band
       do ispinor=1,nspinor
         ibsp=ibsp+1
         do iatom=1,dimcp
           if (iorder==0) then
             iatm=iatom
           else
             iatm=min(atind(iatom),dimcp)
           end if
           do ilmn=1,nlmn(iatm)
             cprj(iatom,ibsp)%cp(1:2,ilmn)=buffer2(jj:jj+1)
             jj=jj+2
           end do
           if (ncpgr>0) then
             do ilmn=1,nlmn(iatm)
               do icpgr=1,ncpgr
                 cprj(iatom,ibsp)%dcp(1:2,icpgr,ilmn)=buffer2(jj:jj+1)
                 jj=jj+2
               end do
             end do
           end if
         end do !iatom
       end do !ispinor
     end do !ii=1,nproc_band
   end do !iband
   LIBPAW_DEALLOCATE(buffer1)
   LIBPAW_DEALLOCATE(buffer2)

 end if ! mode_para=b, nband

end subroutine pawcprj_put
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_reorder
!! NAME
!! pawcprj_reorder
!!
!! FUNCTION
!! Change the order of a cprj datastructure
!!   From unsorted cprj to atom-sorted cprj (atm_indx=atindx)
!!   From atom-sorted cprj to unsorted cprj (atm_indx=atindx1)
!!
!! INPUTS
!!  atm_indx(natom)=index table for atoms
!!   From unsorted cprj to atom-sorted cprj (atm_indx=atindx)
!!   From atom-sorted cprj to unsorted cprj (atm_indx=atindx1)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  cprj(:,:) <type(pawcprj_type)>= cprj datastructure
!!
!! PARENTS
!!      fock2ACE,forstrnps,ks_ddiago,scfcv
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

 subroutine pawcprj_reorder(cprj,atm_indx)

!Arguments ------------------------------------
!scalars
!arrays
 integer,intent(in) :: atm_indx(:)
 type(pawcprj_type),intent(inout) :: cprj(:,:)

!Local variables-------------------------------
!scalars
 integer :: iexit,ii,jj,kk,n1atindx,n1cprj,n2cprj,ncpgr
 character(len=100) :: msg
!arrays
 integer,allocatable :: nlmn(:)
 type(pawcprj_type),allocatable :: cprj_tmp(:,:)

! *************************************************************************

 n1cprj=size(cprj,dim=1);n2cprj=size(cprj,dim=2)
 n1atindx=size(atm_indx,dim=1)
 if (n1cprj==0.or.n2cprj==0.or.n1atindx<=1) return

 if (n1cprj/=n1atindx) then
   msg='wrong sizes (pawcprj_reorder)!'
   MSG_BUG(msg)
 end if

!Nothing to do when the atoms are already sorted
 iexit=1;ii=0
 do while (iexit==1.and.ii<n1atindx)
   ii=ii+1
   if (atm_indx(ii)/=ii) iexit=0
 end do
 if (iexit==1) return

 LIBPAW_ALLOCATE(nlmn,(n1cprj))
 do ii=1,n1cprj
   nlmn(ii)=cprj(ii,1)%nlmn
 end do
 ncpgr=cprj(1,1)%ncpgr
 LIBPAW_DATATYPE_ALLOCATE(cprj_tmp,(n1cprj,n2cprj))
 call pawcprj_alloc(cprj_tmp,ncpgr,nlmn)
 call pawcprj_copy(cprj,cprj_tmp)
 call pawcprj_free(cprj)

 do jj=1,n2cprj
   do ii=1,n1cprj
     kk=atm_indx(ii)
     cprj(kk,jj)%nlmn=nlmn(ii)
     cprj(kk,jj)%ncpgr=ncpgr
     LIBPAW_ALLOCATE(cprj(kk,jj)%cp,(2,nlmn(ii)))
     cprj(kk,jj)%cp(:,:)=cprj_tmp(ii,jj)%cp(:,:)
     if (ncpgr>0) then
       LIBPAW_ALLOCATE(cprj(kk,jj)%dcp,(2,ncpgr,nlmn(ii)))
       cprj(kk,jj)%dcp(:,:,:)=cprj_tmp(ii,jj)%dcp(:,:,:)
     end if
   end do
 end do

 call pawcprj_free(cprj_tmp)
 LIBPAW_DATATYPE_DEALLOCATE(cprj_tmp)
 LIBPAW_DEALLOCATE(nlmn)

end subroutine pawcprj_reorder
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_mpi_exch
!! NAME
!! pawcprj_mpi_exch
!!
!! FUNCTION
!! Exchange a pawcprj_type between two processors inside a MPI communicator.
!!
!! INPUTS
!!  natom=Number of atoms (size of first dimension of Cprj_send and Cprj_recv).
!!  n2dim=Size of the second dimension.
!!  nlmn(natom)=Number of nlm partial waves for each atom.
!!  Cprj_send= The datatype to be transmitted.
!!  receiver=ID of the receiver in spaceComm.
!!  sender=ID of the sender in spaceComm.
!!  spaceComm=MPI Communicator.
!!
!! OUTPUT
!!  ierr=Error status.
!!  Cprj_recv=The datatype copied on proc. receiver.
!!
!! NOTES
!!  If sender==receiver, Cprj_send is copied into Cprj_recv.
!!  It should be easy to avoid this additional copy in the calling routine.
!!
!! PARENTS
!!      outkss
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine pawcprj_mpi_exch(natom,n2dim,nlmn,ncpgr,Cprj_send,Cprj_recv,sender,receiver,spaceComm,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,n2dim,ncpgr
 integer,intent(in) :: sender,receiver,spaceComm
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: nlmn(natom)
 type(pawcprj_type),intent(in) :: Cprj_send(:,:)
 type(pawcprj_type),intent(inout) :: Cprj_recv(:,:)

!Local variables-------------------------------
!scalars
 integer :: iat,jj,t2dim,tcpgr,n1dim,nn
 integer :: ntotcp,ipck,rank
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: buffer_cp(:,:),buffer_cpgr(:,:,:)

! *************************************************************************

 n1dim=0
 t2dim=0
 tcpgr=0
 ierr=0
 if (sender==receiver) then
   call pawcprj_copy(Cprj_send,Cprj_recv)
   return
 end if

 rank = xmpi_comm_rank(spaceComm)

 nn=size(nlmn,dim=1)
 if (rank==sender) then
   n1dim=size(Cprj_send,dim=1)
   t2dim=size(Cprj_send,dim=2)
   tcpgr=Cprj_send(1,1)%ncpgr
 end if
 if (rank==receiver) then
   n1dim=size(Cprj_recv,dim=1)
   t2dim=size(Cprj_recv,dim=2)
   tcpgr=Cprj_recv(1,1)%ncpgr
 end if
 if (rank/=sender.and.rank/=receiver) then
   write(msg,'(a,3i0)') &
&   'rank is not equal to sender or receiver (pawcprj_mpi_exch): ',rank, sender, receiver
   MSG_BUG(msg)
 end if

 ntotcp=n2dim*SUM(nlmn(:))

 LIBPAW_ALLOCATE(buffer_cp,(2,ntotcp))
 if (ncpgr/=0)  then
   LIBPAW_ALLOCATE(buffer_cpgr,(2,ncpgr,ntotcp))
 end if

!=== Pack Cprj_send ===
 if (rank==sender) then
   ipck=0
   do jj=1,n2dim
     do iat=1,natom
       nn=nlmn(iat)
       buffer_cp(:,ipck+1:ipck+nn)=Cprj_send(iat,jj)%cp(:,1:nn)
       if (ncpgr/=0) buffer_cpgr(:,:,ipck+1:ipck+nn)=Cprj_send(iat,jj)%dcp(:,:,1:nn)
       ipck=ipck+nn
     end do
   end do
 end if

!=== Transmit data ===
 call xmpi_exch(buffer_cp,2*ntotcp,sender,buffer_cp,receiver,spaceComm,ierr)
 if (ncpgr/=0) then
   call xmpi_exch(buffer_cpgr,2*ncpgr*ntotcp,sender,buffer_cpgr,receiver,spaceComm,ierr)
 end if

!=== UnPack buffers into Cprj_recv ===
 if (rank==receiver) then
   ipck=0
   do jj=1,n2dim
     do iat=1,natom
       nn=nlmn(iat)
       Cprj_recv(iat,jj)%cp(:,1:nn)=buffer_cp(:,ipck+1:ipck+nn)
       if (ncpgr/=0) Cprj_recv(iat,jj)%dcp(:,:,1:nn)=buffer_cpgr(:,:,ipck+1:ipck+nn)
       ipck=ipck+nn
     end do
   end do
 end if

 LIBPAW_DEALLOCATE(buffer_cp)
 if (ncpgr/=0)  then
   LIBPAW_DEALLOCATE(buffer_cpgr)
 end if

end subroutine pawcprj_mpi_exch
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_mpi_send
!! NAME
!! pawcprj_mpi_send
!!
!! FUNCTION
!! Send a pawcprj_type inside a MPI communicator.
!!
!! INPUTS
!!  natom=Number of atoms (size of first dimension of cprj_out).
!!  n2dim=Size of the second dimension.
!!  nlmn(natom)=Number of nlm partial waves for each atom.
!!  ncpgr = number of gradients in cprj_out
!!  cprj_out= The datatype to be transmitted.
!!  receiver=ID of the receiver in spaceComm.
!!  spaceComm=MPI Communicator.
!!
!! OUTPUT
!!  ierr=Error status.
!!
!! NOTES
!!   perhaps in general it is more efficient to use pawcprj_mpi_exch but it is
!!   convenient for coding to have separate send and recieve routines.
!!
!! PARENTS
!!      berryphase_new,posdoppler
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine pawcprj_mpi_send(natom,n2dim,nlmn,ncpgr,cprj_out,receiver,spaceComm,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,n2dim,ncpgr
 integer,intent(in) :: receiver,spaceComm
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: nlmn(natom)
 type(pawcprj_type),intent(in) :: cprj_out(:,:)

!Local variables-------------------------------
!scalars
 integer :: iat,jj,t2dim,tcpgr,n1dim,nn
 integer :: ntotcp,ipck,tag
 character(len=100) :: msg
!arrays
 real(dp),allocatable :: buffer_cp(:,:),buffer_cpgr(:,:,:)

! *************************************************************************

 n1dim=0
 t2dim=0
 tcpgr=0
 ierr=0

 nn=size(nlmn,dim=1)
 n1dim=size(cprj_out,dim=1)
 t2dim=size(cprj_out,dim=2)
 tcpgr=cprj_out(1,1)%ncpgr

 if (nn/=n1dim) then
   msg='size mismatch in natom (pawcprj_mpi_send)!'
   MSG_BUG(msg)
 end if
 if (t2dim/=n2dim) then
   msg='size mismatch in dim=2 (pawcprj_mpi_send)!'
   MSG_BUG(msg)
 end if
 if (tcpgr/=ncpgr) then
   msg='size mismatch in ncpgr (pawcprj_mpi_send)!'
   MSG_BUG(msg)
 end if

 ntotcp=n2dim*SUM(nlmn(:))

 LIBPAW_ALLOCATE(buffer_cp,(2,ntotcp))
 if (ncpgr/=0)  then
   LIBPAW_ALLOCATE(buffer_cpgr,(2,ncpgr,ntotcp))
 end if

!=== Pack cprj_out ====
 ipck=0
 do jj=1,n2dim
   do iat=1,natom
     nn=nlmn(iat)
     buffer_cp(:,ipck+1:ipck+nn)=cprj_out(iat,jj)%cp(:,1:nn)
     if (ncpgr/=0) buffer_cpgr(:,:,ipck+1:ipck+nn)=cprj_out(iat,jj)%dcp(:,:,1:nn)
     ipck=ipck+nn
   end do
 end do

!=== Transmit data ===
 tag = 2*ntotcp
 call xmpi_send(buffer_cp,receiver,tag,spaceComm,ierr)
 if (ncpgr/=0) then
   tag=tag*ncpgr
   call xmpi_send(buffer_cpgr,receiver,tag,spaceComm,ierr)
 end if

!=== Clean up ===
 LIBPAW_DEALLOCATE(buffer_cp)
 if (ncpgr/=0)  then
   LIBPAW_DEALLOCATE(buffer_cpgr)
 end if

end subroutine pawcprj_mpi_send
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_mpi_recv
!! NAME
!! pawcprj_mpi_recv
!!
!! FUNCTION
!! Receive a pawcprj_type inside a MPI communicator.
!!
!! INPUTS
!!  natom=Number of atoms (size of first dimension of Cprj_in).
!!  n2dim=Size of the second dimension.
!!  nlmn(natom)=Number of nlm partial waves for each atom.
!!  ncpgr = number of gradients in cprj_in
!!  sender=ID of the sender in spaceComm.
!!  spaceComm=MPI Communicator.
!!
!! OUTPUT
!!  ierr=Error status.
!!  cprj_in=The datatype copied on proc. receiver.
!!
!! NOTES
!!   Perhaps in general it is more efficient to use pawcprj_mpi_exch but it is
!!   convenient for coding to have separate send and receive routines.
!!
!! PARENTS
!!      berryphase_new,posdoppler
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine pawcprj_mpi_recv(natom,n2dim,nlmn,ncpgr,cprj_in,sender,spaceComm,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,n2dim,ncpgr
 integer,intent(in) :: sender,spaceComm
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: nlmn(natom)
 type(pawcprj_type),intent(inout) :: cprj_in(:,:)

!Local variables-------------------------------
!scalars
 integer :: iat,jj,t2dim,tcpgr,n1dim,nn
 integer :: ntotcp,ipck,tag
 character(len=100) :: msg
!arrays
 real(dp),allocatable :: buffer_cp(:,:),buffer_cpgr(:,:,:)

! *************************************************************************

 n1dim=0
 t2dim=0
 tcpgr=0
 ierr=0

 nn=size(nlmn,dim=1)
 n1dim=size(cprj_in,dim=1)
 t2dim=size(cprj_in,dim=2)
 tcpgr=cprj_in(1,1)%ncpgr

 if (nn/=n1dim) then
   msg='size mismatch in natom (pawcprj_mpi_recv)!'
   MSG_BUG(msg)
 end if
 if (t2dim/=n2dim) then
   msg='size mismatch in dim=2 (pawcprj_mpi_recv)!'
   MSG_BUG(msg)
 end if
 if (tcpgr/=ncpgr) then
   msg='size mismatch in ncpgr (pawcprj_mpi_recv)!'
   MSG_BUG(msg)
 end if

 ntotcp=n2dim*SUM(nlmn(:))

 LIBPAW_ALLOCATE(buffer_cp,(2,ntotcp))
 if (ncpgr/=0)  then
   LIBPAW_ALLOCATE(buffer_cpgr,(2,ncpgr,ntotcp))
 end if

!=== Receive data ===
 tag = 2*ntotcp
 call xmpi_recv(buffer_cp,sender,tag,spaceComm,ierr)
 if (ncpgr/=0) then
   tag=tag*ncpgr
   call xmpi_recv(buffer_cpgr,sender,tag,spaceComm,ierr)
 end if

!=== UnPack buffers into cprj_in ===
 ipck=0
 do jj=1,n2dim
   do iat=1,natom
     nn=nlmn(iat)
     cprj_in(iat,jj)%cp(:,1:nn)=buffer_cp(:,ipck+1:ipck+nn)
     if (ncpgr/=0) cprj_in(iat,jj)%dcp(:,:,1:nn)=buffer_cpgr(:,:,ipck+1:ipck+nn)
     ipck=ipck+nn
   end do
 end do

!=== Clean up ===
 LIBPAW_DEALLOCATE(buffer_cp)
 if (ncpgr/=0)  then
   LIBPAW_DEALLOCATE(buffer_cpgr)
 end if

end subroutine pawcprj_mpi_recv
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_mpi_sum
!! NAME
!! pawcprj_mpi_sum
!!
!! FUNCTION
!! Perform MPI_SUM on a pawcprj_type inside a MPI communicator.
!!
!! INPUTS
!!  spaceComm=MPI Communicator.
!!
!! SIDE EFFECTS
!!  cprj=the cprj datastructure
!!  ierr=Error status.
!!
!! PARENTS
!!      ctocprj
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine pawcprj_mpi_sum(cprj,spaceComm,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spaceComm
 integer,intent(out) :: ierr
!arrays
 type(pawcprj_type),intent(inout) :: cprj(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: maxBytes=100*1024*1024 ! 100 MBytes
 integer :: ii,ipck,jj,ncpgr,nlmn,nn,n1dim,n2dim,n2dim1,n2dim2,sizeBytes,step
 logical,parameter :: save_memory=.true.
!arrays
 real(dp),allocatable :: buffer_cprj(:,:,:)

! *************************************************************************

 if (xmpi_comm_size(spaceComm)<2) return

 n1dim=size(cprj,1);n2dim=size(cprj,2)
 nlmn=sum(cprj(:,:)%nlmn)
 ncpgr=maxval(cprj(:,:)%ncpgr)

 step=n2dim
 if (save_memory) then
   sizeBytes=2*(1+ncpgr)*nlmn *8
   step=n2dim/max(1,sizeBytes/maxBytes)
   if (step==0) step=1
 end if

 do n2dim1=1,n2dim,step

   n2dim2=min(n2dim1+step-1,n2dim)
   nlmn=sum(cprj(:,n2dim1:n2dim2)%nlmn)
   LIBPAW_ALLOCATE(buffer_cprj,(2,1+ncpgr,nlmn))

   ipck=0 ; buffer_cprj=zero
   do jj=n2dim1,n2dim2
     do ii=1,n1dim
       nn=cprj(ii,jj)%nlmn
       buffer_cprj(:,1,ipck+1:ipck+nn)=cprj(ii,jj)%cp(:,1:nn)
       if (cprj(ii,jj)%ncpgr/=0) buffer_cprj(:,2:1+ncpgr,ipck+1:ipck+nn)=cprj(ii,jj)%dcp(:,1:ncpgr,1:nn)
       ipck=ipck+nn
     end do
   end do

   call xmpi_sum(buffer_cprj,spaceComm,ierr)

   ipck=0
   do jj=n2dim1,n2dim2
     do ii=1,n1dim
       nn=cprj(ii,jj)%nlmn
       cprj(ii,jj)%cp(:,1:nn)=buffer_cprj(:,1,ipck+1:ipck+nn)
       if (cprj(ii,jj)%ncpgr/=0) cprj(ii,jj)%dcp(:,1:ncpgr,1:nn)=buffer_cprj(:,2:1+ncpgr,ipck+1:ipck+nn)
       ipck=ipck+nn
     end do
   end do

   LIBPAW_DEALLOCATE(buffer_cprj)

 end do

end subroutine pawcprj_mpi_sum
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_mpi_allgather
!! NAME
!! pawcprj_mpi_allgather
!!
!! FUNCTION
!! Perform MPI_ALLGATHER on a pawcprj_type inside a MPI communicator.
!!
!! INPUTS
!!  cprj_loc= The cprj on the local proc being all-gathered
!!  natom=Number of atoms (size of first dimension of cprj_loc).
!!  n2dim=Size of the second dimension of cprj_loc.
!!  nlmn(natom)=Number of nlm partial waves for each atom.
!!  ncpgr = number of gradients in cprj_loc
!!  nproc=number of processors being gathered
!!  spaceComm=MPI Communicator.
!!  [rank_ordered]= optional, default=FALSE
!!                  TRUE: second dimension of gathered datastructure is rank-ordered
!!                  FALSE: second dimension of gathered datastructure is not rank-ordered
!!
!! OUTPUT
!!  cprj_gat=the gathered cprjs
!!  ierr=Error status.
!!
!! PARENTS
!!      berryphase_new,cgwf,optics_paw,optics_paw_core,suscep_stat
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine pawcprj_mpi_allgather(cprj_loc,cprj_gat,natom,n2dim,nlmn,ncpgr,nproc,spaceComm,ierr,&
&                                rank_ordered)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,n2dim,ncpgr,nproc,spaceComm
 integer,intent(out) :: ierr
 logical,optional,intent(in) :: rank_ordered
!arrays
 integer,intent(in) :: nlmn(natom)
 type(pawcprj_type),intent(in) :: cprj_loc(:,:)
 type(pawcprj_type),intent(inout) :: cprj_gat(:,:) !vz_i

!Local variables-------------------------------
!scalars
 integer :: iat,jj,t2dim,tcpgr,tg2dim,n1dim,nn
 integer :: ntotcp,ibuf,ipck,iproc
 logical :: rank_ordered_
 character(len=100) :: msg
!arrays
 real(dp),allocatable :: buffer_cpgr(:,:,:),buffer_cpgr_all(:,:,:)

! *************************************************************************

 n1dim=0
 t2dim=0
 tg2dim=0
 tcpgr=0
 ierr=0

 nn=size(nlmn,dim=1)
 n1dim=size(cprj_loc,dim=1)
 t2dim=size(cprj_loc,dim=2)
 tg2dim=size(cprj_gat,dim=2)
 tcpgr=cprj_loc(1,1)%ncpgr

 if (nn/=n1dim) then
   msg='size mismatch in natom (pawcprj_mpi_allgather)!'
   MSG_BUG(msg)
 end if
 if (t2dim/=n2dim) then
   msg='size mismatch in dim=2 (pawcprj_mpi_allgather)!'
   MSG_BUG(msg)
 end if
 if (tg2dim/=n2dim*nproc) then
   msg='size mismatch in dim=2 (pawcprj_mpi_allgather)!'
   MSG_BUG(msg)
 end if
 if (tcpgr/=ncpgr) then
   msg='size mismatch in ncpgr (pawcprj_mpi_allgather)!'
   MSG_BUG(msg)
 end if

 rank_ordered_=.false.;if(present(rank_ordered)) rank_ordered_=rank_ordered

 ntotcp=n2dim*SUM(nlmn(:))
 LIBPAW_ALLOCATE(buffer_cpgr,(2,1+ncpgr,ntotcp))
 LIBPAW_ALLOCATE(buffer_cpgr_all,(2,1+ncpgr,nproc*ntotcp))

!=== Pack cprj_loc ====
 ipck=0
 do jj=1,n2dim
   do iat=1,natom
     nn=nlmn(iat)
     buffer_cpgr(:,1,ipck+1:ipck+nn)=cprj_loc(iat,jj)%cp(:,1:nn)
     if (ncpgr/=0) buffer_cpgr(:,2:1+ncpgr,ipck+1:ipck+nn)=cprj_loc(iat,jj)%dcp(:,:,1:nn)
     ipck=ipck+nn
   end do
 end do

!=== allgather data ===
 call xmpi_allgather(buffer_cpgr,2*(ncpgr+1)*ntotcp,buffer_cpgr_all,spaceComm,ierr)

!=== unpack gathered data into cprj(natom,n2dim*nproc)
!=== second dimension is rank-ordered if rank_ordered_=true
 ipck=0
 do iproc=1,nproc
   do jj=1,n2dim
     if (rank_ordered_) then
       ibuf=jj+(iproc-1)*n2dim
     else
       ibuf=iproc+(jj-1)*nproc
     end if
     do iat=1,natom
       nn=nlmn(iat)
       cprj_gat(iat,ibuf)%cp(:,1:nn)=buffer_cpgr_all(:,1,ipck+1:ipck+nn)
       if (ncpgr/=0) cprj_gat(iat,ibuf)%dcp(:,1:ncpgr,1:nn)=&
&          buffer_cpgr_all(:,2:1+ncpgr,ipck+1:ipck+nn)
       ipck=ipck+nn
     end do
   end do
 end do

!=== Clean up ===
 LIBPAW_DEALLOCATE(buffer_cpgr)
 LIBPAW_DEALLOCATE(buffer_cpgr_all)

end subroutine pawcprj_mpi_allgather
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_bcast
!! NAME
!! pawcprj_bcast
!!
!! FUNCTION
!! Broadcast a pawcprj_type from master to all nodes inside a MPI communicator.
!!
!! INPUTS
!!  natom=Number of atoms (size of the first dimension of Cprj).
!!  n2dim=Size of the second dimension of Cprj.
!!  ncpgr=Number of gradients that have to be cast. It is a bit redundant but, it can be used to
!!   broad cast only the %cp"s without caring about the gradients. Just set it to 0 but be careful!
!!  nlmn(natom)=Number of nlm partial waves for each atom.
!!  master=ID of the sending node in spaceComm.
!!  spaceComm=MPI Communicator.
!!
!! OUTPUT
!!  ierr=Error status.
!!  Cprj(natom,n2dim)<pawcprj_type>=The datatype to be transmitted by master and received by the others nodes.
!!
!! PARENTS
!!      m_fock,posdoppler
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine pawcprj_bcast(Cprj,natom,n2dim,nlmn,ncpgr,master,spaceComm,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,n2dim,ncpgr,master,spaceComm
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: nlmn(natom)
 type(pawcprj_type),intent(inout) :: Cprj(natom,n2dim)

!Local variables-------------------------------
!scalars
 integer :: iat,jj,n1dim,nn
 integer :: ntotcp,ipck,rank,nprocs
 character(len=100) :: msg
!arrays
 real(dp),allocatable :: buffer_cp(:,:),buffer_cpgr(:,:,:)

! *************************************************************************

 ierr=0
 nprocs = xmpi_comm_size(spaceComm)
 if (nprocs==1) return

 rank = xmpi_comm_rank(spaceComm)

 nn=size(nlmn,dim=1)
 n1dim=size(Cprj,dim=1)
 if (nn/=n1dim) then
   msg='size mismatch in natom (pawcprj_bcast)!'
   MSG_BUG(msg)
 end if

 ntotcp=n2dim*SUM(nlmn(:))

 LIBPAW_ALLOCATE(buffer_cp,(2,ntotcp))
 if (ncpgr/=0)  then
   LIBPAW_ALLOCATE(buffer_cpgr,(2,ncpgr,ntotcp))
 end if

!=== Master packs Cprj ===
!Write a routine to pack/unpack?
 if (rank==master) then
   ipck=0
   do jj=1,n2dim
     do iat=1,natom
       nn=nlmn(iat)
       buffer_cp(:,ipck+1:ipck+nn)=Cprj(iat,jj)%cp(:,1:nn)
       if (ncpgr/=0) buffer_cpgr(:,:,ipck+1:ipck+nn)=Cprj(iat,jj)%dcp(:,:,1:nn)
       ipck=ipck+nn
     end do
   end do
 end if

!=== Transmit data ===
 call xmpi_bcast(buffer_cp,master,spaceComm,ierr)
 if (ncpgr/=0) then
   call xmpi_bcast(buffer_cpgr,master,spaceComm,ierr)
 end if

!=== UnPack the received buffer ===
 if (rank/=master) then
   ipck=0
   do jj=1,n2dim
     do iat=1,natom
       nn=nlmn(iat)
       Cprj(iat,jj)%cp(:,1:nn)=buffer_cp(:,ipck+1:ipck+nn)
       if (ncpgr/=0) Cprj(iat,jj)%dcp(:,:,1:nn)=buffer_cpgr(:,:,ipck+1:ipck+nn)
       ipck=ipck+nn
     end do
   end do
 end if

 LIBPAW_DEALLOCATE(buffer_cp)
 if (ncpgr/=0)  then
   LIBPAW_DEALLOCATE(buffer_cpgr)
 end if

end subroutine pawcprj_bcast
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_transpose
!! NAME
!! pawcprj_transpose
!!
!! FUNCTION
!! Transpose a cprj datastructure FOR A GIVEN (K,SPIN)
!! in order to change the parallel distribution from atom to band (or the contrary).
!! At input, cprj is distributed over bands (or atoms); at output, it is distributed over atoms (or bands)
!!
!! INPUTS
!!  cprjin(n1indim,n2indim)<pawcprj_type>=the input cprj datastructure
!!  cprj_bandpp=number of bands to be treated simultaneoulsy by a processor
!!  natom=number of atoms in cell
!!  nband=number of bands
!!  nspinor=number of spinorial components
!!  spaceComm=MPI Communicator.
!!
!! OUTPUT
!!  cprjout(n1outdim,n2outdim)<pawcprj_type>=the output cprj datastructure with another distribution
!!
!! NOTES
!!  On the dimensions:
!!   To transfer cprj from band distribution to atom distribution, dimensions should be:
!!    n1indim =natom       n2indim =nband/nproc*nspinor
!!    n1outdim=natom/nproc n2outdim=nband*nspinor
!!   To transfer cprj from atom distribution to band distribution, dimensions should be:
!!    n1indim =natom       n2indim =nband/nproc*nspinor
!!    n1outdim=natom/nproc n2outdim=nband*nspinor
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

 subroutine pawcprj_transpose(cprjin,cprjout,cprj_bandpp,natom,nband,nspinor,spaceComm)

!Arguments-------------------------------------
!scalars
 integer :: cprj_bandpp,natom,nband,nspinor,spaceComm
!arrays
 type(pawcprj_type),intent(in) ::  cprjin(:,:)
 type(pawcprj_type),intent(out) :: cprjout(:,:)

!Local variables-------------------------------
!scalars
 integer :: bpp,buf_indx
 integer :: iashft,iatom,iatom_max_sd,iatom_max_rc,iatom_1,iatom_2,iatm1_sd,iatm1_rc,iatm2_sd,iatm2_rc
 integer :: ib,iband,iband_1,iband_2,iband_shift,iblock_atom,iblock_band,ibshft
 integer :: ierr,ip,ispinor,me,nba,nbb,nbnp_sd,nbnp_rc,ncpgr,nlmn,np
 integer :: rbufsize,sbufsize,size11,size12,size21,size22,transpose_mode
 character(len=100) :: msg
!arrays
 integer,allocatable :: cprjsz_atom(:),cprjsz_block(:,:)
 integer,allocatable,target :: count_atom(:),count_band(:),displ_atom(:),displ_band(:)
 integer,pointer :: scount(:),sdispl(:),rcount(:),rdispl(:)
 real(dp),allocatable :: rbuf(:),sbuf(:)

! *************************************************************************

!MPI data
 me = xmpi_comm_rank(spaceComm)
 np = xmpi_comm_size(spaceComm)

!Nothing to do if nprocs=1
 if (np==1) then
   call pawcprj_copy(cprjin,cprjout)
   return
 end if

!Compute bloc sizes
 bpp=cprj_bandpp
 nba=natom/np;if (mod(natom,np)/=0) nba=nba+1
 nbb=nband/(np*bpp)

!Check sizes, select direction of transposition
 transpose_mode=0
 size11=size(cprjin,1);size12=size(cprjin,2)
 size21=size(cprjout,1);size22=size(cprjout,2)
 if (size11==natom.and.size12==nbb*bpp*nspinor.and.&
& size21==nba.and.size22==nband*nspinor) then
   transpose_mode=1
 else if (size11==nba.and.size12==nband*nspinor.and.&
&   size21==natom.and.size22==nbb*bpp*nspinor) then
 else
   msg='wrong cprjin/cprjout sizes (pawcprj_transpose)!'
   MSG_BUG(msg)
 end if

!Compute size of atom bloc (wr to cprj)
 LIBPAW_ALLOCATE(cprjsz_atom,(natom))
 LIBPAW_ALLOCATE(cprjsz_block,(np,nba))
 cprjsz_atom=0;cprjsz_block=0
 if (transpose_mode==1) then
   do iatom=1,natom
     cprjsz_atom(iatom)=2*cprjin(iatom,1)%nlmn*(1+cprjin(iatom,1)%ncpgr)
   end do
 else
   do iblock_atom=1,nba
     iatom=(iblock_atom-1)*np+1+me
     if (iatom<=natom) cprjsz_atom(iatom)=2*cprjin(iblock_atom,1)%nlmn*(1+cprjin(iblock_atom,1)%ncpgr)
   end do
   call xmpi_sum(cprjsz_atom,spaceComm,ierr)
 end if
 do iblock_atom=1,nba
   iashft=(iblock_atom-1)*np
   iatom_1=iashft+1;iatom_2=iashft+np
   if (iatom_1>natom) cycle
   if (iatom_2>natom) iatom_2=natom
   do iatom=iatom_1,iatom_2
     cprjsz_block(iatom-iashft,iblock_atom)=cprjsz_atom(iatom)+2  ! +2 for nlmn et ncpgr
   end do
 end do
 LIBPAW_DEALLOCATE(cprjsz_atom)

!Allocations for MPI_ALLTOALL
 LIBPAW_ALLOCATE(count_atom,(np))
 LIBPAW_ALLOCATE(displ_atom,(np))
 LIBPAW_ALLOCATE(count_band,(np))
 LIBPAW_ALLOCATE(displ_band,(np))

!Loop on blocks of bands
 do iblock_band=1,nbb !(note: np divides nband)
   ibshft=(iblock_band-1)*np*bpp
   iband_1=ibshft+1;iband_2=ibshft+np*bpp
   if (iband_1>nband.or.iband_2>nband) cycle ! for security

!  Loop on blocks of atoms
   do iblock_atom=1,nba
     iashft=(iblock_atom-1)*np
     iatom_1=iashft+1;iatom_2=iashft+np
     if (iatom_1>natom) cycle
     if (iatom_2>natom) iatom_2=natom

!    Computation of displacements and sizes of blocks when data are band-distributed
     count_band(1)=cprjsz_block(1,iblock_atom)*nspinor*bpp;displ_band(1)=0
     do ip=2,np
       count_band(ip)=cprjsz_block(ip,iblock_atom)*nspinor*bpp
       displ_band(ip)=displ_band(ip-1)+count_band(ip-1)
     end do

!    Computation of displacements and sizes of blocks when data are atom-distributed
     count_atom(1)=cprjsz_block(1+me,iblock_atom)*bpp*nspinor;displ_atom(1)=0
     do ip=2,np
       count_atom(ip)=count_atom(1)
       displ_atom(ip)=displ_atom(ip-1)+count_atom(ip-1)
     end do

!    According to transposition mode, select
!    - displacements and sizes of blocks
!    - shifts in arrays
     if (transpose_mode==1) then
       scount => count_band ; sdispl => displ_band
       rcount => count_atom ; rdispl => displ_atom
       nbnp_sd=bpp;nbnp_rc=np*bpp
       iatm1_sd=iatom_1;iatm2_sd=iatom_2
       iatm1_rc=iblock_atom;iatm2_rc=iatm1_rc
       iatom_max_sd=iatom_2;iatom_max_rc=iashft+1+me
     else
       scount => count_atom ; sdispl => displ_atom
       rcount => count_band ; rdispl => displ_band
       nbnp_sd=np*bpp;nbnp_rc=bpp
       iatm1_sd=iblock_atom;iatm2_sd=iatm1_sd
       iatm1_rc=iatom_1;iatm2_rc=iatom_2
       iatom_max_sd=iashft+1+me;iatom_max_rc=iatom_2
     end if

!    Allocation of buffers
     sbufsize=sdispl(np)+scount(np)
     rbufsize=rdispl(np)+rcount(np)
     LIBPAW_ALLOCATE(sbuf,(sbufsize))
     LIBPAW_ALLOCATE(rbuf,(rbufsize))

!    Coying of input cprj to buffer for sending
     buf_indx=0
     iband_shift=(iblock_band-1)*nbnp_sd-1
     if (iatom_max_sd<=natom) then
       do iatom=iatm1_sd,iatm2_sd
         do ib=1,nbnp_sd
           iband=(iband_shift+ib)*nspinor
           do ispinor=1,nspinor
             iband=iband+1
             nlmn=cprjin(iatom,iband)%nlmn;ncpgr=cprjin(iatom,iband)%ncpgr
             sbuf(buf_indx+1)=dble(nlmn) ;buf_indx=buf_indx+1
             sbuf(buf_indx+1)=dble(ncpgr);buf_indx=buf_indx+1
             sbuf(buf_indx+1:buf_indx+2*nlmn)=reshape(cprjin(iatom,iband)%cp(1:2,1:nlmn),(/2*nlmn/))
             buf_indx=buf_indx+2*nlmn
             if (ncpgr>0) then
               sbuf(buf_indx+1:buf_indx+2*ncpgr*nlmn)=reshape(cprjin(iatom,iband)%dcp(1:2,1:ncpgr,1:nlmn),(/2*ncpgr*nlmn/))
               buf_indx=buf_indx+2*ncpgr*nlmn
             end if
           end do
         end do
       end do
     end if
     if (buf_indx/=sbufsize) then
       msg='wrong buffer size for sending (pawcprj_transpose)!'
       MSG_BUG(msg)
     end if

!    Main call to MPI_ALLTOALL
     call xmpi_alltoallv(sbuf,scount,sdispl,rbuf,rcount,rdispl,spaceComm,ierr)

!    Retrieving of output cprj for received buffer
     buf_indx=0
     iband_shift=(iblock_band-1)*nbnp_rc-1
     if (iatom_max_rc<=natom) then
       do iatom=iatm1_rc,iatm2_rc
         do ib=1,nbnp_rc
           iband=(iband_shift+ib)*nspinor
           do ispinor=1,nspinor
             iband=iband+1
             nlmn =int(rbuf(buf_indx+1));buf_indx=buf_indx+1
             ncpgr=int(rbuf(buf_indx+1));buf_indx=buf_indx+1
             cprjout(iatom,iband)%nlmn=nlmn;cprjout(iatom,iband)%ncpgr=ncpgr
             cprjout(iatom,iband)%cp(1:2,1:nlmn)=reshape(rbuf(buf_indx+1:buf_indx+2*nlmn),(/2,nlmn/))
             buf_indx=buf_indx+2*nlmn
             if (ncpgr>0) then
               cprjout(iatom,iband)%dcp(1:2,1:ncpgr,1:nlmn)=reshape(rbuf(buf_indx+1:buf_indx+2*nlmn*ncpgr),(/2,ncpgr,nlmn/))
               buf_indx=buf_indx+2*nlmn*ncpgr
             end if
           end do
         end do
       end do
     else
       cprjout(iatom,iband)%nlmn=0;cprjout(iatom,iband)%ncpgr=0
     end if
     if (buf_indx/=rbufsize) then
       msg='wrong buffer size for receiving (pawcprj_transpose)!'
       MSG_BUG(msg)
     end if

!    Deallocation of buffers
     LIBPAW_DEALLOCATE(sbuf)
     LIBPAW_DEALLOCATE(rbuf)

!    End of loops
   end do ! do iblock_atom
 end do ! do iblock_atom

!Free memory
 LIBPAW_DEALLOCATE(count_atom)
 LIBPAW_DEALLOCATE(displ_atom)
 LIBPAW_DEALLOCATE(count_band)
 LIBPAW_DEALLOCATE(displ_band)
 LIBPAW_DEALLOCATE(cprjsz_block)
 nullify(scount,rcount,sdispl,rdispl)

 end subroutine pawcprj_transpose
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_gather_spin
!! NAME
!! pawcprj_gather_spin
!!
!! FUNCTION
!!
!! INPUTS
!!  cprj(:,:)=the input cprj datastructure
!!  n2size=number of cprj datastructures to be gathered (second dim)
!!  nspinor : number of spinorial component (on current proc)
!!  nspinortot : total number of spinorial component
!!
!! OUTPUT
!!  cprj_gat(:,:) = the cprj containing all nspinor componants
!!
!! NOTES
!! The cprj has been built like the following:
!!   loop on nsppol
!!   loop on k point
!!   loop over band or block of band
!! These quantities were build only if treated by the current proc
!! the inner quantities being nspinor
!!
!! PARENTS
!!      energy,pawmkrhoij
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE
 subroutine pawcprj_gather_spin(cprj,cprj_gat,natom,n2size,nspinor,nspinortot,&
&                            spaceComm_spin,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nspinor,nspinortot,n2size
 integer,intent(in) :: spaceComm_spin
 integer,intent(out) :: ierr
!arrays
 type(pawcprj_type),intent(in) :: cprj(:,:)
 type(pawcprj_type),intent(inout) :: cprj_gat(:,:)

!Local variables-------------------------------
!scalars
 integer :: i1,iatom,ibsp,icpgr,ilmn,isp,ispinor,jj,lmndim,n2dim,n2dim_gat,ncpgr
 character(len=100) :: msg
!arrays
 integer :: nlmn(natom)
 real(dp),allocatable :: buffer1(:),buffer2(:)

! *************************************************************************

 n2dim    =size(cprj,dim=2)
 n2dim_gat=size(cprj_gat,dim=2)
 if (n2dim_gat/=(nspinortot/nspinor)*n2dim) then
   msg='wrong dims (pawcprj_gather_spin)!'
   MSG_BUG(msg)
 end if

 do iatom=1,natom
   nlmn(iatom)=size(cprj(iatom,1)%cp(1,:))
 end do
 ncpgr=cprj(1,1)%ncpgr
 lmndim=2*n2size*sum(nlmn(1:natom))*(1+ncpgr)
 LIBPAW_ALLOCATE(buffer1,(lmndim))
 LIBPAW_ALLOCATE(buffer2,(lmndim*nspinortot))

 isp=0;ibsp=0
 jj=1
 do i1=1,n2size
   isp=isp+1
   do iatom=1,natom
     do ilmn=1,nlmn(iatom)
       buffer1(jj:jj+1)=cprj(iatom,isp)%cp(1:2,ilmn)
       jj=jj+2
     end do
     if (ncpgr>0) then
       do ilmn=1,nlmn(iatom)
         do icpgr=1,ncpgr
           buffer1(jj:jj+1)=cprj(iatom,isp)%dcp(1:2,icpgr,ilmn)
           jj=jj+2
         end do
       end do
     end if
   end do
 end do

 call xmpi_allgather(buffer1,lmndim,buffer2,spaceComm_spin,ierr)

 jj=1
 do ispinor=1,nspinortot
   do i1 =1,n2size
     ibsp=(i1-1)*nspinortot + ispinor
     do iatom=1,natom
       do ilmn=1,nlmn(iatom)
         cprj_gat(iatom,ibsp)%cp(1:2,ilmn)=buffer2(jj:jj+1)
         jj=jj+2
       end do
       if (ncpgr>0) then
         do ilmn=1,nlmn(iatom)
           do icpgr=1,ncpgr
             cprj_gat(iatom,ibsp)%dcp(1:2,icpgr,ilmn)=buffer2(jj:jj+1)
             jj=jj+2
           end do
         end do
       end if
     end do
   end do
 end do

 LIBPAW_DEALLOCATE(buffer1)
 LIBPAW_DEALLOCATE(buffer2)

 end subroutine pawcprj_gather_spin
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_getdim
!! NAME
!! pawcprj_getdim
!!
!! FUNCTION
!!  Helper function returning the number of lmn components in the <p_{lmn}^i|\psi> for the i-th atom.
!!  Used to initialize the dimensioning array that is passed to the pawcprj_alloc routines when the
!!  pawcprj_type structure is allocated and initialized.
!!
!! INPUTS
!! natom=number of atoms in the unit cell
!! nattyp(ntypat)=number of atoms of each type
!! ntypat=number of atom types
!! typat(natom-= type of each atom
!! Pawtab(ntypat)<pawtab_type>=PAW tabulated starting data.
!! sort_mode(len=*)=String defining the sorting of the atoms in the Cprj arrays.
!!   Two modes are possible:
!!   -- "O[rdered]", if atoms are sorted by atom type.
!!   -- "R[andom]", if atoms are sorted randomly i.e. according the values of typat specified in the input file.
!!
!! OUTPUT
!!  dimcprj(natom)=Number of nlm elements in the <p_{lmn}^i|\psi> matrix elements for i=1,...,natom.
!!
!! PARENTS
!!      afterscfloop,berryphase_new,chern_number,dfpt_looppert,dfpt_scfcv
!!      extrapwf,forstr,getghc,initberry,m_fock,m_hamiltonian,mlwfovlp_qp
!!      outkss,scfcv,smatrix_pawinit,wf_mixing
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine pawcprj_getdim(dimcprj,natom,nattyp,ntypat,typat,Pawtab,sort_mode)

!Arguments ------------------------------------
 integer,intent(in) :: natom,ntypat
 character(len=*),intent(in) :: sort_mode
!arrays
 integer,intent(in) :: nattyp(:),typat(natom)
 integer,intent(inout) :: dimcprj(natom)
 type(Pawtab_type),intent(in) :: Pawtab(ntypat)

!Local variables-------------------------------
 integer :: iatom,itypat
 character(len=500) :: msg

! *************************************************************************

 SELECT CASE (sort_mode(1:1))

 CASE ("o","O") ! Ordered by atom-type

  iatom=0
  do itypat=1,ntypat
   dimcprj(iatom+1:iatom+nattyp(itypat))=Pawtab(itypat)%lmn_size
   iatom=iatom+nattyp(itypat)
  end do

 CASE ("r","R") ! Randomly ordered (typat from input file)

  do iatom=1,natom
   itypat=typat(iatom)
   dimcprj(iatom)=Pawtab(itypat)%lmn_size
  end do

 CASE DEFAULT
  msg='Wrong value for sort_mode: '//TRIM(sort_mode)
  MSG_ERROR(msg)
 END SELECT

end subroutine pawcprj_getdim
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/paw_overlap
!! NAME
!! paw_overlap
!!
!! FUNCTION
!!  Helper function returning the onsite contribution to the overlap between two states.
!!
!! INPUTS
!!   spinor_comm= (optional) communicator over spinorial components
!!   typat(:)=The type of each atom.
!!   Pawtab(ntypat)<type(pawtab_type)>=paw tabulated starting data.
!!   cprj1,cprj2<pawcprj_type>
!!     Projected wave functions <Proj_i|Cnk> with all NL projectors for the left and the right wavefunction,respectively.
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!  xmpi_sum
!!
!! SOURCE

function paw_overlap(cprj1,cprj2,typat,pawtab,spinor_comm) result(onsite)

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: spinor_comm
!arrays
 integer,intent(in) :: typat(:)
 real(dp) :: onsite(2)
 type(pawcprj_type),intent(in) :: cprj1(:,:),cprj2(:,:)
 type(pawtab_type),intent(in) :: pawtab(:)

!Local variables-------------------------------
!scalars
 integer :: iatom,ilmn,itypat,j0lmn,jlmn,klmn,natom,nspinor,isp
 real(dp) :: sij
 character(len=500) :: msg
!arrays

! *************************************************************************

 natom=SIZE(typat)

 if (SIZE(cprj1,DIM=1)/=SIZE(cprj2,DIM=1) .or. SIZE(cprj1,DIM=1)/=natom) then
   write(msg,'(a,3i4)')' Wrong size in typat, cprj1, cprj2 : ',natom,SIZE(cprj1),SIZE(cprj2)
   MSG_ERROR(msg)
 end if

 nspinor = SIZE(cprj1,DIM=2)

 onsite=zero
 do iatom=1,natom
   itypat=typat(iatom)
   do jlmn=1,pawtab(itypat)%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       sij=pawtab(itypat)%sij(klmn); if (jlmn==ilmn) sij=sij*half
       if (ABS(sij)>tol16) then
         do isp=1,nspinor

           onsite(1)=onsite(1) + sij*(                                 &
&            cprj1(iatom,isp)%cp(1,ilmn) * cprj2(iatom,isp)%cp(1,jlmn) &
&           +cprj1(iatom,isp)%cp(2,ilmn) * cprj2(iatom,isp)%cp(2,jlmn) &
&           +cprj1(iatom,isp)%cp(1,jlmn) * cprj2(iatom,isp)%cp(1,ilmn) &
&           +cprj1(iatom,isp)%cp(2,jlmn) * cprj2(iatom,isp)%cp(2,ilmn) &
&           )

           onsite(2)=onsite(2) + sij*(                                 &
&           cprj1(iatom,isp)%cp(1,ilmn) * cprj2(iatom,isp)%cp(2,jlmn)  &
&           -cprj1(iatom,isp)%cp(2,ilmn) * cprj2(iatom,isp)%cp(1,jlmn) &
&           +cprj1(iatom,isp)%cp(1,jlmn) * cprj2(iatom,isp)%cp(2,ilmn) &
&           -cprj1(iatom,isp)%cp(2,jlmn) * cprj2(iatom,isp)%cp(1,ilmn) &
&           )
         end do
       end if
     end do
   end do
 end do

 if (present(spinor_comm)) then
   call xmpi_sum(onsite,spinor_comm,isp)
 end if

end function paw_overlap
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_pack
!! NAME
!! pawcprj_pack
!!
!! FUNCTION
!! Pack structured data into a simple buffer
!!
!! INPUTS
!!  nlmn(natom)=Number of nlm partial waves for each atom.
!!  ncpgr = number of gradients in cprj_out
!!  cprj= The datatype to be packed.
!!
!! OUTPUT
!!  buffer = the data packed, dim : (2, n2dim*sum(nlmn))
!!  [buffer_gr] = if present the gradient data packed, dim : (2, ncpgr, n2dim*sum(nlmn))
!!
!! PARENTS
!!  pawmkrhoij
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawcprj_pack(nlmn,cprj,buffer,buffer_gr)

!Arguments ------------------------------------
!scalars
!arrays
 integer,intent(in) :: nlmn(:)
 type(pawcprj_type),intent(in) :: cprj(:,:)
 real(dp),intent(out) :: buffer(:,:)
 real(dp),intent(out),optional :: buffer_gr(:,:,:)

!Local variables-------------------------------
!scalars
 integer :: natom,n2buffer,ncpgr,n2dim
 integer :: iat,jj,n1dim,nn
 integer :: ipck
 character(len=100) :: msg
!arrays

! *************************************************************************

 natom=size(nlmn,dim=1)
 n2buffer=size(buffer,dim=2)
 n1dim=size(cprj,dim=1)
 n2dim=size(cprj,dim=2)

 if (natom/=n1dim) then
   msg='size mismatch in natom (pawcprj_pack)!'
   MSG_BUG(msg)
 end if
 if (n2dim*SUM(nlmn)/=n2buffer) then
   msg='size mismatch in dim=2 (pawcprj_pack)!'
   MSG_BUG(msg)
 end if
 ncpgr=0
 if (present(buffer_gr)) then
   ncpgr=size(buffer_gr,dim=2)
 end if

!=== Pack cprj ====
 ipck=0
 do jj=1,n2dim
   do iat=1,natom
     nn=nlmn(iat)
     buffer(:,ipck+1:ipck+nn)=cprj(iat,jj)%cp(:,1:nn)
     if (ncpgr/=0) then
       buffer_gr(:,:,ipck+1:ipck+nn)=cprj(iat,jj)%dcp(:,:,1:nn)
     end if
     ipck=ipck+nn
   end do
 end do

end subroutine pawcprj_pack
!!***

!----------------------------------------------------------------------

!!****f* m_pawcprj/pawcprj_unpack
!! NAME
!! pawcprj_unpack
!!
!! FUNCTION
!! Unpack structured data from a simple buffer
!!
!! INPUTS
!!  nlmn(natom)=Number of nlm partial waves for each atom.
!!  ncpgr = number of gradients in cprj_in
!!  buffer = the data to be unpacked, dim : (2, n2dim*sum(nlmn))
!!  [buffer_gr] = if present the gradient data to be unpacked, dim : (2, ncpgr, n2dim*sum(nlmn))
!!
!! OUTPUT
!!  cprj=The datatype unpacked
!!
!! PARENTS
!!  pawmkrhoij
!!
!! CHILDREN
!!
!! SOURCE

subroutine pawcprj_unpack(nlmn,cprj,buffer,buffer_gr)

!Arguments ------------------------------------
!scalars
!arrays
 integer,intent(in) :: nlmn(:)
 real(dp),intent(in) :: buffer(:,:)
 real(dp),intent(in),optional :: buffer_gr(:,:,:)
 type(pawcprj_type),intent(inout) :: cprj(:,:)

!Local variables-------------------------------
!scalars
 integer :: natom,n2buffer,ncpgr,n2dim
 integer :: iat,jj,n1dim,nn
 integer :: ipck
 character(len=100) :: msg
!arrays

! *************************************************************************

 natom=size(nlmn,dim=1)
 n2buffer=size(buffer,dim=2)
 n1dim=size(cprj,dim=1)
 n2dim=size(cprj,dim=2)

 if (natom/=n1dim) then
   msg='size mismatch in natom (pawcprj_unpack)!'
   MSG_BUG(msg)
 end if
 if (n2dim*SUM(nlmn)/=n2buffer) then
   msg='size mismatch in dim=2 (pawcprj_unpack)!'
   MSG_BUG(msg)
 end if
 ncpgr=0
 if (present(buffer_gr)) then
   ncpgr=size(buffer_gr,dim=2)
 end if

!=== Unpack buffers into cprj ===
 ipck=0
 do jj=1,n2dim
   do iat=1,natom
     nn=nlmn(iat)
     cprj(iat,jj)%cp(:,1:nn)=buffer(:,ipck+1:ipck+nn)
     if (ncpgr/=0) then
       cprj(iat,jj)%dcp(:,:,1:nn)=buffer_gr(:,:,ipck+1:ipck+nn)
     end if
     ipck=ipck+nn
   end do
 end do

end subroutine pawcprj_unpack

end module m_pawcprj
!!***
