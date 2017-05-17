!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_abi2big
!! NAME
!!  m_abi2big
!!
!! FUNCTION
!!  Module to copy objects from ABINIT to BigDFT and viceversa.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2017 ABINIT group (TR,DC,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_abi2big
    
 use defs_basis
 use m_errors
 use defs_wvltypes
use m_xmpi
!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'm_abi2big'
!End of the abilint section

 implicit none
 
 private
 
 public :: wvl_vtrial_abi2big 
 !to copy vtrial to wvl_den%rhov and viceversa.

 public :: wvl_rho_abi2big
 !to copy a density from ABINIT to BigDFT or viceversa

 public :: wvl_vhartr_abi2big
 ! to copy V_hartree from ABINIT to BigDFT and viceversa

 public :: wvl_vxc_abi2big
 ! to copy Vxc from ABINIT to BigDFT and viceversa

 public :: wvl_occ_abi2big
 ! to copy occupations from/to ABINIT to/from BigDFT
 
 public :: wvl_eigen_abi2big
 ! to copy eigenvalues from/to ABINIT to/from BigDFT

 public :: wvl_occopt_abi2big
 ! maps occupation method in ABINIT and BigDFT

 public :: wvl_rhov_abi2big
 !generic routine to copy a density or potential from/to 
 !ABINIT to/from BigDFT.
 interface wvl_rhov_abi2big
   module procedure wvl_rhov_abi2big_2D_4D
   module procedure wvl_rhov_abi2big_1D_4D
   module procedure wvl_rhov_abi2big_2D_2D
   module procedure wvl_rhov_abi2big_1D_2D
   module procedure wvl_rhov_abi2big_2D_1D
   module procedure wvl_rhov_abi2big_1D_1D
 end interface wvl_rhov_abi2big

 logical,parameter :: hmem=.false. !high memory
!!  Set hmem=.false. if memory is limited. It will copy element by element.
!!  If  hmem=.true. all elements are copied at once.

contains
!!***

!!****f* m_abi2big/wvl_vtrial_abi2big
!! NAME
!!  wvl_vtrial_abi2big
!!
!! FUNCTION
!!  Copies vtrial in ABINIT to BigDFT objects and viceversa.
!!
!! INPUTS
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!  vtrial(nfft,nspden)= trial potential (ABINIT array)
!!  wvl_den= density-potential BigDFT object
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! vtrial is copied to wvl_den, or viceversa, depending on "opt" (see above).
!! It verifies that (or sets) wvl_den%rhov_is = KS_POTENTIAL.
!!
!! NOTES
!! It uses the generic routine wvl_rhov_abi2big.
!!
!! PARENTS
!!      afterscfloop,newvtr,rhotov,setvtr,wvl_psitohpsi
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_vtrial_abi2big(opt,vtrial,wvl_den)
    
#if defined HAVE_BIGDFT
  use BigDFT_API, only : KS_POTENTIAL
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_vtrial_abi2big'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: opt
 real(dp),intent(inout) :: vtrial(:,:)
 type(wvl_denspot_type),intent(inout) :: wvl_den

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer :: shiftV
 character(len=100) :: message
#endif
 
! *************************************************************************
 
 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT
! write(message,'(2a)') ch10,' wvl_vtrial_abi2big: but why are you copying me :..o('
! call wrtout(std_out,message,'COLL')

 shiftV=wvl_den%denspot%dpbox%ndims(1)*wvl_den%denspot%dpbox%ndims(2) &
&      *wvl_den%denspot%dpbox%i3xcsh

 if(opt==1) then !ABINIT -> BIGDFT

   call wvl_rhov_abi2big(opt,vtrial,wvl_den%denspot%rhov,shift=shiftV)
   wvl_den%denspot%rhov_is = KS_POTENTIAL

 elseif(opt==2) then !BigDFT -> ABINIT

   if(wvl_den%denspot%rhov_is .ne. KS_POTENTIAL) then
     message='wvl_vtrial_abi2big: rhov should contain the KS_POTENTIAL'
     MSG_BUG(message)
   end if

   call wvl_rhov_abi2big(opt,vtrial,wvl_den%denspot%rhov,shift=shiftV)

 else
   message='wvl_vtrial_abi2big: wrong option'
   MSG_BUG(message)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*)  opt,vtrial(1,1),wvl_den%symObj
#endif

 DBG_EXIT("COLL")

end subroutine wvl_vtrial_abi2big
!!***

!!****f* m_abi2big/wvl_rho_abi2big
!! NAME
!!  wvl_rho_abi2big
!!
!! FUNCTION
!!  Copies the density from ABINIT to BigDFT, or viceversa.
!!
!! INPUTS
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!  rhor(nfft,nspden)= trial potential (ABINIT array)
!!  wvl_den= density-potential BigDFT object
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! Density copied from ABINIT to BigDFT or viceversa.
!! It verifies that (or sets) wvl_den%rhov_is= ELECTRONIC_DENSITY.
!!
!! NOTES
!! It uses the generic routine wvl_rhov_abi2big.
!!
!! PARENTS
!!      afterscfloop,mklocl,newrho,vtorho,wvl_mkrho
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_rho_abi2big(opt,rhor,wvl_den)
    
#if defined HAVE_BIGDFT
  use BigDFT_API, only : ELECTRONIC_DENSITY
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_rho_abi2big'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: opt
 real(dp) , intent(inout)  :: rhor(:,:)
 type(wvl_denspot_type), intent(inout) :: wvl_den

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 character(len=100) :: message
#endif
 
! *************************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT
! write(message,'(2a)') ch10,'wvl_rho_abi2big: but why are you copying me :..o('
! call wrtout(std_out,message,'COLL')

 if(opt==1) then !ABINIT -> BIGDFT

   call wvl_rhov_abi2big(opt,rhor,wvl_den%denspot%rhov)
   wvl_den%denspot%rhov_is = ELECTRONIC_DENSITY

 elseif(opt==2) then !BigDFT -> ABINIT

   if(wvl_den%denspot%rhov_is .ne. ELECTRONIC_DENSITY) then
     message='wvl_rho_abi2big: rhov should contain the ELECTRONIC_DENSITY'
     MSG_BUG(message)
   end if
   call wvl_rhov_abi2big(opt,rhor,wvl_den%denspot%rhov)

 else
   message='wvl_rho_abi2big: wrong option'
   MSG_BUG(message)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*)  opt,rhor(1,1),wvl_den%symObj
#endif

 DBG_EXIT("COLL")

end subroutine wvl_rho_abi2big
!!***

!!****f* m_abi2big/wvl_vhartr_abi2big
!! NAME
!!  wvl_vhartr_abi2big
!!
!! FUNCTION
!!  Copies vhartree in ABINIT to BigDFT objects and viceversa.
!!
!! INPUTS
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!  vhartr(nfft)= Hartree potential (ABINIT array)
!!  wvl_den= density-potential BigDFT object
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! vhartr is copied to wvl_den, or viceversa, depending on "opt" (see above).
!! It verifies that (or sets) wvl_den%rhov_is = HARTREE_POTENTIAL
!!
!! NOTES
!! It uses the generic routine wvl_rhov_abi2big.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_vhartr_abi2big(opt,vhartr,wvl_den)
    
#if defined HAVE_BIGDFT
  use BigDFT_API, only : HARTREE_POTENTIAL
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_vhartr_abi2big'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: opt
 real(dp) , intent(inout)  :: vhartr(:)
 type(wvl_denspot_type), intent(inout) :: wvl_den

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer :: shiftV
 character(len=100) :: message
#endif
 
! *************************************************************************
 
 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT
! write(message,'(2a)') ch10, 'wvl_vhartr_abi2big: but why are you copying me :..o('
! call wrtout(std_out,message,'COLL')

 shiftV=wvl_den%denspot%dpbox%ndims(1)*wvl_den%denspot%dpbox%ndims(2) &
&      *wvl_den%denspot%dpbox%i3xcsh

 if(opt==1) then !ABINIT -> BIGDFT

   call wvl_rhov_abi2big(opt,vhartr,wvl_den%denspot%rhov,shift=shiftV)
   wvl_den%denspot%rhov_is = HARTREE_POTENTIAL

 elseif(opt==2) then !BigDFT -> ABINIT

   if(wvl_den%denspot%rhov_is .ne. HARTREE_POTENTIAL) then
     message='wvl_vhartr_abi2big: rhov should contain the HARTREE_POTENTIAL'
     MSG_BUG(message)
   end if
   call wvl_rhov_abi2big(opt,vhartr,wvl_den%denspot%rhov,shift=shiftV)

 else
   message='wvl_vhartr_abi2big: wrong option'
   MSG_BUG(message)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*)  opt,vhartr(1),wvl_den%symObj
#endif

 DBG_EXIT("COLL")

end subroutine wvl_vhartr_abi2big
!!***

!!****f* m_abi2big/wvl_vxc_abi2big
!! NAME
!!  wvl_vxc_abi2big
!!
!! FUNCTION
!!  It copies the Vxc potential from ABINIT to BigDFT or viceversa.
!!
!! INPUTS
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!  vxc(nfft,nspden)= trial potential (ABINIT array)
!!  wvl_den= density-potential BigDFT object
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! vxc is copied to wvl_den, or viceversa, depending on "opt" (see above).
!!
!! NOTES
!! It uses the generic routine wvl_rhov_abi2big.
!!
!! PARENTS
!!      wvl_psitohpsi
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_vxc_abi2big(opt,vxc,wvl_den)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_vxc_abi2big'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: opt
 real(dp),intent(inout) :: vxc(:,:)
 type(wvl_denspot_type), intent(inout) :: wvl_den

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer :: shiftV
#endif
 
! *************************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT
! write(message,'(2a)') ch10, 'wvl_vxc_abi2big: but why are you copying me :..o('
! call wrtout(std_out,message,'COLL')

 shiftV=wvl_den%denspot%dpbox%ndims(1)*wvl_den%denspot%dpbox%ndims(2) &
&      *wvl_den%denspot%dpbox%i3xcsh

 call wvl_rhov_abi2big(opt,vxc,wvl_den%denspot%v_xc,shift=shiftV)

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) opt,vxc(1,1),wvl_den%symObj
#endif

 DBG_EXIT("COLL")

end subroutine wvl_vxc_abi2big
!!***

!!****f* m_abi2big/wvl_occ_abi2big
!! NAME
!!  wvl_occ_abi2big
!!
!! FUNCTION
!!  Copies occupations in ABINIT to BigDFT objects and viceversa.
!!
!! INPUTS
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!  nsppol= number of spin polarization
!!
!! OUTPUT
!!  
!! SIDE EFFECTS
!! occ is copied to wfs%ks%orbs%occup, or viceversa, depending on "opt" (see above).
!!
!! PARENTS
!!      afterscfloop,gstate,vtorho,wvl_wfsinp_disk,wvl_wfsinp_scratch
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_occ_abi2big(mband,nkpt,nsppol,occ,opt,wvl_wfs)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_occ_abi2big'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: mband,nkpt,nsppol,opt
 real(dp) , intent(inout)  :: occ(mband*nkpt*nsppol)
 type(wvl_wf_type), intent(inout) :: wvl_wfs

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer :: norb,norbd,norbu,ii
 character(len=100) :: message
#endif
 
! *************************************************************************
 
 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT
!PENDING: I am not sure this will work for nsppol==2
!check also the parallel case.

 norbu=wvl_wfs%ks%orbs%norbu
 norbd=wvl_wfs%ks%orbs%norbd
 norb =wvl_wfs%ks%orbs%norb
 if(opt==1) then !ABINIT -> BIGDFT
   if (nsppol == 1) then
    do ii=1,norb
      wvl_wfs%ks%orbs%occup(ii)=occ(ii)
    end do
   else
     wvl_wfs%ks%orbs%occup(1:norbu)=occ(1:norbu)
     wvl_wfs%ks%orbs%occup(norbu + 1:norb)= &
&       occ(mband + 1:mband + norbd)
   end if
 elseif(opt==2) then !BigDFT -> ABINIT
   if (nsppol == 1) then
     do ii=1,norb
       occ=wvl_wfs%ks%orbs%occup
     end do
   else
     occ(1:norbu) = wvl_wfs%ks%orbs%occup(1:norbu)
     occ(mband + 1:mband + norbd) = &
&     wvl_wfs%ks%orbs%occup(norbu + 1:norb)
   end if
 else
   message='wvl_occ_abi2big: wrong option'
   MSG_BUG(message)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) mband,nkpt,nsppol,opt,occ(1),wvl_wfs%ks
#endif

 DBG_EXIT("COLL")

end subroutine wvl_occ_abi2big
!!***

!!****f* m_abi2big/wvl_eigen_abi2big
!! NAME
!!  wvl_eigen_abi2big
!!
!! FUNCTION
!!  Copies eigenvalues in ABINIT to BigDFT objects and viceversa.
!!
!! INPUTS
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!  nsppol= number of spin polarization
!!
!! OUTPUT 
!!
!! SIDE EFFECTS
!! occ is copied to wfs%ks%orbs%occup, or viceversa, depending on "opt" (see above).
!!
!! PARENTS
!!      afterscfloop,vtorho
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_eigen_abi2big(mband,nkpt,nsppol,eigen,opt,wvl_wfs)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_eigen_abi2big'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: mband,nkpt,nsppol,opt
 real(dp) , intent(inout)  :: eigen(mband*nkpt*nsppol)
 type(wvl_wf_type), intent(inout) :: wvl_wfs

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer :: ii,norb,norbd,norbu
 character(len=100) :: message
#endif
 
! *************************************************************************
 
 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT
!PENDING: I am not sure this will work for nsppol==2
!check also the parallel case.

 norbu=wvl_wfs%ks%orbs%norbu
 norbd=wvl_wfs%ks%orbs%norbd
 norb =wvl_wfs%ks%orbs%norb
 if(opt==1) then !ABINIT -> BIGDFT
   if (nsppol == 1) then
    wvl_wfs%ks%orbs%eval=eigen
   else
     wvl_wfs%ks%orbs%eval(1:norbu)=eigen(1:norbu)
     wvl_wfs%ks%orbs%eval(norbu + 1:norb)= &
&       eigen(mband + 1:mband + norbd)
   end if
 elseif(opt==2) then !BigDFT -> ABINIT
   if (nsppol == 1) then
     do ii=1,norb
       eigen(ii)=wvl_wfs%ks%orbs%eval(ii)
     end do
   else
     eigen(1:norbu) = wvl_wfs%ks%orbs%eval(1:norbu)
     eigen(mband + 1:mband + norbd) = &
&     wvl_wfs%ks%orbs%eval(norbu + 1:norb)
   end if
 else
   message='wvl_eigen_abi2big: wrong option'
   MSG_BUG(message)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) mband,nkpt,nsppol,opt,eigen(1),wvl_wfs%ks
#endif

 DBG_EXIT("COLL")

end subroutine wvl_eigen_abi2big
!!***

!!****f* m_abi2big/wvl_occopt_abi2big
!! NAME
!!  wvl_occopt_abi2big
!!
!! FUNCTION
!!  Copies occopt in ABINIT to BigDFT objects and viceversa.
!!
!! INPUTS
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Several smearing schemes do not exists in both codes such 
!! as the SMEARING_DIST_ERF in BigDFT.
!!
!! PARENTS
!!      vtorho,wvl_wfsinp_scratch
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_occopt_abi2big(occopt_abi,occopt_big,opt)
    
#if defined HAVE_BIGDFT
 use BigDFT_API, only : &
&  SMEARING_DIST_FERMI, SMEARING_DIST_COLD1, SMEARING_DIST_COLD2,&
&  SMEARING_DIST_METPX
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_occopt_abi2big'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer , intent(inout)  :: occopt_abi,occopt_big
 integer , intent(in)     :: opt

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 character(len=500) :: message
#endif
 
! *************************************************************************
 
 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT

 if(opt==1) then !ABINIT -> BIGDFT
   if(occopt_abi==3) then
     occopt_big=SMEARING_DIST_FERMI
   elseif(occopt_abi==4) then
     occopt_big=SMEARING_DIST_COLD1 
   elseif(occopt_abi==5) then
     occopt_big=SMEARING_DIST_COLD2
   elseif(occopt_abi==6) then
     occopt_big=SMEARING_DIST_METPX
   else
     write(message,'(4a)') ch10,&
&     ' wvl_occopt_abi2big: occopt does not have a corresponding option in BigDFT.',ch10,&
&     ' Action: change the value of occopt to a number between 3 and 6'
     MSG_ERROR(message)
   end if
 elseif(opt==2) then !BigDFT -> ABINIT
   if(occopt_big==SMEARING_DIST_FERMI) then
     occopt_abi=3
   elseif(occopt_big==SMEARING_DIST_COLD1) then
     occopt_abi=4
   elseif(occopt_big==SMEARING_DIST_COLD2) then
     occopt_abi=5
   elseif(occopt_big==SMEARING_DIST_METPX) then
     occopt_abi=6
   else
!    One should never get here.
     write(message,'(4a)') ch10,&
&     ' wvl_occopt_abi2big: occopt in BigDFT does not have a corresponding option in ABINIT.',ch10,&
&     ' Action: contact the ABINIT group'
     MSG_ERROR(message)
   end if
 else
   message='wvl_occopt_abi2big: wrong option'
   MSG_BUG(message)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) occopt_abi,occopt_big,opt
#endif

 DBG_EXIT("COLL")

end subroutine wvl_occopt_abi2big
!!***

!!****f* m_abi2big/wvl_rhov_abi2big_2D_4D
!! NAME
!!  wvl_rhov_abi2big_2D_4D
!!
!! FUNCTION
!!  Copies a density/potential from ABINIT to BigDFT or viceversa.
!!  Target : ABINIT 2D arrays (with spin), BigDFT 4D arrays (with spin)
!!
!! INPUTS
!!  opt= 1: copy from ABINIT to BigDFT, 2: copy from BigDFT to ABINIT
!!  rhov_abi(:,:)     = density/potential array in ABINIT
!!  rhov_big(:,:,:,:) = density/potential array in BigDFT
!!  [shift] = shift to be applied in rhov_abi array (parallelism)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  At output rhov_abi is copied to rhov_abinit, or viceversa
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_rhov_abi2big_2D_4D(opt,rhov_abi,rhov_big,shift)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_rhov_abi2big_2D_4D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: opt
 integer,intent(in),optional :: shift
 real(dp) :: rhov_abi(:,:),rhov_big(:,:,:,:)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer :: nfft_abi,nfft_big,nspden,shift_
 character(len=100) :: message
#endif
 
! *************************************************************************

#if defined HAVE_BIGDFT
 nspden=size(rhov_abi,2)
 if (size(rhov_big,4)/=nspden) then
   message='wvl_rhov_abi2big: ABINIT and BigDFT objects do not have the same nspden'
   MSG_BUG(message)
 end if
 nfft_abi=size(rhov_abi)/nspden ; nfft_big=size(rhov_big)/nspden
 shift_=0;if (present(shift)) shift_=shift
 call wvl_rhov_abi2big_gen(nfft_abi,nfft_big,nspden,opt,rhov_abi,rhov_big,shift_)
#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*)  opt,rhov_abi(1,1),rhov_big(1,1,1,1)
#endif

end subroutine wvl_rhov_abi2big_2D_4D
!!***

!!****f* m_abi2big/wvl_rhov_abi2big_1D_4D
!! NAME
!!  wvl_rhov_abi2big_1D_4D
!!
!! FUNCTION
!!  Copies a density/potential from ABINIT to BigDFT or viceversa.
!!  Target : ABINIT 1D arrays (without spin), BigDFT 4D arrays (with spin)
!!
!! INPUTS
!!  opt= 1: copy from ABINIT to BigDFT, 2: copy from BigDFT to ABINIT
!!  rhov_abi(:)       = density/potential array in ABINIT
!!  rhov_big(:,:,:,:) = density/potential array in BigDFT
!!  [shift] = shift to be applied in rhov_abi array (parallelism)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  At output rhov_abi is copied to rhov_abinit, or viceversa
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_rhov_abi2big_1D_4D(opt,rhov_abi,rhov_big,shift)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_rhov_abi2big_1D_4D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: opt
 integer,intent(in),optional :: shift
 real(dp) :: rhov_abi(:),rhov_big(:,:,:,:)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer :: nfft_abi,nfft_big,shift_
 character(len=100) :: message
#endif
 
! *************************************************************************

#if defined HAVE_BIGDFT
 if (size(rhov_big,4)/=1) then
   message='wvl_rhov_abi2big: ABINIT and BigDFT objects do not have the same nspden'
   MSG_BUG(message)
 end if
 nfft_abi=size(rhov_abi) ; nfft_big=size(rhov_big)
 shift_=0;if (present(shift)) shift_=shift
 call wvl_rhov_abi2big_gen(nfft_abi,nfft_big,1,opt,rhov_abi,rhov_big,shift_)
#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*)  opt,rhov_abi(1),rhov_big(1,1,1,1)
#endif

end subroutine wvl_rhov_abi2big_1D_4D
!!***

!!****f* m_abi2big/wvl_rhov_abi2big_2D_2D
!! NAME
!!  wvl_rhov_abi2big_2D_2D
!!
!! FUNCTION
!!  Copies a density/potential from ABINIT to BigDFT or viceversa.
!!  Target : ABINIT 2D arrays (with spin), BigDFT 2D arrays (with spin)
!!
!! INPUTS
!!  opt= 1: copy from ABINIT to BigDFT, 2: copy from BigDFT to ABINIT
!!  rhov_abi(:,:) = density/potential array in ABINIT
!!  rhov_big(:,:) = density/potential array in BigDFT
!!  [shift] = shift to be applied in rhov_abi array (parallelism)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  At output rhov_abi is copied to rhov_abinit, or viceversa
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_rhov_abi2big_2D_2D(opt,rhov_abi,rhov_big,shift)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_rhov_abi2big_2D_2D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: opt
 integer,intent(in),optional :: shift
 real(dp) :: rhov_abi(:,:),rhov_big(:,:)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer :: nfft_abi,nfft_big,nspden,shift_
 character(len=100) :: message
#endif
 
! *************************************************************************

#if defined HAVE_BIGDFT
 nspden=size(rhov_abi,2)
 if (size(rhov_big,2)/=nspden) then
   message='wvl_rhov_abi2big: ABINIT and BigDFT objects do not have the same nspden'
   MSG_BUG(message)
 end if
 nfft_abi=size(rhov_abi)/nspden ; nfft_big=size(rhov_big)/nspden
 shift_=0;if (present(shift)) shift_=shift
 call wvl_rhov_abi2big_gen(nfft_abi,nfft_big,nspden,opt,rhov_abi,rhov_big,shift_)
#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*)  opt,rhov_abi(1,1),rhov_big(1,1)
#endif

end subroutine wvl_rhov_abi2big_2D_2D
!!***

!!****f* m_abi2big/wvl_rhov_abi2big_1D_2D
!! NAME
!!  wvl_rhov_abi2big_1D_2D
!!
!! FUNCTION
!!  Copies a density/potential from ABINIT to BigDFT or viceversa.
!!  Target : ABINIT 1D arrays (without spin), BigDFT 2D arrays (with spin)
!!
!! INPUTS
!!  opt= 1: copy from ABINIT to BigDFT, 2: copy from BigDFT to ABINIT
!!  rhov_abi(:)   = density/potential array in ABINIT
!!  rhov_big(:,:) = density/potential array in BigDFT
!!  [shift] = shift to be applied in rhov_abi array (parallelism)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  At output rhov_abi is copied to rhov_abinit, or viceversa
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_rhov_abi2big_1D_2D(opt,rhov_abi,rhov_big,shift)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_rhov_abi2big_1D_2D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: opt
 integer,intent(in),optional :: shift
 real(dp) :: rhov_abi(:),rhov_big(:,:)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer :: nfft_abi,nfft_big,shift_
 character(len=100) :: message
#endif
 
! *************************************************************************

#if defined HAVE_BIGDFT
 if (size(rhov_big,2)/=1) then
   message='wvl_rhov_abi2big: ABINIT and BigDFT objects do not have the same nspden'
   MSG_BUG(message)
 end if
 nfft_abi=size(rhov_abi) ; nfft_big=size(rhov_big)
 shift_=0;if (present(shift)) shift_=shift
 call wvl_rhov_abi2big_gen(nfft_abi,nfft_big,1,opt,rhov_abi,rhov_big,shift_)
#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*)  opt,rhov_abi(1),rhov_big(1,1)
#endif

end subroutine wvl_rhov_abi2big_1D_2D
!!***

!!****f* m_abi2big/wvl_rhov_abi2big_2D_1D
!! NAME
!!  wvl_rhov_abi2big_2D_1D
!!
!! FUNCTION
!!  Copies a density/potential from ABINIT to BigDFT or viceversa.
!!  Target : ABINIT 2D arrays (with spin), BigDFT 1D arrays (with or without spin)
!!
!! INPUTS
!!  opt= 1: copy from ABINIT to BigDFT, 2: copy from BigDFT to ABINIT
!!  rhov_abi(:,:) = density/potential array in ABINIT
!!  rhov_big(:)   = density/potential array in BigDFT
!!  [shift] = shift to be applied in rhov_abi array (parallelism)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  At output rhov_abi is copied to rhov_abinit, or viceversa
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_rhov_abi2big_2D_1D(opt,rhov_abi,rhov_big,shift)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_rhov_abi2big_2D_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: opt
 integer,intent(in),optional :: shift
 real(dp) :: rhov_abi(:,:),rhov_big(:)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer :: nfft_abi,nfft_big,nspden,shift_
#endif
 
! *************************************************************************

#if defined HAVE_BIGDFT
 nspden=size(rhov_abi,2)
 nfft_abi=size(rhov_abi)/nspden ; nfft_big=size(rhov_big)/nspden
 shift_=0;if (present(shift)) shift_=shift
 call wvl_rhov_abi2big_gen(nfft_abi,nfft_big,nspden,opt,rhov_abi,rhov_big,shift_)
#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*)  opt,rhov_abi(1,1),rhov_big(1)
#endif

end subroutine wvl_rhov_abi2big_2D_1D
!!***

!!****f* m_abi2big/wvl_rhov_abi2big_1D_1D
!! NAME
!!  wvl_rhov_abi2big_1D_1D
!!
!! FUNCTION
!!  Copies a density/potential from ABINIT to BigDFT or viceversa.
!!  Target : ABINIT 2D arrays (without spin), BigDFT 1D arrays (with or without spin)
!!
!! INPUTS
!!  opt= 1: copy from ABINIT to BigDFT, 2: copy from BigDFT to ABINIT
!!  rhov_abi(:) = density/potential array in ABINIT
!!  rhov_big(:) = density/potential array in BigDFT
!!  [shift] = shift to be applied in rhov_abi array (parallelism)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  At output rhov_abi is copied to rhov_abinit, or viceversa
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_rhov_abi2big_1D_1D(opt,rhov_abi,rhov_big,shift)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_rhov_abi2big_1D_1D'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: opt
 integer,intent(in),optional :: shift
 real(dp) :: rhov_abi(:),rhov_big(:)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer :: nfft_abi,nfft_big,shift_
#endif
 
! *************************************************************************

#if defined HAVE_BIGDFT
 nfft_abi=size(rhov_abi) ; nfft_big=size(rhov_big)
 shift_=0;if (present(shift)) shift_=shift
 call wvl_rhov_abi2big_gen(nfft_abi,nfft_big,1,opt,rhov_abi,rhov_big,shift_)
#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*)  opt,rhov_abi(1),rhov_big(1)
#endif

end subroutine wvl_rhov_abi2big_1D_1D
!!***

!!****f* m_abi2big/wvl_rhov_abi2big_gen
!! NAME
!!  wvl_rhov_abi2big_gen
!!
!! FUNCTION
!!  Copies a density/potential from ABINIT to BigDFT or viceversa.
!!  This is a generic routine to copy objects.
!!
!! INPUTS
!!  nfft_abi = size of rhov_abi
!!  nfft_big = size of rhov_big
!!  nspden =  number of spin components
!!  opt= 1) copy from ABINIT to BigDFT
!!       2) copy from BigDFT to ABINIT
!!  rhov_abi(nfft_abi,nspden) = density/potential array in ABINIT
!!  rhov_big(nfft_big,nspden) = density/potential array in BigDFT
!!  shift = shift to be applied in rhov_abi array (parallelism)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! At output rhov_abi is copied to rhov_abinit, or viceversa
!!
!! NOTES
!! This routine is duplicated:
!! This option is faster but it requires more memory.
!! Notice that we cannot point the variables since the spin convention is not
!! the same in BigDFT and ABINIT.
!! In ABINIT: index 1 is for the total spin (spin up + spin down) and index 2 is for spin up.
!! In BigDFT: indices 1 and 2 are for spin up and down, respectively.
!!
!! PARENTS
!!      m_abi2big
!!
!! CHILDREN
!!
!! SOURCE

subroutine wvl_rhov_abi2big_gen(nfft_abi,nfft_big,nspden,opt,rhov_abi,rhov_big,shift)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_rhov_abi2big_gen'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: nfft_abi,nfft_big,nspden,opt,shift
 real(dp) :: rhov_abi(nfft_abi,nspden),rhov_big(nfft_big,nspden)

!Local variables-------------------------------
#if defined HAVE_BIGDFT
 integer :: ifft,jfft,nfft
 real(dp) :: tmpUp,tmpDown,tmpTot
 character(len=100) :: message
 real(dp),allocatable :: rhoup(:),rhodn(:),rhotot(:)
#endif
 
! *************************************************************************

 DBG_ENTER("COLL")

#if defined HAVE_BIGDFT
!No objects to copy; in BigDFT by default they have size of 1!
 if(size(rhov_big)==1.and.size(rhov_abi)==0) return

 nfft=nfft_big;if (nfft_big+shift>nfft_abi) nfft=nfft-shift

 if (nfft_abi<nfft+shift) then
   message='wvl_rhov_abi2big: cannot handle nfft(abi)<nfft(big)+shift case'
   MSG_BUG(message)
 end if
 if(nspden==4) then
   message='wvl_rhov_abi2big: nspden=4 not yet supported'
   MSG_ERROR(message)
 end if

 if (hmem.and.nspden==2) then
   if (opt==1) then
     ABI_ALLOCATE(rhoup,(nfft))
     ABI_ALLOCATE(rhodn,(nfft))
   else if (opt==2)  then
     ABI_ALLOCATE(rhotot,(nfft))
   end if
 end if
 
 if (opt==1) then !ABINIT -> BIGDFT
   if (nspden==2) then
     if (hmem) then
       rhoup(1:nfft)=rhov_abi(shift+1:shift+nfft,2)
       rhodn(1:nfft)=rhov_abi(shift+1:shift+nfft,1)-rhoup(1:nfft)
       rhov_big(:,1)=rhoup(:)
       rhov_big(:,2)=rhodn(:)
     else
       do ifft=1,nfft
         jfft=shift+ifft
!        We change convention for BigDFT
         tmpDown=rhov_abi(jfft,1)-rhov_abi(jfft,2)
         tmpUp  =rhov_abi(jfft,2)
         rhov_big(ifft,1)=tmpUp
         rhov_big(ifft,2)=tmpDown
       end do
     end if !hmem
   else !nspden==1
     if (hmem) then
       rhov_big(1:nfft,1)=rhov_abi(shift+1:shift+nfft,1)
     else
       do ifft=1,nfft
         rhov_big(ifft,1)=rhov_abi(shift+ifft,1)
       end do
     end if!hmem
   end if !nspden

 else if (opt==2) then !BigDFT -> ABINIT
   if (nspden==2) then
     if (hmem) then
       rhotot(:)=rhov_big(:,1)+rhov_big(:,2)
       rhov_abi(shift+1:shift+nfft,1)=rhotot(1:nfft)
       rhov_abi(shift+1:shift+nfft,2)=rhov_big(1:nfft,1)
     else 
       do ifft=1,nfft
         jfft=shift+ifft
!        We change convention for BigDFT
         tmpTot=rhov_big(ifft,1)+rhov_big(ifft,2)
         rhov_abi(jfft,1)=tmpTot
         rhov_abi(jfft,2)=rhov_big(ifft,1) !Spin Up
       end do
     end if !hmem
   else if (nspden==1) then
     if (hmem) then
       rhov_abi(shift+1:shift+nfft,1)=rhov_big(1:nfft,1)
     else
       do ifft=1,nfft
         rhov_abi(shift+ifft,1)=rhov_big(ifft,1)
       end do
     end if !hmem
   end if !nspden

 else
   message='wvl_rhov_abi2big_gen: wrong option'
   MSG_BUG(message)
 end if

 if (hmem.and.nspden==2) then
   if (opt==1) then
     ABI_DEALLOCATE(rhoup)
     ABI_DEALLOCATE(rhodn)
   else if (opt==2) then
     ABI_DEALLOCATE(rhotot)
   end if !opt
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*)  nfft_abi,nfft_big,nspden,opt,shift,rhov_big(1,1),rhov_abi(1,1)
#endif

 DBG_EXIT("COLL")

end subroutine wvl_rhov_abi2big_gen
!!***

end module m_abi2big
!!***

