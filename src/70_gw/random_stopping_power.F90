!{\src2tex{textfont=tt}}
!!****f* ABINIT/random_stopping_power
!! NAME
!! random_stopping_power
!!
!! FUNCTION
!! Calculate the electronic random stopping power
!!
!! COPYRIGHT
!! Copyright (C) 2001-2017 ABINIT group (AS,FB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!    
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!      get_bz_item,spline,splint,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine random_stopping_power(iqibz,npvel,pvelmax,Ep,Gsph_epsG0,Qmesh,Vcp,Cryst,Dtfil,epsm1,rspower)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_splines          

 use m_io_tools,      only : open_file
 use m_bz_mesh,       only : get_BZ_item, kmesh_t
 use m_gsphere,       only : gsphere_t
 use m_crystal,       only : crystal_t
 use m_gwdefs,        only : GW_TOLQ0, em1params_t
 use m_vcoul,         only : vcoul_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'random_stopping_power'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit  none

!Arguments ------------------------------------
!scalars
 integer,intent(in)                    :: iqibz,npvel
 real(dp),intent(in)                   :: pvelmax(3)

 type(em1params_t),intent(in) :: Ep
 type(gsphere_t),intent(in)            :: Gsph_epsG0
 type(kmesh_t),intent(in)              :: Qmesh
 type(vcoul_t),intent(in)              :: Vcp
 type(crystal_t),intent(in)            :: Cryst
 type(Datafiles_type),intent(in)       :: Dtfil

 complex(gwpc),intent(in)              :: epsm1(Ep%npwe,Ep%npwe,Ep%nomega)     

 real(dp),intent(inout)                :: rspower(npvel) 

!Local variables ------------------------------
 integer :: ipvel,ig
 integer :: iq_bz,iq_ibz,isym_q,itim_q
 integer :: iomega,iomegap,nomega_re
 integer :: unt_rsp
 integer,allocatable :: iomega_re(:)

 real(dp),parameter :: zp=1.0_dp              ! Hard-coded charge of the impinging particle
 real(dp) :: omega_p
 real(dp) :: im_epsm1_int(1) 
 real(dp) :: qbz(3),qpgcart(3),qpgred(3)
 real(dp) :: pvel(3,npvel),pvel_norm(npvel)
 real(dp) :: ypp_i(Ep%nomega) 
 real(dp) :: vcoul(Ep%npwe)
 real(dp),allocatable :: im_epsm1_diag_qbz(:,:),tmp_data(:)
 real(dp),allocatable :: omega_re(:)

 character(len=500)     :: msg
 character(len=fnlen+4) :: fname

!************************************************************************

 !
 ! First set up the velocities array from the input variables npvel and pvelmax(3)
 ! Remember pvelmax is in Cartesian coordinates and so is pvel
 do ipvel=1,npvel
   pvel(:,ipvel)    = REAL(ipvel,dp) / REAL(npvel,dp) * pvelmax(:)
   pvel_norm(ipvel) = SQRT( SUM( pvel(:,ipvel)**2 ) )
 enddo
 !
 ! Select the purely real frequency in Ep%omega
 nomega_re=0
 do iomega=1,Ep%nomega
   if( AIMAG(Ep%omega(iomega)) < 1.0e-4_dp ) then
     nomega_re=nomega_re+1
   endif
 enddo
 ABI_ALLOCATE(omega_re,(nomega_re))
 ABI_ALLOCATE(iomega_re,(nomega_re))
 ABI_ALLOCATE(im_epsm1_diag_qbz,(Ep%npwe,Ep%nomega))
 ABI_ALLOCATE(tmp_data,(Ep%nomega))

 iomegap=0
 do iomega=1,Ep%nomega
   if( AIMAG(Ep%omega(iomega)) < 1.0e-4_dp ) then
     iomegap=iomegap+1
     iomega_re(iomegap)=iomega
     omega_re(iomegap)=REAL(Ep%omega(iomega),dp)
   endif
 enddo

 !
 ! Loop over all the q-points in the full Brillouin zone and select only the
 ! ones that corresponds to the correct q-point in the irreducible wedge we are
 ! currently treating (index iqibz)
 do iq_bz=1,Qmesh%nbz

   ! Perform the check and obtain the symmetry information
   call get_BZ_item(Qmesh,iq_bz,qbz,iq_ibz,isym_q,itim_q)
   if( iqibz /= iq_ibz ) cycle

   ! Apply the symmetry operation to the diagonal of epsm1 
   do iomega=1,nomega_re
     do ig=1,Ep%npwe
       im_epsm1_diag_qbz(Gsph_epsG0%rottb(ig,itim_q,isym_q),iomega)= AIMAG( epsm1(ig,ig,iomega_re(iomega)) )
     enddo
   enddo
   ! Apply the symmetry operation to the Coulomb interaction
   do ig=1,Ep%npwe
     vcoul(Gsph_epsG0%rottb(ig,itim_q,isym_q))=Vcp%vc_sqrt(ig,iqibz)**2 
   enddo

   !
   ! Sum over G vectors
   do ig=1,Ep%npwe
     !
     ! Loop over velocities
     do ipvel=1,npvel

       qpgred(:) = qbz(:) + Gsph_epsG0%gvec(:,ig)
       ! Transform q + G from reduced to cartesian with the symmetry operation
       qpgcart(:) = two_pi * Cryst%gprimd(:,1) * qpgred(1) & 
&                 + two_pi * Cryst%gprimd(:,2) * qpgred(2) & 
&                 + two_pi * Cryst%gprimd(:,3) * qpgred(3)   

       ! omega_p = ( q + G ) . v
       omega_p =  DOT_PRODUCT( qpgcart(:) , pvel(:,ipvel) )
       
       ! Check that the calculated frequency omega_p is within the omega
       ! range of epsm1 and thus that the interpolation will go fine
       if ( ABS(omega_p) > omega_re(nomega_re) ) then
         write(msg,'(a,e16.4,2a,e16.4)') ' freqremax is currently ',omega_re(nomega_re),ch10,&
&                                        ' increase it to at least ',omega_p
         MSG_WARNING(msg)
       endif
       
       ! Perform the spline interpolation to obtain epsm1 at the desired omega = omega_p
       tmp_data = im_epsm1_diag_qbz(ig,:)

       call spline( omega_re, tmp_data, nomega_re, 1.0e+32_dp, 1.0e+32_dp, ypp_i)
       call splint( nomega_re, omega_re, tmp_data, ypp_i, 1, (/ ABS(omega_p) /),  im_epsm1_int )

       !
       ! Apply the odd parity of Im epsm1 in  omega to recover the causal response function
       if (omega_p<zero) im_epsm1_int(1)=-im_epsm1_int(1)

       !
       ! Calculate 4 * pi / |q+G|**2 * omega_p * Im{ epsm1_GG(q,omega_p) }
       !
       im_epsm1_int(1) = omega_p * vcoul(ig) * im_epsm1_int(1)

       ! Accumulate the final result without the prefactor
       ! (It will be included at the very end)
       rspower(ipvel) = rspower(ipvel) + im_epsm1_int(1)

     end do ! end of velocity loop 
   end do ! end G-loop 

 enddo ! end of q loop in the full BZ

 ! If it is the last q, write down the result in the main output file and in a
 ! separate file _RSP (for Random Stopping Power)
 if (iqibz == Qmesh%nibz ) then

   ! Multiply by the prefactors
   ! Note that this expression differs from Eq. (3.11) in Campillo PRB 58, 10309 (1998).
   ! A factor one half is missing in the paper.
   rspower(:) = - zp**2 / ( Cryst%ucvol * Qmesh%nbz * pvel_norm(:) ) * rspower(:)

   write(msg,'(2a)')         ch10,' ==== Random stopping power along Cartesian direction  === '
   call wrtout(std_out,msg,'COLL') 
   call wrtout(ab_out,msg,'COLL') 
   write(msg,'(a,3(f12.4,2x),a)') ' ====  ',pvelmax(:),'===='
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL') 
   write(msg,'(a)')               '#  |v| (a.u.) , RSP (a.u.) '
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL') 
   do ipvel=1,npvel
     write(msg,'(f16.8,4x,f16.8)') pvel_norm(ipvel),rspower(ipvel)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL') 
   enddo
   write(msg,'(2a)')              ' ========================================================= ',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL') 

   fname=TRIM(Dtfil%filnam_ds(4))//'_RSP'
   if (open_file(fname,msg,newunit=unt_rsp,status='unknown',form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if

   write(msg,'(a)')               '# ==== Random stopping power along Cartesian direction  === '
   call wrtout(unt_rsp,msg,'COLL')
   write(msg,'(a,3(f12.4,2x))')   '# ====  ',pvelmax(:)
   call wrtout(unt_rsp,msg,'COLL')
   write(msg,'(a)')               '#  |v| (a.u.) , RSP (a.u.) '
   call wrtout(unt_rsp,msg,'COLL')
   do ipvel=1,npvel
     write(msg,'(f16.8,4x,f16.8)') pvel_norm(ipvel),rspower(ipvel)
     call wrtout(unt_rsp,msg,'COLL')
   enddo
   close(unt_rsp)
 end if

 ABI_DEALLOCATE(omega_re)
 ABI_DEALLOCATE(iomega_re)
 ABI_DEALLOCATE(im_epsm1_diag_qbz)
 ABI_DEALLOCATE(tmp_data)

end subroutine random_stopping_power
!!***
