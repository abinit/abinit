!{\src2tex{textfont=tt}}
!!****f* ABINIT/normsq_gkq
!!
!! NAME
!! normsq_gkq
!!
!! FUNCTION
!! This routine takes the gkq matrix elements for a given qpoint,
!!   does the scalar product with the phonon displacement vector,
!!   squares the gkq matrix elements multiplies by the appropriate weights 
!!   and puts them in a uniform (atom,icart) basis
!!
!! COPYRIGHT
!! Copyright (C) 2004-2016 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   displ_red = phonon mode displacement vectors in reduced coordinated.
!!   eigvec = eigenvectors of phonons (to turn to cartesian coord frame)
!!   elph_ds = datastructure with gkk matrix elements
!!   FSfullpqtofull = mapping of k + q to k
!!   h1_mat_el_sq = matrix elements $<psi_{k+q,m} | H^{1} | psi_{k,n}>$ matrix-squared
!!   iqptirred = index of present qpoint
!!   phfrq_tmp = phonon frequencies
!!   qpt_irred = array of qpoint coordinates
!!
!! OUTPUT
!!   elph_ds%gkq filled
!!   qdata(elph_ds%nbranch,elph_ds%nsppol,3) = array containing the phonon frequency, the linewidth and $\lambda_{q,\nu}$.
!!
!! PARENTS
!!      read_gkk
!!
!! CHILDREN
!!      gam_mult_displ,nmsq_gam,nmsq_gam_sumfs,nmsq_pure_gkk
!!      nmsq_pure_gkk_sumfs,wrtout,xmpi_sum,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine normsq_gkq(displ_red,eigvec,elph_ds,FSfullpqtofull,&
&    h1_mat_el_sq,iqptirred,phfrq_tmp,qpt_irred,qdata)

 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_errors
 use m_xmpi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'normsq_gkq'
 use interfaces_14_hidewrite
 use interfaces_77_ddb, except_this_one => normsq_gkq
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptirred
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(inout) :: &
& h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol)
 real(dp),intent(in) :: phfrq_tmp(elph_ds%nbranch),qpt_irred(3,elph_ds%nqptirred)
 real(dp),intent(out) :: qdata(elph_ds%nbranch,elph_ds%nsppol,3)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,ier,ii,isppol,jbranch,comm
 real(dp) :: lambda_tot
 character(len=500) :: message
!arrays
 real(dp) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp) :: gam_now2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: lambda(elph_ds%nsppol)
 real(dp),allocatable :: matrx(:,:),val(:),vec(:,:,:)
 real(dp),allocatable :: zhpev1(:,:),zhpev2(:)

! *************************************************************************

 DBG_ENTER("COLL")

 accum_mat  = zero
 accum_mat2 = zero
 comm = xmpi_world

 if (elph_ds%ep_scalprod == 1) then
!  
   if (elph_ds%ep_keepbands == 0) then
     call wrtout(std_out,' normsq_gkq : calling nmsq_gam_sumFS',"COLL")
     call nmsq_gam_sumFS (accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&
&     h1_mat_el_sq,iqptirred)

   else if (elph_ds%ep_keepbands == 1) then
     call wrtout(std_out,' normsq_gkq : calling nmsq_gam',"COLL")
     call nmsq_gam (accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&
&     h1_mat_el_sq,iqptirred)

   else
     write (message,'(a,i0)')' Wrong value for elph_ds%ep_keepbands = ',elph_ds%ep_keepbands
     MSG_BUG(message)
   end if
!  
 else if (elph_ds%ep_scalprod == 0) then  ! Interpolate on the pure "matrix of matrix elements" and do the scalar products later.
!  
   if (elph_ds%ep_keepbands == 0) then
     call wrtout(std_out,' normsq_gkq : calling nmsq_pure_gkk_sumFS',"COLL")
     call nmsq_pure_gkk_sumFS (accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,&
&     h1_mat_el_sq,iqptirred)

   else if (elph_ds%ep_keepbands == 1) then
     call wrtout(std_out,' normsq_gkq : calling nmsq_pure_gkk',"COLL")

     call nmsq_pure_gkk (accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,&
&     h1_mat_el_sq,iqptirred)
   else
     write (message,'(a,i0)')' Wrong value for elph_ds%ep_keepbands = ',elph_ds%ep_keepbands
     MSG_BUG(message)
   end if
!  
 else
   write (message,'(a,i0)')' Wrong value for elph_ds%ep_scalprod = ',elph_ds%ep_scalprod
   MSG_BUG(message)
 end if
!end if flag for doing scalar product now.


!MG: values without the good prefactor
 accum_mat = accum_mat * elph_ds%occ_factor/elph_ds%k_phon%nkpt

!MG: accum_mat2 contains the line-widhts before the Fourier interpolation
 accum_mat2 = accum_mat2 * elph_ds%occ_factor/elph_ds%k_phon%nkpt

!mpi sum over procs for accum_mat2
 call xmpi_sum (accum_mat, comm, ier)
 call xmpi_sum (accum_mat2, comm, ier)

!MG20060531i
!write e-ph quantities before Fourier interpolation
!save e-ph values in the temporary array qdata that will be copied into elph_ds%qgrid_data

 write (message,'(4a,3es16.6,63a)')ch10,                  &
& ' Phonon linewidths before interpolation ',ch10,        &
& ' Q point = ',qpt_irred(:,iqptirred),ch10,('=',ii=1,60),ch10,&
& ' Mode          Frequency (Ha)  Linewidth (Ha)  Lambda '
 call wrtout(std_out,message,'COLL')

 lambda_tot = zero
 do isppol=1,elph_ds%nsppol
   do ii=1,elph_ds%nbranch
     lambda(isppol)=zero
!    MG: the tolerance factor is somehow arbitrary
     if (abs(phfrq_tmp(ii)) > tol10) lambda(isppol)=accum_mat2(1,ii,ii,isppol)/&
&     (pi*elph_ds%n0(isppol)*phfrq_tmp(ii)**2)
     lambda_tot=lambda_tot+lambda(isppol)
     write(message,'(i8,es20.6,2es16.6)' )ii,phfrq_tmp(ii),accum_mat2(1,ii,ii,isppol),lambda(isppol)
     call wrtout(std_out,message,'COLL')
!    save values
     qdata(ii,isppol,1)=phfrq_tmp(ii)
     qdata(ii,isppol,2)=accum_mat2(1,ii,ii,isppol)
     qdata(ii,isppol,3)=lambda(isppol)
   end do !loop over branch
 end do !loop over sppol

!normalize for number of spins
 lambda_tot = lambda_tot / elph_ds%nsppol

 write(message,'(61a,44x,es16.6,62a)' )('=',ii=1,60),ch10,lambda_tot,ch10,('=',ii=1,60),ch10
 call wrtout(std_out,message,'COLL')
!ENDMG20060531

!immediately calculate linewidths:
 write(std_out,*) 'summed accum_mat = '
 write(std_out,'(3(2E18.6,1x))') accum_mat(:,:,:,1)
 write(std_out,*) 'summed accum_mat2 = '
 write(std_out,'(3(2E18.6,1x))')  (accum_mat2(:,ii,ii,1),ii=1,elph_ds%nbranch)
 write(std_out,*) 'displ_red  = '
 write(std_out,'(3(2E18.6,1x))') displ_red

 if (elph_ds%ep_scalprod == 1) then
   do isppol=1,elph_ds%nsppol
!    Diagonalize gamma matrix at qpoint (complex matrix). Copied from dfpt_phfrq
     ier=0
     ii=1
     ABI_ALLOCATE(matrx,(2,(elph_ds%nbranch*(elph_ds%nbranch+1))/2))
     do i2=1,elph_ds%nbranch
       do i1=1,i2
         matrx(1,ii)=accum_mat2(1,i1,i2,isppol)
         matrx(2,ii)=accum_mat2(2,i1,i2,isppol)
         ii=ii+1
       end do
     end do
     ABI_ALLOCATE(zhpev1,(2,2*elph_ds%nbranch-1))
     ABI_ALLOCATE(zhpev2,(3*elph_ds%nbranch-2))
     ABI_ALLOCATE(val,(elph_ds%nbranch))
     ABI_ALLOCATE(vec,(2,elph_ds%nbranch,elph_ds%nbranch))
     call ZHPEV ('V','U',elph_ds%nbranch,matrx,val,vec,elph_ds%nbranch,zhpev1,zhpev2,ier)

     write (std_out,*) ' normsq_gkq : accumulated eigenvalues isppol ',isppol, ' = '
     write (std_out,'(3E18.6)') val
     ABI_DEALLOCATE(matrx)
     ABI_DEALLOCATE(zhpev1)
     ABI_DEALLOCATE(zhpev2)
     ABI_DEALLOCATE(vec)
     ABI_DEALLOCATE(val)
   end do ! isppol

 else if (elph_ds%ep_scalprod == 0) then


   do isppol=1,elph_ds%nsppol
     call gam_mult_displ(elph_ds%nbranch, displ_red, accum_mat(:,:,:,isppol), gam_now2)

     write (std_out,*) ' normsq_gkq : accumulated eigenvalues isppol ', isppol, ' = '
     write (std_out,'(3(E14.6,1x))') (gam_now2(1,jbranch,jbranch), jbranch=1,elph_ds%nbranch)
     write (std_out,*) ' normsq_gkq : imag part = '
     write (std_out,'(3(E14.6,1x))') (gam_now2(2,jbranch,jbranch), jbranch=1,elph_ds%nbranch)
   end do ! isppol

 end if

 DBG_EXIT("COLL")

end subroutine normsq_gkq
!!***

