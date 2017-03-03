!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkph_linwid
!!
!! NAME
!! mkph_linwid
!!
!! FUNCTION
!!  Calculate the phonon linewidths on a trajectory in q space
!!
!! COPYRIGHT
!! Copyright (C) 2004-2017 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  Cryst<crystal_t>=Info on the unit cell and symmetries.
!!  Ifc<ifc_type>=Object containing the interatomic force constants.
!!  elph_ds = datastructure with phonon matrix elements
!!  nqpath = dimension of qpath_vertices
!!  qpath_vertices = vertices of reciprocal space trajectory
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      ftgam,ftgam_init,gam_mult_displ,ifc_fourq,make_path,phdispl_cart2red
!!      wrap2_pmhalf,wrtout,zgemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mkph_linwid(Cryst,ifc,elph_ds,nqpath,qpath_vertices)

 use defs_basis
 use defs_elphon
 use m_profiling_abi
 use m_errors
 use m_bz_mesh
 use m_ifc

 use m_io_tools,        only : open_file
 use m_geometry,        only : phdispl_cart2red
 use m_dynmat,          only : ftgam_init, ftgam
 use m_crystal,         only : crystal_t
 use m_numeric_tools,   only : wrap2_pmhalf

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkph_linwid'
 use interfaces_14_hidewrite
 use interfaces_77_ddb, except_this_one => mkph_linwid
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqpath
 type(crystal_t),intent(in) :: Cryst
 type(ifc_type),intent(in) :: ifc
 type(elph_type),intent(inout) :: elph_ds
!arrays
 real(dp),intent(in) :: qpath_vertices(3,nqpath)

!Local variables-------------------------------
!scalars
 integer :: ibranch,natom,ii,indx,ipoint,nbranch,nqbz,nsppol,nrpt
 integer :: isppol,jbranch,qtor,unit_bs,unit_lambda,unit_lwd,npt_tot
 real(dp) :: diagerr,res
 character(len=500) :: msg
 character(len=fnlen) :: fname,base_name
!arrays
 integer :: ndiv(nqpath-1)
 integer, allocatable :: indxprtqpt(:)
 real(dp),parameter :: c0(2)=(/0._dp,0._dp/),c1(2)=(/1._dp,0._dp/)
 real(dp) :: displ_cart(2,3*Cryst%natom,3*Cryst%natom)
 real(dp) :: displ_red(2,3*Cryst%natom,3*Cryst%natom)
 real(dp) :: eigval(3*Cryst%natom)
 real(dp) :: gam_now(2,(3*Cryst%natom)**2)
 real(dp) :: imeigval(3*Cryst%natom)
 real(dp) :: lambda(3*Cryst%natom)
 real(dp) :: pheigvec(2*3*Cryst%natom*3*Cryst%natom),phfrq_tmp(3*Cryst%natom)
 real(dp) :: qpt(3),redkpt(3)
 real(dp) :: tmpgam1(2,3*Cryst%natom,3*Cryst%natom)
 real(dp) :: tmpgam2(2,3*Cryst%natom,3*Cryst%natom)
 real(dp), allocatable :: coskr(:,:), sinkr(:,:),finepath(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 natom     = Cryst%natom
 nbranch   = elph_ds%nbranch 
 nsppol    = elph_ds%nsppol
 base_name = elph_ds%elph_base_name
 nrpt = ifc%nrpt

!===================================================================
!Definition of the q path along which ph linwid will be interpolated
!===================================================================
 call make_path(nqpath,qpath_vertices,Cryst%gmet,'G',20,ndiv,npt_tot,finepath)
 ABI_ALLOCATE(indxprtqpt,(npt_tot))
 indxprtqpt = 0

!==========================================================
!Open _LWD file and write header
!==========================================================
 fname=trim(base_name) // '_LWD'
 if (open_file(fname,msg,newunit=unit_lwd,status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if 

 write (unit_lwd,'(a)')       '#'
 write (unit_lwd,'(a)')       '# ABINIT package : Phonon linewidth file'
 write (unit_lwd,'(a)')       '#'
 write (unit_lwd,'(a,i10,a)') '#  Phonon linewidths calculated on ',npt_tot,' points along the qpath'
 write (unit_lwd,'(a)')       '#  Description of the Q-path :'
 write (unit_lwd, '(a,i10)')  '#  Number of line segments = ',nqpath-1
 write (unit_lwd,'(a)')       '#  Vertices of the Q-path and corresponding index = '

 indx=1
 indxprtqpt(1) = 1
 indxprtqpt(npt_tot) = 1

 do ii=1,nqpath
   write (unit_lwd,'(a,3(e16.6,1x),i8)')'#  ',qpath_vertices(:,ii),indx
   if (ii<nqpath) then
     indx=indx+ndiv(ii)
     indxprtqpt(indx) = 1
   end if
 end do

 write (unit_lwd,'(a)')'#'

!==========================================================
!Open _BST file and write header
!==========================================================
 fname=trim(base_name) // '_BST'
 if (open_file(fname,msg,newunit=unit_bs,status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if 

 write (unit_bs, '(a)')      '#'
 write (unit_bs, '(a)')      '# ABINIT package : Phonon band structure file'
 write (unit_bs, '(a)')      '#'
 write (unit_bs, '(a,i10,a)')'# Phonon BS calculated on ', npt_tot,' points along the qpath'
 write (unit_bs, '(a,i10)')  '# Number of line segments = ', nqpath-1
 indx=1
 do ii=1,nqpath
   write (unit_bs,'(a,3(E16.6,1x),i8)')'#  ',qpath_vertices(:,ii),indx
   if (ii<nqpath) indx=indx+ndiv(ii)
 end do
 write (unit_bs,'(a)')'#'

!MG20060606
!==========================================================
!open _LAMBDA file and write header
!contains \omega(q,n) and \lambda(q,n) and can be plotted using xmgrace
!==========================================================
 fname=trim(base_name) // '_LAMBDA'
 if (open_file(fname,msg,newunit=unit_lambda,status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if 

 write (unit_lambda,'(a)')      '#'
 write (unit_lambda,'(a)')      '# ABINIT package : Lambda file'
 write (unit_lambda,'(a)')      '#'
 write (unit_lambda,'(a,i10,a)')'#  Lambda(q,nu) calculated on ',npt_tot,' Q-points'
 write (unit_lambda,'(a)')      '# Description of the Q-path :'
 write (unit_lambda,'(a,i10)')  '# Number of line segments = ',nqpath-1
 write (unit_lambda,'(a)')      '# Vertices of the Q-path and corresponding index = '

 indx=1
 do ii=1,nqpath
   write (unit_lambda,'(a,3(E16.6,1x),i8)')'#  ',qpath_vertices(:,ii),indx
   if (ii<nqpath) indx=indx+ndiv(ii)
 end do
 write (unit_lambda,'(a)')'#'
 write (unit_lambda,'(a)')'# index frequency lambda(q,n) frequency lambda(q,n) .... lambda_tot'
 write (unit_lambda,'(a)')'#'

!real space to q space
 qtor=0

!initialize the maximum phonon frequency
 elph_ds%omega_min = zero
 elph_ds%omega_max = zero

 ABI_ALLOCATE(coskr, (npt_tot,nrpt))
 ABI_ALLOCATE(sinkr, (npt_tot,nrpt))
 call ftgam_init(ifc%gprim, npt_tot, nrpt, finepath, ifc%rpt, coskr, sinkr)

 write (std_out,*) ' mkph_linwid : shape(elph_ds%gamma_qpt) = ',shape(elph_ds%gamma_qpt)
 nqbz =  SIZE(elph_ds%gamma_qpt,DIM=4)
 write(std_out,*) " nqbz =  SIZE(elph_ds%gamma_qpt,DIM=4) = ",nqbz
!
!Big do loop over spin polarizations
!could put in locally, so phonon stuff is not done twice...
!
 do isppol=1,nsppol
   indx=1

!  Output to the main output file
   write(msg,'(a,a)')ch10,&
&   ' Output of the linewidths for the first point of each segment. Linewidths are given in Hartree.'
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write (std_out,*) ' mkph_linwid : elph_ds%ep_scalprod = ', elph_ds%ep_scalprod

   qtor = 0

!  Interpolation along specified path in q space
   do ipoint=1,npt_tot

!    Get qpoint along the path from qpath_vertices
     qpt(:) = finepath(:,ipoint)

     call wrap2_pmhalf(qpt(1),redkpt(1),res)
     call wrap2_pmhalf(qpt(2),redkpt(2),res)
     call wrap2_pmhalf(qpt(3),redkpt(3),res)
     qpt(:) = redkpt(:)
!    
!    This reduced version of ftgkk supposes the kpoints have been integrated
!    in integrate_gamma. Do FT from real-space gamma grid to 1 qpt.
     call ftgam(ifc%wghatm,gam_now,elph_ds%gamma_rpt(:,:,isppol,:),natom,1,ifc%nrpt,qtor, &
&     coskr(ipoint,:), sinkr(ipoint,:))
!    
!    get phonon freqs and eigenvectors anyway
!    
     call ifc_fourq(ifc,cryst,qpt,phfrq_tmp,displ_cart,out_eigvec=pheigvec)
!    
!    additional frequency factor for some cases
!    
!    If the matrices do not contain the scalar product with the displ_cart vectors yet do it now
     if (elph_ds%ep_scalprod == 0) then

       call phdispl_cart2red(natom,Cryst%gprimd,displ_cart,displ_red)

       tmpgam2 = reshape (gam_now, (/2,nbranch,nbranch/))
       call gam_mult_displ(nbranch, displ_red, tmpgam2, tmpgam1)

       do jbranch=1,nbranch
         eigval(jbranch) = tmpgam1(1, jbranch, jbranch)
         imeigval(jbranch) = tmpgam1(2, jbranch, jbranch)

         if (abs(imeigval(jbranch)) > tol8) then
           write (msg,'(a,i0,a,es16.8)')' imaginary values for branch = ',jbranch,' imeigval = ',imeigval(jbranch)
           MSG_WARNING(msg)
         end if
       end do

     else if (elph_ds%ep_scalprod == 1) then
!      
!      Diagonalize gamma matrix at qpoint (complex matrix).
!      MJV NOTE: gam_now is recast implicitly here to matrix 
       call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, c1, gam_now, 3*natom,&
&       pheigvec, 3*natom, c0, tmpgam1, 3*natom)

       call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, c1, pheigvec, 3*natom,&
&       tmpgam1, 3*natom, c0, tmpgam2, 3*natom)

       diagerr = zero
       do ibranch=1,nbranch

         eigval(ibranch) = tmpgam2(1,ibranch,ibranch)

         do jbranch=1,ibranch-1
           diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))+abs(tmpgam2(2,jbranch,ibranch))
         end do
         do jbranch=ibranch+1,nbranch
           diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))+abs(tmpgam2(2,jbranch,ibranch))
         end do
         diagerr = diagerr + abs(tmpgam2(2,ibranch,ibranch))
       end do

       if (diagerr > tol12) then
         write (msg,'(a,es14.6)')' Numerical error in diagonalization of gamma with phon eigenvectors: ', diagerr
         MSG_WARNING(msg)
       end if

     else
       write (msg,'(a,i0)')' Wrong value for elph_ds%ep_scalprod = ',elph_ds%ep_scalprod
       MSG_BUG(msg)
     end if ! end elph_ds%ep_scalprod if
!    
!    ==========================================================
!    write data to files for each q point
!    ==========================================================
     write (unit_lwd,'(i5)', advance='no') indx
     do ii=1, nbranch
       write (unit_lwd,'(E16.5)',advance='no') eigval(ii)
     end do
     write (unit_lwd,*)

!    only print phonon BS for isppol 1: independent of electron spins
     if (isppol==1) then
       write (unit_bs,'(i5)', advance='no') indx
       do ii=1, nbranch
         write (unit_bs,'(E16.5)',advance='no') phfrq_tmp(ii)
       end do
       write (unit_bs,*)
     end if

     write (unit_lambda,'(i5)', advance='no') indx
     do ii=1,nbranch
       lambda(ii)=zero
       if (abs(phfrq_tmp(ii)) > tol10) lambda(ii)=eigval(ii)/(pi*elph_ds%n0(isppol)*phfrq_tmp(ii)**2)
       write (unit_lambda,'(es16.8)',advance='no')phfrq_tmp(ii),lambda(ii)
     end do
     write (unit_lambda,'(es16.8)',advance='no') sum(lambda)
     write (unit_lambda,*)

!    MG NOTE: I wrote a piece of code to output all these quantities using units
!    chosen by the user, maybe in version 5.2?
!    In this version the output of lambda(q,\nu) has been added

!    Output to the main output file, for first point in segment
     if(indxprtqpt(ipoint)==1)then
       write(msg,'(a,a,3es16.6,a,i4,a,a)')ch10,&
&       ' Q point =',qpt(:),'   isppol = ',isppol,ch10,&
&       ' Mode number    Frequency (Ha)  Linewidth (Ha)  Lambda(q,n)'
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       do ii=1,nbranch
         write(msg,'(i8,es20.6,2es16.6)' )ii,phfrq_tmp(ii),eigval(ii),lambda(ii)
         call wrtout(std_out,msg,'COLL')
         call wrtout(ab_out,msg,'COLL')
       end do
     end if

!    find max/min phonon frequency along path chosen
!    presumed to be representative of full BZ to within 10 percent
     elph_ds%omega_min = min(elph_ds%omega_min,1.1_dp*phfrq_tmp(1))
     elph_ds%omega_max = max(elph_ds%omega_max,1.1_dp*phfrq_tmp(nbranch))

     indx = indx+1
   end do !  end ipoint do

!  add blank lines to output files between sppol
   write(msg,'(a)' ) ''
   call wrtout(unit_lwd,msg,'COLL')
   call wrtout(unit_lambda,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end do ! isppol

 ABI_DEALLOCATE(coskr)
 ABI_DEALLOCATE(sinkr)

 close(unit=unit_lwd)
 close(unit=unit_bs)
 close(unit=unit_lambda)

 ABI_DEALLOCATE(finepath)
 ABI_DEALLOCATE(indxprtqpt)

 write(std_out,*) ' elph_linwid : omega_min, omega_max = ',elph_ds%omega_min, elph_ds%omega_max

 DBG_EXIT("COLL")

end subroutine mkph_linwid
!!***
