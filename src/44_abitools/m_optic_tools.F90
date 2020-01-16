!!****m* ABINIT/m_optic_tools
!! NAME
!! m_optic_tools
!!
!! FUNCTION
!!  Helper functions used in the optic code
!!
!! COPYRIGHT
!! Copyright (C) 2002-2019 ABINIT group (SSharma,MVer,VRecoules,TD,YG, NAP)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! COMMENTS
!!
!!  Right now the routine sums over the k-points. In future linear tetrahedron method might be useful.
!!  Reference articles:
!!  1. S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl, Phys. Rev. B {\bf 67} 165332 2003 [[cite:Sharma2003]]
!!  2. J. L. P. Hughes and J. E. Sipe, Phys. Rev. B {\bf 53} 10 751 1996 [[cite:Hughes1996]]
!!  3. S. Sharma and C. Ambrosch-Draxl, Physica Scripta T 109 2004 [[cite:Sharma2004]]
!!  4. J. E. Sipe and Ed. Ghahramani, Phys. Rev. B {\bf 48} 11 705 1993 [[cite:Sipe1993]]
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_optic_tools

 use defs_basis
 use m_errors
 use m_abicore
 use m_linalg_interfaces
 use m_xmpi
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_datatypes,    only : ebands_t
 use m_numeric_tools,   only : wrap2_pmhalf, c2r
 use m_io_tools,        only : open_file

 implicit none

 private

 public :: sym2cart
 public :: getwtk
 public :: pmat2cart
 public :: pmat_renorm
 public :: linopt           ! Compute dielectric function for semiconductors
 public :: nlinopt          ! Second harmonic generation susceptibility for semiconductors
 public :: linelop          ! Linear electro-optic susceptibility for semiconductors
 public :: nonlinopt        ! nonlinear electro-optic susceptibility for semiconductors

CONTAINS  !===========================================================
!!***

!!****f* m_optic_tools/sym2cart
!! NAME
!! sym2cart
!!
!! FUNCTION
!! Routine called by the program optic
!! Convert to symmetry matrice in cartesian coordinates
!!
!! INPUTS
!!	gprimd(3,3)=dimensional primitive translations for reciprocal space
!!	nsym=number of symmetries in group
!!	rprimd(3,3)=dimensional real space primitive translations (bohr)
!!	symrel(3,3,nsym)=symmetry matrices in terms of real space
!!
!! OUTPUT
!!	symcart(3,3)=symmetry matrice in cartesian coordinates (reals)
!!
!! PARENTS
!!      optic
!!
!! CHILDREN
!!      xmpi_max,xmpi_min,xmpi_split_work,xmpi_sum
!!
!! SOURCE


subroutine sym2cart(gprimd,nsym,rprimd,symrel,symcart)

!Arguments -----------------------------------------------
!scalars
 integer,intent(in) :: nsym
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3)
 real(dp),intent(out) :: symcart(3,3,nsym)

!Local variables-------------------------------
!scalars
 integer :: isym
!arrays
 real(dp) :: rsym(3,3),rsymcart(3,3),tmp(3,3)

! *************************************************************************

 do isym=1,nsym
   rsym(:,:) = dble(symrel(:,:,isym))
!  write(std_out,*) 'rsym = ',rsym
   call dgemm('N','N',3,3,3,one,rprimd,3,rsym,  3,zero,tmp,     3)
   call dgemm('N','N',3,3,3,one,tmp,   3,gprimd,3,zero,rsymcart,3)
!  write(std_out,*) 'rsymcart = ',rsymcart
   symcart(:,:,isym) = rsymcart(:,:)
! purify symops in cartesian dp coordinates
   where( abs(symcart(:,:,isym))<tol14)
     symcart(:,:,isym) = zero
   end where
 end do

end subroutine sym2cart
!!***

!----------------------------------------------------------------------

!!****f* m_optic_tools/getwtk
!! NAME
!! getwtk
!!
!! FUNCTION
!! Routine called by the program optic
!! Presumes kpts are the irreducible ones of a good uniform grid
!!
!! INPUTS
!!  kpt(3,nkpt)=reduced coordinates of k points.
!!  nkpt = number of k points
!!  nsym=Number of symmetry operations.
!!  symrel(3,3,nsym)=symmetry operations
!!
!! OUTPUT
!!  wtk(nkpt)=weight assigned to each k point.
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_max,xmpi_min,xmpi_split_work,xmpi_sum
!!
!! SOURCE

subroutine getwtk(kpt,nkpt,nsym,symrel,wtk)

!Arguments -----------------------------------------------
!scalars
 integer,intent(in) :: nkpt,nsym
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: kpt(3,nkpt)
 real(dp),intent(out) :: wtk(nkpt)

!Local variables -----------------------------------------
!scalars
 integer :: ikpt,istar,isym,itim,new,nkpt_tot
 real(dp) :: shift,timsign,tmp
!arrays
 integer :: nstar(nkpt)
 real(dp) :: dkpt(3),kptstar(3,2*nkpt*nsym),rsymrel(3,3,nsym),symkpt(3)
 real(dp) :: tsymkpt(3)

! *************************************************************************

 do isym=1,nsym
   rsymrel(:,:,isym) = dble(symrel(:,:,isym))
 end do

!for each kpt find star and accumulate nkpts
 do ikpt=1,nkpt
   write(std_out,*) ' getwtk : ikpt = ', ikpt
   nstar(ikpt) = 0
   kptstar(:,:) = zero
   do isym=1,nsym

     call dgemv('N',3,3,one,rsymrel(:,:,isym),3,kpt(:,ikpt),1,zero,symkpt,1)

!    is symkpt already in star?
     do itim=0,1
       timsign=one-itim*two
       tsymkpt(:) = timsign*symkpt(:)
       call wrap2_pmhalf(tsymkpt(1),tmp,shift) ;  tsymkpt(1) = tmp
       call wrap2_pmhalf(tsymkpt(2),tmp,shift) ;  tsymkpt(2) = tmp
       call wrap2_pmhalf(tsymkpt(3),tmp,shift) ;  tsymkpt(3) = tmp
       new=1
       do istar=1,nstar(ikpt)
         dkpt(:) = abs(tsymkpt(:)-kptstar(:,istar))
         if ( sum(dkpt) < 1.0d-6) then
           new=0
           exit
         end if
       end do
       if (new==1) then
         nstar(ikpt) = nstar(ikpt)+1
         kptstar(:,nstar(ikpt)) = tsymkpt(:)
       end if
     end do

   end do
!  end do nsym
!  DEBUG
!  write(std_out,*) ' getwtk : nstar = ', nstar(ikpt)
!  write(std_out,*) ' getwtk : star = '
!  write(std_out,*)  kptstar(:,1:nstar(ikpt))
!  ENDDEBUG
 end do
!end do nkpt

 nkpt_tot = sum(nstar)
!write(std_out,*) ' getwtk : nkpt_tot = ', nkpt_tot
 do ikpt=1,nkpt
   wtk(ikpt) = dble(nstar(ikpt))/dble(nkpt_tot)
 end do

end subroutine getwtk
!!***

!----------------------------------------------------------------------

!!****f* m_optic_tools/pmat2cart
!! NAME
!! pmat2cart
!!
!! FUNCTION
!!  turn momentum matrix elements to cartesian axes. To be used in optic calculation of linear
!!  and non-linear RPA dielectric matrices
!!
!! INPUTS
!!  eigen11,eigen12,eigen13 = first order ddk eigen values = d eig_i,k / dk for 3 reduced directions
!!  mband=maximum number of bands
!!  nkpt = number of k-points
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!
!! OUTPUT
!!  pmat(mband,mband,nkpt,3,nsppol) = matrix elements of momentum operator, in cartesian coordinates
!!
!! PARENTS
!!      optic
!!
!! CHILDREN
!!      xmpi_max,xmpi_min,xmpi_split_work,xmpi_sum
!!
!! SOURCE

subroutine pmat2cart(eigen11,eigen12,eigen13,mband,nkpt,nsppol,pmat,rprimd)

!Arguments -----------------------------------------------
!scalars
 integer,intent(in) :: mband,nkpt,nsppol
!arrays
 real(dp),intent(in) :: eigen11(2,mband,mband,nkpt,nsppol)
 real(dp),intent(in) :: eigen12(2,mband,mband,nkpt,nsppol)
 real(dp),intent(in) :: eigen13(2,mband,mband,nkpt,nsppol),rprimd(3,3)
!no_abirules
 complex(dpc),intent(out) :: pmat(mband,mband,nkpt,3,nsppol)

!Local variables -----------------------------------------
!scalars
 integer :: iband1,iband2,ikpt,isppol
!arrays
 real(dp) :: rprim(3,3)

! *************************************************************************

!rescale the rprim
 rprim(:,:) = rprimd(:,:) / two_pi

 do isppol=1,nsppol
   do ikpt=1,nkpt
     do iband1=1,mband
       do iband2=1,mband
         pmat(iband2,iband1,ikpt,:,isppol) =             &
&         rprim(:,1)*cmplx(eigen11(1,iband2,iband1,ikpt,isppol),eigen11(2,iband2,iband1,ikpt,isppol),kind=dp) &
&         +rprim(:,2)*cmplx(eigen12(1,iband2,iband1,ikpt,isppol),eigen12(2,iband2,iband1,ikpt,isppol),kind=dp) &
&         +rprim(:,3)*cmplx(eigen13(1,iband2,iband1,ikpt,isppol),eigen13(2,iband2,iband1,ikpt,isppol),kind=dp)
       end do
     end do
   end do
 end do

end subroutine pmat2cart
!!***

!----------------------------------------------------------------------

!!****f* m_optic_tools/pmat_renorm
!! NAME
!! pmat_renorm
!!
!! FUNCTION
!! Renormalize the momentum matrix elements according to the scissor shift which is imposed
!!
!! INPUTS
!!  mband= number of bands
!!  nkpt = number of k-points
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  efermi = Fermi level
!!  sc = scissor shift for conduction bands
!!  evalv = eigenvalues for ground state
!!
!! OUTPUT
!!  pmat(mband,mband,nkpt,3,nsppol) = momentum matrix elements, renormalized by denominator change with scissor shift
!!
!! PARENTS
!!      optic
!!
!! CHILDREN
!!      xmpi_max,xmpi_min,xmpi_split_work,xmpi_sum
!!
!! SOURCE

subroutine pmat_renorm(efermi, evalv, mband, nkpt, nsppol, pmat, sc)

!Arguments -----------------------------------------------
!scalars
 integer, intent(in) :: nsppol
 integer, intent(in) :: nkpt
 integer, intent(in) :: mband
 real(dp), intent(in) :: efermi
 real(dp), intent(in) :: sc
!arrays
 real(dp), intent(in) :: evalv(mband,nkpt,nsppol)
 complex(dpc), intent(inout) :: pmat(mband,mband,nkpt,3,nsppol)

!Local variables -----------------------------------------
!scalars
 integer :: iband1,iband2,ikpt,isppol
 real(dp) :: corec, e1, e2

! *************************************************************************

 if (abs(sc) < tol8) then
   call wrtout(std_out,' No scissor shift to be applied. Returning to main optic routine.',"COLL")
   return
 end if

 do isppol=1,nsppol
   do ikpt=1,nkpt
     do iband1=1,mband ! valence states
       e1 = evalv(iband1,ikpt,isppol)
       if (e1 > efermi) cycle
       do iband2=1,mband ! conduction states
         e2 = evalv(iband2,ikpt,isppol)
         if (e2 < efermi) cycle
         corec = (e2+sc-e1)/(e2-e1)
         pmat(iband2,iband1,ikpt,:,isppol) = corec * pmat(iband2,iband1,ikpt,:,isppol)
         pmat(iband1,iband2,ikpt,:,isppol) = corec * pmat(iband1,iband2,ikpt,:,isppol)
       end do
     end do
   end do
 end do

end subroutine pmat_renorm
!!***

!----------------------------------------------------------------------

!!****f* m_optic_tools/linopt
!! NAME
!! linopt
!!
!! FUNCTION
!! Compute optical frequency dependent dielectric function for semiconductors
!!
!! INPUTS
!!  icomp=Sequential index associated to computed tensor components (used for netcdf output)
!!  itemp=Temperature index (used for netcdf output)
!!  nspin=number of spins(integer)
!!  omega=crystal volume in au (real)
!!  nkpt=total number of kpoints (integer)
!!  wkpt(nkpt)=weights of kpoints (real)
!!  nsymcrys=number of crystal symmetry operations(integer)
!!  symcrys(3,3,nsymcrys)=symmetry operations in cartisian coordinates(real)
!!  nstval=total number of valence states(integer)
!!  occv(nstval,nkpt,nspin)=occupation number for each band(real)
!!  evalv(nstval,nkpt,nspin)=eigen value for each band in Ha(real)
!!  efermi=Fermi energy in Ha(real)
!!  pmat(nstval,nstval,nkpt,3,nspin)=momentum matrix elements in cartesian coordinates(complex)
!!  v1,v2=desired component of the dielectric function(integer) 1=x,2=y,3=z
!!  nmesh=desired number of energy mesh points(integer)
!!  de=desired step in energy(real); nmesh*de=maximum energy
!!  sc=scissors shift in Ha(real)
!!  brod=broadening in Ha(real)
!!  fnam=root for filename that will contain the output filename will be trim(fnam)//'-linopt.out'
!!  ncid=Netcdf id to save output data.
!!
!! SIDE EFFECTS
!!  Dielectric function for semiconductors, on a desired energy mesh and for a desired
!!  direction of polarisation is written to file.
!!  The output is in a file named trim(fnam)//'-linopt.out' and contains
!!  Im(\epsilon_{v1v2}(\omega), Re(\epsilon_{v1v2}(\omega) and abs(\epsilon_{v1v2}(\omega).
!!  Comment:
!!  Right now the routine sums over the kpoints. In future linear tetrahedron method should be useful.
!!
!! PARENTS
!!      optic
!!
!! CHILDREN
!!      xmpi_max,xmpi_min,xmpi_split_work,xmpi_sum
!!
!! SOURCE

subroutine linopt(icomp,itemp,nspin,omega,nkpt,wkpt,nsymcrys,symcrys,nstval,KSBSt,EPBSt,efermi,pmat, &
  v1,v2,nmesh,de,sc,brod,fnam,ncid,comm)

!Arguments ------------------------------------
integer, intent(in) :: icomp,itemp,nspin,ncid
real(dp), intent(in) :: omega
integer, intent(in) :: nkpt
real(dp), intent(in) :: wkpt(nkpt)
integer, intent(in) :: nsymcrys
real(dp), intent(in) :: symcrys(3,3,nsymcrys)
integer, intent(in) :: nstval
type(ebands_t),intent(in) :: KSBSt,EPBSt
real(dp), intent(in) :: efermi
complex(dpc), intent(in) :: pmat(nstval,nstval,nkpt,3,nspin)
integer, intent(in) :: v1
integer, intent(in) :: v2
integer, intent(in) :: nmesh
real(dp), intent(in) :: de
real(dp), intent(in) :: sc
real(dp), intent(in) :: brod
character(len=*), intent(in) :: fnam
integer, intent(in) :: comm

!Local variables -------------------------
!no_abirules
integer :: isp
integer :: i,j,isym,lx,ly,ik
integer :: ist1,ist2,iw
! Parallelism
integer :: my_rank, nproc
integer,parameter :: master=0
integer :: my_k1, my_k2
#ifdef HAVE_NETCDF
integer :: ncerr
#endif
integer :: ierr
integer :: fout1
logical :: do_linewidth
complex(dpc) :: e1,e2,e12
complex(dpc) :: e1_ep,e2_ep,e12_ep
real(dp) :: deltav1v2
real(dp) :: ha2ev
real(dp) :: tmpabs
real(dp) :: renorm_factor,emin,emax
real(dp) :: ene
complex(dpc) :: b11,b12
complex(dpc) :: ieta,w
character(len=fnlen) :: fnam1
character(len=500) :: msg
! local allocatable arrays
real(dp) :: s(3,3),sym(3,3)
complex(dpc), allocatable :: chi(:)
complex(dpc), allocatable :: eps(:)

! *********************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 if (my_rank == master) then
  !check polarisation
   if (v1.le.0.or.v2.le.0.or.v1.gt.3.or.v2.gt.3) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    Error in linopt:                         '
     write(std_out,*) '    the polarisation directions incorrect    '
     write(std_out,*) '    1=x and 2=y and 3=z                      '
     write(std_out,*) '---------------------------------------------'
     MSG_ERROR("Aborting now")
   end if
  !number of energy mesh points
   if (nmesh.le.0) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    Error in linopt:                         '
     write(std_out,*) '    number of energy mesh points incorrect   '
     write(std_out,*) '    number has to integer greater than 0     '
     write(std_out,*) '    nmesh*de = max energy for calculation    '
     write(std_out,*) '---------------------------------------------'
     MSG_ERROR("Aborting now")
   end if
  !step in energy
   if (de.le.0._dp) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    Error in linopt:                         '
     write(std_out,*) '    energy step is incorrect                 '
     write(std_out,*) '    number has to real greater than 0.0      '
     write(std_out,*) '    nmesh*de = max energy for calculation    '
     write(std_out,*) '---------------------------------------------'
     MSG_ERROR("Aborting now")
   end if
  !broadening
   if (brod.gt.0.009) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    ATTENTION: broadening is quite high      '
     write(std_out,*) '    ideally should be less than 0.005        '
     write(std_out,*) '---------------------------------------------'
   else if (brod.gt.0.015) then
     write(std_out,*) '----------------------------------------'
     write(std_out,*) '    ATTENTION: broadening is too high   '
     write(std_out,*) '    ideally should be less than 0.005   '
     write(std_out,*) '----------------------------------------'
   end if
  !fermi energy
   if(efermi<-1.0d4) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    ATTENTION: Fermi energy seems extremely  '
     write(std_out,*) '    low                                      '
     write(std_out,*) '---------------------------------------------'
   end if
  !scissors operator
   if (sc.lt.0._dp) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    Error in linopt:                         '
     write(std_out,*) '    scissors shift is incorrect              '
     write(std_out,*) '    number has to be greater than 0.0      '
     write(std_out,*) '---------------------------------------------'
     MSG_ERROR("Aborting now")
   end if
!fool proof end
 end if

 ABI_CHECK(KSBSt%mband==nstval, "The number of bands in the BSt should be equal to nstval !")

 do_linewidth = allocated(EPBSt%linewidth)
! TODO: activate this, and remove do_linewidth - always add it in even if 0.
! if (.not. allocated(EPBSt%linewidth)) then
!   ABI_ALLOCATE(EPBSt%linewidth, (1, nstval, my_k2-my_k1+1, nspin))
!   EPBSt%linewidth = zero
! end if

!allocate local arrays
 ABI_ALLOCATE(chi,(nmesh))
 ABI_ALLOCATE(eps,(nmesh))
 ieta=(0._dp,1._dp)*brod
 renorm_factor=1._dp/(omega*dble(nsymcrys))
 ha2ev=13.60569172*2._dp
!output file names
 fnam1=trim(fnam)//'-linopt.out'
!construct symmetrisation tensor
 sym(:,:)=0._dp
 do isym=1,nsymcrys
   s(:,:)=symcrys(:,:,isym)
   do i=1,3
     do j=1,3
       sym(i,j)=sym(i,j)+s(i,v1)*s(j,v2)
     end do
   end do
 end do

!calculate the energy window
 emin=0._dp
 emax=0._dp
 do ik=1,nkpt
   do isp=1,nspin
     do ist1=1,nstval
       emin=min(emin,EPBSt%eig(ist1,ik,isp))
       emax=max(emax,EPBSt%eig(ist1,ik,isp))
     end do
   end do
 end do

 ! Split work
 call xmpi_split_work(nkpt,comm,my_k1,my_k2)

!start calculating linear optical response
 chi(:)=0._dp
! TODO: this loop should be outside the ik one, for speed and cache.
 do isp=1,nspin
   !do ik=1,nkpt
   do ik=my_k1,my_k2
     write(std_out,*) "P-",my_rank,": ",ik,'of',nkpt
     do ist1=1,nstval
       e1=KSBSt%eig(ist1,ik,isp)
       e1_ep=EPBSt%eig(ist1,ik,isp)
! TODO: unless memory is a real issue, should set lifetimes to 0 and do this sum systematically
! instead of putting an if statement in a loop! See above
       if(do_linewidth) then
         e1_ep = e1_ep + EPBSt%linewidth(1,ist1,ik,isp)*(0.0_dp,1.0_dp)
       end if
!      if (e1.lt.efermi) then
!      do ist2=ist1,nstval
       do ist2=1,nstval
         e2=KSBSt%eig(ist2,ik,isp)
         e2_ep=EPBSt%eig(ist2,ik,isp)
         if(do_linewidth) then
           e2_ep = e2_ep - EPBSt%linewidth(1,ist2,ik,isp)*(0.0_dp,1.0_dp)
         end if
!        if (e2.gt.efermi) then
         if (ist1.ne.ist2) then
!          scissors correction of momentum matrix
           if(REAL(e1) > REAL(e2)) then
             e12 = e1-e2+sc
           else
             e12 = e1-e2-sc
           end if
           if(REAL(e1_ep) > REAL(e2_ep)) then
             e12_ep = e1_ep-e2_ep+sc
           else
             e12_ep = e1_ep-e2_ep-sc
           end if
!          e12=e1-e2-sc
           b11=0._dp
!          symmetrization of momentum matrix
           do lx=1,3
             do ly=1,3
               b11=b11+(sym(lx,ly)*pmat(ist1,ist2,ik,lx,isp)* &
               conjg(pmat(ist1,ist2,ik,ly,isp)))
             end do
           end do
           b12=b11*renorm_factor*(1._dp/(e12**2))
!          calculate on the desired energy grid
           do iw=2,nmesh
             w=(iw-1)*de+ieta
             chi(iw)=chi(iw)+(wkpt(ik)*(KSBSt%occ(ist1,ik,isp)-KSBSt%occ(ist2,ik,isp))* &
             (b12/(-e12_ep-w)))
           end do
         end if
       end do ! states
!      end if
     end do
   end do ! spin
 end do ! k-points

 call xmpi_sum(chi,comm,ierr)

 ! calculate epsilon
 eps(1) = zero
 deltav1v2=zero; if (v1 == v2) deltav1v2=one
 do iw=2,nmesh
   ene=(iw-1)*de
   ene=ene*ha2ev
   eps(iw)=deltav1v2+4._dp*pi*chi(iw)
 end do

 if (my_rank == master) then
   !  open the output files
   if (open_file(fnam1,msg,newunit=fout1,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   ! write the output
   write(fout1, '(a,2i3,a)' )' #calculated the component:',v1,v2,'  of dielectric function'
   write(std_out,*) 'calculated the component:',v1,v2,'  of dielectric function'
   write(fout1, '(a,2es16.6)' ) ' #broadening:', real(ieta),aimag(ieta)
   write(std_out,*) ' with broadening:',ieta
   write(fout1, '(a,es16.6)' ) ' #scissors shift:',sc
   write(std_out,*) 'and scissors shift:',sc
   write(fout1, '(a,es16.6,a,es16.6,a)' ) ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(std_out,*) 'energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout1,*)
   write(fout1, '(a)' ) ' # Energy(eV)         Im(eps(w))'
   do iw=2,nmesh
     ene=(iw-1)*de
     ene=ene*ha2ev
     write(fout1, '(2es16.6)' ) ene,aimag(eps(iw))
   end do
   write(fout1,*)
   write(fout1,*)
   write(fout1, '(a)' ) ' # Energy(eV)         Re(eps(w))'
   do iw=2,nmesh
     ene=(iw-1)*de
     ene=ene*ha2ev
     write(fout1, '(2es16.6)' ) ene,dble(eps(iw))
   end do
   write(fout1,*)
   write(fout1,*)
   write(fout1, '(a)' )' # Energy(eV)         abs(eps(w))'
   do iw=2,nmesh
     ene=(iw-1)*de
     ene=ene*ha2ev
     write(fout1, '(2es16.6)' ) ene,abs(eps(iw))
   end do
   write(fout1,*)
   write(fout1,*)
   write(fout1, '(a)' )' # Energy(eV)         Im(refractive index(w)) aka kappa'
   do iw=2,nmesh
     ene=(iw-1)*de
     ene=ene*ha2ev
     write(fout1, '(2es16.6)' ) ene,sqrt(half*(abs(eps(iw)) - dble(eps(iw)) ))
   end do
   write(fout1,*)
   write(fout1,*)
   write(fout1, '(a)' )' # Energy(eV)         Re(refractive index(w)) aka n'
   do iw=2,nmesh
     ene=(iw-1)*de
     ene=ene*ha2ev
     write(fout1, '(2es16.6)' ) ene,sqrt(half*(abs(eps(iw)) + dble(eps(iw)) ))
   end do
   write(fout1,*)
   write(fout1,*)
   write(fout1, '(a)' )' # Energy(eV)         Reflectivity(w) from vacuum, at normal incidence'
   do iw=2,nmesh
     ene=(iw-1)*de
     ene=ene*ha2ev
     write(fout1, '(2es16.6)' ) ene, sqrt(half*(abs(eps(iw)) + dble(eps(iw)) ))
   end do
   write(fout1,*)
   write(fout1,*)
   write(fout1, '(a)' )' # Energy(eV)         absorption coeff (in m-1) = omega Im(eps) / c n(eps)'
   do iw=2,nmesh
     ene=(iw-1)*de
     tmpabs=zero
     if (abs(eps(iw)) + dble(eps(iw)) > zero) then
       tmpabs = aimag(eps(iw))*ene / sqrt(half*( abs(eps(iw)) + dble(eps(iw)) )) / Sp_Lt / Bohr_meter
     end if
     write(fout1, '(2es16.6)' ) ha2ev*ene, tmpabs
   end do

   ! close output file
   close(fout1)

#ifdef HAVE_NETCDF
   if (ncid /= nctk_noid) then
     ncerr = nf90_put_var(ncid, nctk_idname(ncid, "linopt_epsilon"), c2r(eps), start=[1, 1, icomp, itemp])
     NCF_CHECK(ncerr)
   end if
#endif
 end if

 ABI_DEALLOCATE(chi)
 ABI_FREE(eps)

end subroutine linopt
!!***

!----------------------------------------------------------------------

!!****f* m_optic_tools/nlinopt
!! NAME
!! nlinopt
!!
!! FUNCTION
!! Compute optical frequency dependent second harmonic generation susceptibility for semiconductors
!!
!! INPUTS
!!  icomp=Sequential index associated to computed tensor components (used for netcdf output)
!!  itemp=Temperature index (used for netcdf output)
!!  nspin = number of spins(integer)
!!  omega = crystal volume in au (real)
!!  nkpt  = total number of kpoints (integer)
!!  wkpt(nkpt) = weights of kpoints (real)
!!  nsymcrys = number of crystal symmetry operations(integer)
!!  symcrys(3,3,nsymcrys) = symmetry operations in cartisian coordinates(real)
!!  nstval = total number of valence states(integer)
!!  evalv(nstval,nspin,nkpt) = eigen value for each band in Ha(real)
!!  efermi = Fermi energy in Ha(real)
!!  pmat(nstval,nstval,nkpt,3,nspin) = momentum matrix elements in cartesian coordinates(complex)
!!  v1,v2,v3 = desired component of the dielectric function(integer) 1=x,2=y,3=z
!!  nmesh = desired number of energy mesh points(integer)
!!  de = desired step in energy(real); nmesh*de=maximum energy for plotting
!!  sc = scissors shift in Ha(real)
!!  brod = broadening in Ha(real)
!!  tol = tolerance:how close to the singularity exact exact is calculated(real)
!!  fnam=root for filenames that will contain the output  :
!!   fnam1=trim(fnam)//'-ChiTotIm.out'
!!   fnam2=trim(fnam)//'-ChiTotRe.out'
!!   fnam3=trim(fnam)//'-ChiIm.out'
!!   fnam4=trim(fnam)//'-ChiRe.out'
!!   fnam5=trim(fnam)//'-ChiAbs.out'
!!  ncid=Netcdf id to save output data.
!!
!! OUTPUT
!!  Calculates the second harmonic generation susceptibility on a desired energy mesh and
!!  for desired direction of polarisation. The output is in files named
!!  ChiTot.out : Im\chi_{v1v2v3}(2\omega,\omega,-\omega) and Re\chi_{v1v2v3}(2\omega,\omega,-\omega)
!!  ChiIm.out  : contributions to the Im\chi_{v1v2v3}(2\omega,\omega,-\omega) from various terms
!!  ChiRe.out  : contributions to Re\chi_{v1v2v3}(2\omega,\omega,-\omega) from various terms
!!  ChiAbs.out : abs\chi_{v1v2v3}(2\omega,\omega,-\omega). The headers in these files contain
!!  information about the calculation.
!!
!! PARENTS
!!      optic
!!
!! CHILDREN
!!      xmpi_max,xmpi_min,xmpi_split_work,xmpi_sum
!!
!! SOURCE

subroutine nlinopt(icomp,itemp,nspin,omega,nkpt,wkpt,nsymcrys,symcrys,nstval,evalv,efermi, &
  pmat,v1,v2,v3,nmesh,de,sc,brod,tol,fnam,ncid,comm)

!Arguments ------------------------------------
integer, intent(in) :: icomp,itemp,nspin, ncid
real(dp), intent(in) :: omega
integer, intent(in) :: nkpt
real(dp), intent(in) :: wkpt(nkpt)
integer, intent(in) :: nsymcrys
real(dp), intent(in) :: symcrys(3,3,nsymcrys)
integer, intent(in) :: nstval
real(dp), intent(in) :: evalv(nstval,nkpt,nspin)
real(dp), intent(in) :: efermi
complex(dpc), intent(in) :: pmat(nstval,nstval,nkpt,3,nspin)
integer, intent(in) :: v1
integer, intent(in) :: v2
integer, intent(in) :: v3
integer, intent(in) :: nmesh
integer, intent(in) :: comm
real(dp), intent(in) :: de
real(dp), intent(in) :: sc
real(dp), intent(in) :: brod
real(dp), intent(in) :: tol
character(len=*), intent(in) :: fnam

!Local variables -------------------------
integer :: iw
integer :: i,j,k,lx,ly,lz
integer :: isp,isym,ik
integer :: ist1,ist2,istl,istn,istm
integer,parameter :: master=0
integer :: my_rank, nproc
integer :: my_k1, my_k2
integer :: ierr
integer :: fout1,fout2,fout3,fout4,fout5,fout6,fout7
real(dp) :: f1,f2,f3
real(dp) :: ha2ev
real(dp) :: t1,t2,t3,tst
real(dp) :: ene,totre,totabs,totim
real(dp) :: e1,e2,el,en,em
real(dp) :: emin,emax,my_emin,my_emax
real(dp) :: const_esu,const_au,au2esu
real(dp) :: wmn,wnm,wln,wnl,wml,wlm
complex(dpc) :: idel,w,zi
complex(dpc) :: mat2w,mat1w1,mat1w2,mat2w_tra,mat1w3_tra
complex(dpc) :: b111,b121,b131,b112,b122,b132,b113,b123,b133
complex(dpc) :: b241,b242,b243,b221,b222,b223,b211,b212,b213,b231
complex(dpc) :: b311,b312,b313,b331
complex(dpc) :: b24,b21_22,b11,b12_13,b31_32
character(len=fnlen) :: fnam1,fnam2,fnam3,fnam4,fnam5,fnam6,fnam7
character(500) :: msg
! local allocatable arrays
integer :: start4(4),count4(4)
real(dp) :: s(3,3),sym(3,3,3)
complex(dpc), allocatable :: px(:,:,:,:,:)
complex(dpc), allocatable :: py(:,:,:,:,:)
complex(dpc), allocatable :: pz(:,:,:,:,:)
complex(dpc), allocatable :: delta(:,:,:)
complex(dpc), allocatable :: inter2w(:)
complex(dpc), allocatable :: inter1w(:)
complex(dpc), allocatable :: intra2w(:)
complex(dpc), allocatable :: intra1w(:)
complex(dpc), allocatable :: intra1wS(:),chi2tot(:)

! *********************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

!calculate the constant
 zi=(0._dp,1._dp)
 idel=zi*brod
 const_au=-2._dp/(omega*dble(nsymcrys))
 au2esu=5.8300348177d-8
 const_esu=const_au*au2esu
 ha2ev=13.60569172*2._dp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!5.8300348177d-8 : au2esu : bohr*c*10^4/4pi*2*ry2ev
!bohr: 5.2917ifc nlinopt.f907E-11
!c: 2.99792458   velocity of sound
!ry2ev: 13.60569172
!au2esu=(5.29177E-11*2.99792458*1.0E4)/(13.60569172*2)
!this const includes (e^3*hbar^3*hbar^3)/(vol*hbar^5*m_e^3)
!mass comes from converting P_mn to r_mn
!hbar^3 comes from converting all frequencies to energies in denominator
!hbar^3 comes from operator for momentum (hbar/i nabla)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!output file names
 fnam1=trim(fnam)//'-ChiTotIm.out'
 fnam2=trim(fnam)//'-ChiTotRe.out'
 fnam3=trim(fnam)//'-ChiIm.out'
 fnam4=trim(fnam)//'-ChiRe.out'
 fnam5=trim(fnam)//'-ChiAbs.out'
 fnam6=trim(fnam)//'-ChiImDec.out'
 fnam7=trim(fnam)//'-ChiReDec.out'

 if(my_rank == master) then
  !If there exists inversion symmetry exit with a message.
   tst=1.d-09
   do isym=1,nsymcrys
     t1=symcrys(1,1,isym)+1
     t2=symcrys(2,2,isym)+1
     t3=symcrys(3,3,isym)+1
  !  test if diagonal elements are -1
     if (abs(t1).lt.tst.and.abs(t2).lt.tst.and.abs(t3).lt.tst) then
  !    test if off-diagonal elements are zero
       if (abs(symcrys(1,2,isym)).lt.tst.and.abs(symcrys(1,3,isym)).lt.tst &
       .and.abs(symcrys(2,1,isym)).lt.tst.and.abs(symcrys(2,3,isym)).lt.tst.and.  &
       abs(symcrys(3,1,isym)).lt.tst.and.abs(symcrys(3,2,isym)).lt.tst) then
         write(std_out,*) '-----------------------------------------'
         write(std_out,*) '    the crystal has inversion symmetry   '
         write(std_out,*) '    the SHG susceptibility is zero       '
         write(std_out,*) '-----------------------------------------'
         MSG_ERROR("Aborting now")
       end if
     end if
   end do
  !check polarisation
   if (v1.le.0.or.v2.le.0.or.v3.le.0.or.v1.gt.3.or.v2.gt.3.or.v3.gt.3) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    Error in nlinopt:                        '
     write(std_out,*) '    the polarisation directions incorrect    '
     write(std_out,*) '    1=x,  2=y  and 3=z                       '
     write(std_out,*) '---------------------------------------------'
     MSG_ERROR("Aborting now")
   end if
  !number of energy mesh points
   if (nmesh.le.0) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    Error in nlinopt:                        '
     write(std_out,*) '    number of energy mesh points incorrect   '
     write(std_out,*) '    number has to be integer greater than 0  '
     write(std_out,*) '    nmesh*de = max energy for calculation    '
     write(std_out,*) '---------------------------------------------'
     MSG_ERROR("Aborting now")
   end if
  !step in energy
   if (de.le.0._dp) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    Error in nlinopt:                        '
     write(std_out,*) '    energy step is incorrect                 '
     write(std_out,*) '    number has to real greater than 0.0      '
     write(std_out,*) '    nmesh*de = max energy for calculation    '
     write(std_out,*) '---------------------------------------------'
     MSG_ERROR("Aborting now")
   end if
  !broadening
   if (brod.gt.0.009) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    ATTENTION: broadening is quite high      '
     write(std_out,*) '    ideally should be less than 0.005        '
     write(std_out,*) '---------------------------------------------'
   else if (brod.gt.0.015) then
     write(std_out,*) '----------------------------------------'
     write(std_out,*) '    ATTENTION: broadening is too high   '
     write(std_out,*) '    ideally should be less than 0.005   '
     write(std_out,*) '----------------------------------------'
   end if
  !tolerance
   if (tol.gt.0.006) then
     write(std_out,*) '----------------------------------------'
     write(std_out,*) '    ATTENTION: tolerance is too high    '
     write(std_out,*) '    ideally should be less than 0.004   '
     write(std_out,*) '----------------------------------------'
   end if
 end if

 !allocate local arrays
 ABI_ALLOCATE(px,(nstval,nstval,3,3,3))
 ABI_ALLOCATE(py,(nstval,nstval,3,3,3))
 ABI_ALLOCATE(pz,(nstval,nstval,3,3,3))
 ABI_ALLOCATE(inter2w,(nmesh))
 ABI_ALLOCATE(inter1w,(nmesh))
 ABI_ALLOCATE(intra2w,(nmesh))
 ABI_ALLOCATE(intra1w,(nmesh))
 ABI_ALLOCATE(intra1wS,(nmesh))
 ABI_ALLOCATE(delta,(nstval,nstval,3))

!generate the symmetrizing tensor
 sym(:,:,:)=0._dp
 do isym=1,nsymcrys
   s(:,:)=symcrys(:,:,isym)
   do i=1,3
     do j=1,3
       do k=1,3
         sym(i,j,k)=sym(i,j,k)+(s(i,v1)*s(j,v2)*s(k,v3))
       end do
     end do
   end do
 end do
 !DBYG
 ! Disable symmetries for now
 !sym(:,:,:) = 0._dp
 !sym(v1,v2,v3) = nsymcrys
 !ENDDBYG

 ! Split work
 call xmpi_split_work(nkpt,comm,my_k1,my_k2)

!initialise
 inter2w(:)=0._dp
 inter1w(:)=0._dp
 intra2w(:)=0._dp
 intra1w(:)=0._dp
 intra1wS(:)=0._dp
 delta(:,:,:)=0._dp

 my_emin=HUGE(0._dp)
 my_emax=-HUGE(0._dp)
! loop over kpts
 do ik=my_k1,my_k2
   write(std_out,*) "P-",my_rank,": ",ik,'of',nkpt
! loop over spins
   do isp=1,nspin
!  loop over states
     do ist1=1,nstval
       e1=evalv(ist1,ik,isp)
       if (e1.lt.efermi) then   ! ist1 is a valence state
         do ist2=1,nstval
           e2=evalv(ist2,ik,isp)
           if (e2.gt.efermi) then ! ist2 is a conduction state
!            symmetrize the momentum matrix elements
             do lx=1,3
               do ly=1,3
                 do lz=1,3
                   f1=sym(lx,ly,lz)+sym(lx,lz,ly)
                   f2=sym(ly,lx,lz)+sym(ly,lz,lx)
                   f3=sym(lz,lx,ly)+sym(lz,ly,lx)
                   px(ist1,ist2,lx,ly,lz)=f1*pmat(ist1,ist2,ik,lx,isp)
                   py(ist2,ist1,lx,ly,lz)=f2*pmat(ist2,ist1,ik,lx,isp)
                   pz(ist2,ist1,lx,ly,lz)=f3*pmat(ist2,ist1,ik,lx,isp)
                 end do
               end do
             end do
!            end loop over states
           end if
         end do
       end if
     end do
!    calculate the energy window and \Delta_nm
     do ist1=1,nstval
       my_emin=min(my_emin,evalv(ist1,ik,isp))
       my_emax=max(my_emax,evalv(ist1,ik,isp))
       do ist2=1,nstval
         delta(ist1,ist2,1:3)=pmat(ist1,ist1,ik,1:3,isp)-pmat(ist2,ist2,ik,1:3,isp)
       end do
     end do
!    initialise the factors
!    factors are named according to the Ref. article 2.
     b111=0._dp
     b121=0._dp
     b131=0._dp
     b112=0._dp
     b122=0._dp
     b132=0._dp
     b113=0._dp
     b123=0._dp
     b133=0._dp
     b211=0._dp
     b221=0._dp
     b212=0._dp
     b222=0._dp
     b213=0._dp
     b223=0._dp
     b231=0._dp
     b241=0._dp
     b242=0._dp
     b243=0._dp
     b311=0._dp
     b312=0._dp
     b313=0._dp
     b331=0._dp
!    start the calculation
     do istn=1,nstval
       en=evalv(istn,ik,isp)
       if (en.lt.efermi) then    ! istn is a valence state
         do istm=1,nstval
           em=evalv(istm,ik,isp)
           if (em.gt.efermi) then   ! istm is a conduction state
             em = em + sc ! Should add the scissor to conduction energies
             wmn=em-en
             wnm=-wmn
!            calculate the matrix elements for two band intraband term
             mat2w_tra=0._dp
             mat1w3_tra=0._dp
             do lx=1,3
               do ly=1,3
                 do lz=1,3
                   mat2w_tra=mat2w_tra+px(istn,istm,lx,ly,lz)*pmat(istm,istn,ik,lz,isp)    &
                   *delta(istm,istn,ly)
                   mat1w3_tra=mat1w3_tra+px(istn,istm,lx,ly,lz)*pmat(istm,istn,ik,ly,isp)  &
                   *delta(istm,istn,lz)
!                  NOTE:: lx to ly m to n in pmat matrices respectively
!                  Changes are made so that this (b3) term is according to paper
!                  [[cite:Sipe1993]] (Ref. 4) rather than [[cite:Hughes1996]] (Ref 2) in which this term is incorrect
                 end do
               end do
             end do
             b331=mat1w3_tra/wnm
             b11=0._dp
             b12_13=0._dp
             b24=0._dp
             b31_32=0._dp
             b21_22=0._dp

             b231=8._dp*mat2w_tra/wmn
             b331=mat1w3_tra/(wnm)
!            !!!!!!!!!!!!!!!!!!!
!            istl < istn   !
!            !!!!!!!!!!!!!!!!!!!
             do istl=1,istn-1           ! istl is a valence state below istn
               el=evalv(istl,ik,isp)
               wln=el-en                ! do not add sc to the valence bands!
               wml=em-el
               wnl=-wln
               wlm=-wml
!              calculate the matrix elements for three band terms
               mat2w=0._dp
               mat1w1=0._dp
               mat1w2=0._dp
               do lx=1,3
                 do ly=1,3
                   do lz=1,3

                     mat2w=mat2w+(px(istn,istm,lx,ly,lz)*pmat(istm,istl,ik,ly,isp)   &
                     *pmat(istl,istn,ik,lz,isp))

                     mat1w1=mat1w1+(py(istm,istn,lx,ly,lz)*pmat(istl,istm,ik,lz,isp) &
                     *pmat(istn,istl,ik,ly,isp))

                     mat1w2=mat1w2+(pz(istm,istn,lx,ly,lz)*pmat(istl,istm,ik,lz,isp) &
                     *pmat(istn,istl,ik,ly,isp))
                   end do
                 end do
               end do
               b111=mat2w*(1._dp/(wln+wlm))*(1._dp/wlm)
               b121=mat1w1*(1._dp/(wnm+wlm))*(1._dp/wlm)
               b131=mat1w2*(1._dp/wlm)
!
               b221=0._dp
               b211=mat1w1/wml
               b241=-mat2w/wml
!
               b311=mat1w2/wlm
               if (abs(wln).gt.tol) then
                 b111=b111/wln
                 b121=b121/wln
                 b131=b131/wln
                 b221=mat1w2/wln
                 b241=b241+(mat2w/wln)
                 b311=b311+(mat1w1/wln)
               else
                 b111=0._dp
                 b121=0._dp
                 b131=0._dp
                 b221=0._dp
               end if
               t1=wln-wnm
               if (abs(t1).gt.tol) then
                 b131=b131/t1
               else
                 b131=0._dp
               end if
               b11=b11-2._dp*b111
               b12_13=b12_13+b121+b131
               b21_22=b21_22-b211+b221
               b24=b24+2._dp*b241
               b31_32=b31_32+b311
!              end loop over istl
             end do

!            !!!!!!!!!!!!!!!!!!!!!!!!!!!
!            istn < istl < istm    !
!            !!!!!!!!!!!!!!!!!!!!!!!!!!!
             do istl=istn+1,istm-1
               el=evalv(istl,ik,isp)
!              calculate the matrix elements for three band terms
               mat2w=0._dp
               mat1w1=0._dp
               mat1w2=0._dp
               do lx=1,3
                 do ly=1,3
                   do lz=1,3

                     mat2w=mat2w+(px(istn,istm,lx,ly,lz)*pmat(istm,istl,ik,ly,isp)   &
                     *pmat(istl,istn,ik,lz,isp))

                     mat1w1=mat1w1+(py(istm,istn,lx,ly,lz)*pmat(istl,istm,ik,lz,isp) &
                     *pmat(istn,istl,ik,ly,isp))

                     mat1w2=mat1w2+(pz(istm,istn,lx,ly,lz)*pmat(istl,istm,ik,lz,isp) &
                     *pmat(istn,istl,ik,ly,isp))
                   end do
                 end do
               end do
               if (el.lt.efermi) then
                 wln=el-en
                 wnl=-wln
                 wml=em-el
                 wlm=-wml
               else
                 el=el+sc
                 wln=el-en
                 wnl=-wln
                 wml=em-el
                 wlm=-wml
               end if
!
               b112=0._dp
               b122=mat1w1*(1._dp/(wnm+wlm))
               b132=mat1w2*(1._dp/(wnm+wnl))
               b242=0._dp
               b222=0._dp
               b212=0._dp
               b312=0._dp
               if (abs(wnl).gt.tol) then
                 b112=mat2w/wln
                 b122=b122/wnl
                 b132=b132/wnl
                 b242=mat2w/wln
                 b222=mat1w2/wln
                 b312=mat1w1/wln
               else
                 b122=0._dp
                 b132=0._dp
               end if
               if (abs(wlm).gt.tol) then
                 b112=b112/wml
                 b122=b122/wlm
                 b132=b132/wlm
                 b242=b242-(mat2w/wml)
                 b212=mat1w1/wml
                 b312=b312+(mat1w2/wlm)
               else
                 b112=0._dp
                 b122=0._dp
                 b132=0._dp
                 b212=0._dp
               end if
               t1=wlm-wnl
               if (abs(t1).gt.tol) then
                 b112=b112/t1
               else
                 b112=0._dp
               end if
               b11=b11+2._dp*b112
               b12_13=b12_13-b122+b132
               b24=b24+2._dp*b242
               b21_22=b21_22-b212+b222
               b31_32=b31_32+b312
!              end loop over istl
             end do

!            !!!!!!!!!!!!!!!!!!!!!
!            istl > istm    !
!            !!!!!!!!!!!!!!!!!!!!!
             do istl=istm+1,nstval
               el=evalv(istl,ik,isp)+sc
               wln=el-en
               wnl=-wln
               wml=em-el
               wlm=-wml
!              calculate the matrix elements for three band terms
               mat2w=0._dp
               mat1w1=0._dp
               mat1w2=0._dp
               do lx=1,3
                 do ly=1,3
                   do lz=1,3
                     mat2w=mat2w+px(istn,istm,lx,ly,lz)*pmat(istm,istl,ik,ly,isp) &
                     *pmat(istl,istn,ik,lz,isp)

                     mat1w1=mat1w1+(py(istm,istn,lx,ly,lz)*pmat(istl,istm,ik,lz,isp) &
                     *pmat(istn,istl,ik,ly,isp))

                     mat1w2=mat1w2+(pz(istm,istn,lx,ly,lz)*pmat(istl,istm,ik,lz,isp) &
                     *pmat(istn,istl,ik,ly,isp))
                   end do
                 end do
               end do
!
               b113=mat2w*(1._dp/(wnl+wml))*(1._dp/wnl)
               b123=mat1w1*(1._dp/wnl)
               b133=mat1w2*(1._dp/wnl)*(1._dp/(wnl+wnm))
               b243=mat2w/wln
               b223=mat1w2/wln
               b213=0._dp
               b313=-1._dp*mat1w1/wnl
               if (abs(wml).gt.tol) then
                 b113=b113/wml
                 b123=b123/wml
                 b133=b133/wml
                 b243=b243-(mat2w/wml)
                 b213=mat1w1/wml
                 b313=b313+(mat1w2/wlm)
               else
                 b113=0._dp
                 b123=0._dp
                 b133=0._dp
               end if
               t1=wnm-wml
               if (abs(t1).gt.tol) then
                 b123=b123/t1
               else
                 b123=0._dp
               end if
               b11=b11+2._dp*b113
               b12_13=b12_13+b123-b133
               b21_22=b21_22-b213+b223
               b24=b24+2._dp*b243
               b31_32=b31_32+b313
!              end loop over istl
             end do

             b11=b11*zi*(1._dp/wnm)*const_esu
             b12_13=b12_13*zi*(1._dp/wnm)*const_esu
             b24=(b24+b231)*zi*(1._dp/(wnm**3))*const_esu
             b21_22=(b21_22)*zi*(1._dp/(wnm**3))*const_esu
             b31_32=(b31_32-b331)*zi*(1._dp/(wmn**3))*const_esu*0.5_dp
!            calculate over the desired energy mesh and sum over k-points
             do iw=1,nmesh
               w=(iw-1)*de+idel
               inter2w(iw)=inter2w(iw)+(wkpt(ik)*(b11/(wmn-2._dp*w))) ! Inter(2w) from chi
               inter1w(iw)=inter1w(iw)+(wkpt(ik)*(b12_13/(wmn-w))) ! Inter(1w) from chi
               intra2w(iw)=intra2w(iw)+(wkpt(ik)*(b24/(wmn-2._dp*w))) ! Intra(2w) from eta
               intra1w(iw)=intra1w(iw)+(wkpt(ik)*((b21_22)/(wmn-w))) ! Intra(1w) from eta
               intra1wS(iw)=intra1wS(iw)+(wkpt(ik)*((b31_32)/(wmn-w))) ! Intra(1w) from sigma
             end do
           end if
         end do ! istn and istm
       end if
     end do
   end do  ! spins
 end do ! k-points

 call xmpi_sum(inter2w,comm,ierr)
 call xmpi_sum(inter1w,comm,ierr)
 call xmpi_sum(intra2w,comm,ierr)
 call xmpi_sum(intra1w,comm,ierr)
 call xmpi_sum(intra1wS,comm,ierr)
 call xmpi_min(my_emin,emin,comm,ierr)
 call xmpi_max(my_emax,emax,comm,ierr)

 if (my_rank == master) then
   ! write output in SI units and esu (esu to SI(m/v)=(value_esu)*(4xpi)/30000)

   if (ncid /= nctk_noid) then
     start4 = [1, 1, icomp, itemp]
     count4 = [2, nmesh, 1, 1]
     ABI_MALLOC(chi2tot, (nmesh))
     chi2tot = inter2w + inter1w + intra2w + intra1w + intra1wS
#ifdef HAVE_NETCDF
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "shg_inter2w"), c2r(inter2w), start=start4, count=count4))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "shg_inter1w"), c2r(inter1w), start=start4, count=count4))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "shg_intra2w"), c2r(intra2w), start=start4, count=count4))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "shg_intra1w"), c2r(intra1w), start=start4, count=count4))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "shg_intra1wS"), c2r(intra1wS), start=start4, count=count4))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "shg_chi2tot"), c2r(chi2tot), start=start4, count=count4))
#endif
     ABI_FREE(chi2tot)
   end if

   if (open_file(fnam1,msg,newunit=fout1,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam2,msg,newunit=fout2,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam3,msg,newunit=fout3,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam4,msg,newunit=fout4,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam5,msg,newunit=fout5,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam6,msg,newunit=fout6,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam7,msg,newunit=fout7,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
!  write headers
   write(fout1, '(a,3i3)' ) ' #calculated the component:',v1,v2,v3
   write(fout1, '(a,es16.6)' ) ' #tolerance:',tol
   write(fout1, '(a,es16.6,a)' ) ' #broadening:',brod,'Ha'
   write(fout1, '(a,es16.6,a)' ) ' #scissors shift:',sc,'Ha'
   write(fout1, '(a,es16.6,a,es16.6,a)' ) ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout1, '(a)' )' # Energy      Tot-Im Chi(-2w,w,w)  Tot-Im Chi(-2w,w,w)'
   write(fout1, '(a)' )' # eV          *10^-7 esu        *10^-12 m/V SI units '
   write(fout1, '(a)' )' # '

   write(fout2, '(a,3i3)' ) ' #calculated the component:',v1,v2,v3
   write(fout2, '(a,es16.6)') ' #tolerance:',tol
   write(fout2, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout2, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout2, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout2, '(a)')' # Energy      Tot-Re Chi(-2w,w,w)  Tot-Re Chi(-2w,w,w)'
   write(fout2, '(a)')' # eV          *10^-7 esu        *10^-12 m/V SI units '
   write(fout2, '(a)')' # '

   write(fout3, '(a,3i3)') ' #calculated the component:',v1,v2,v3
   write(fout3, '(a,es16.6)') ' #tolerance:',tol
   write(fout3, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout3, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout3, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout3, '(a)')' # Energy(eV) Inter(2w) inter(1w) intra(2w) intra(1w)'
   write(fout3, '(a)')' # in esu'
   write(fout3, '(a)')' # '

   write(fout4, '(a,3i3)') ' #calculated the component:',v1,v2,v3
   write(fout4, '(a,es16.6)') ' #tolerance:',tol
   write(fout4, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout4, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout4, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout4, '(a)')' # Energy(eV) Inter(2w) inter(1w) intra(2w) intra(1w)'
   write(fout4, '(a)')' # in esu'
   write(fout4, '(a)')' # '

   write(fout5, '(a,3i3)') ' #calculated the component:',v1,v2,v3
   write(fout5, '(a,es16.6)') ' #tolerance:',tol
   write(fout5, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout5, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout5, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout5, '(a)')' # Energy(eV)  |TotChi(-2w,w,w)|   |Tot Chi(-2w,w,w)|'
   write(fout5, '(a)')' # eV          *10^-7 esu        *10^-12 m/V SI units '
   write(fout5, '(a)')' # '

   write(fout6, '(a,3i3)') ' #calculated the component:',v1,v2,v3
   write(fout6, '(a,es16.6)') ' #tolerance:',tol
   write(fout6, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout6, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout6, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout6, '(a)')' # Energy(eV) Chi(w) Eta(w) Sigma(w)'
   write(fout6, '(a)')' # in esu'
   write(fout6, '(a)')' # '

   write(fout7, '(a,3i3)') ' #calculated the component:',v1,v2,v3
   write(fout7, '(a,es16.6)') ' #tolerance:',tol
   write(fout7, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout7, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout7, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout7, '(a)')' # Energy(eV) Chi(w) Eta(w) Sigma(w)'
   write(fout7, '(a)')' # in esu'
   write(fout7, '(a)')' # '

   totim=0._dp
   totre=0._dp
   totabs=0._dp
   do iw=2,nmesh
     ene=(iw-1)*de
     ene=ene*ha2ev

     totim=aimag(inter2w(iw)+inter1w(iw)+intra2w(iw)+intra1w(iw)+intra1wS(iw))/1.d-7
     write(fout1,'(f15.6,2es15.6)') ene,totim,totim*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totim=0._dp

     totre=dble(inter2w(iw)+inter1w(iw)+intra2w(iw)+intra1w(iw)+intra1wS(iw))/1.d-7
     write(fout2,'(f15.6,2es15.6)') ene,totre,totre*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totre=0._dp

     write(fout3,'(f15.6,4es15.6)') ene,aimag(inter2w(iw))/1.d-7,      &
     aimag(inter1w(iw))/1.d-7,aimag(intra2w(iw))/1.d-7, aimag(intra1w(iw)+intra1wS(iw))/1.d-7

     write(fout4,'(f15.6,4es15.6)') ene,dble(inter2w(iw))/1.d-7,       &
     dble(inter1w(iw))/1.d-7,dble(intra2w(iw))/1.d-7,dble(intra1w(iw)+intra1wS(iw))/1.d-7

     totabs=abs(inter2w(iw)+inter1w(iw)+intra2w(iw)+intra1w(iw)+intra1wS(iw))/1.d-7
     write(fout5,'(f15.6,2es15.6)') ene,totabs,totabs*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totabs=0._dp

     write(fout6,'(f15.6,4es15.6)') ene,aimag(inter2w(iw)+inter1w(iw))/1.d-7,      &
     aimag(intra2w(iw)+intra1w(iw))/1.d-7,aimag(intra1wS(iw))/1.d-7

     write(fout7,'(f15.6,4es15.6)') ene,dble(inter2w(iw)+inter1w(iw))/1.d-7,       &
     dble(intra2w(iw)+intra1w(iw))/1.d-7,dble(intra1wS(iw))/1.d-7
   end do

   close(fout1)
   close(fout2)
   close(fout3)
   close(fout4)
   close(fout5)
   close(fout6)
   close(fout7)
!  print information
   write(std_out,*) ' '
   write(std_out,*) 'information about calculation just performed:'
   write(std_out,*) ' '
   write(std_out,*) 'calculated the component:',v1,v2,v3 ,'of second order susceptibility'
   write(std_out,*) 'tolerance:',tol
   if (tol.gt.0.008) write(std_out,*) 'ATTENTION: tolerance is too high'
   write(std_out,*) 'broadening:',brod,'Hartree'
   if (brod.gt.0.009) then
     write(std_out,*) ' '
     write(std_out,*) 'ATTENTION: broadening is quite high'
     write(std_out,*) ' '
   else if (brod.gt.0.015) then
     write(std_out,*) ' '
     write(std_out,*) 'ATTENTION: broadening is too high'
     write(std_out,*) ' '
   end if
   write(std_out,*) 'scissors shift:',sc,'Hartree'
   write(std_out,*) 'energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Hartree'
 end if

 ! deallocate local arrays
 ABI_DEALLOCATE(px)
 ABI_DEALLOCATE(py)
 ABI_DEALLOCATE(pz)
 ABI_DEALLOCATE(inter2w)
 ABI_DEALLOCATE(inter1w)
 ABI_DEALLOCATE(intra2w)
 ABI_DEALLOCATE(intra1w)
 ABI_DEALLOCATE(intra1wS)
 ABI_DEALLOCATE(delta)

end subroutine nlinopt
!!***

!----------------------------------------------------------------------

!!****f* m_optic_tools/linelop
!! NAME
!! linelop
!!
!! FUNCTION
!! Compute optical frequency dependent linear electro-optic susceptibility for semiconductors
!!
!! INPUTS
!!  icomp=Sequential index associated to computed tensor components (used for netcdf output)
!!  itemp=Temperature index (used for netcdf output)
!!  nspin = number of spins(integer)
!!  omega = crystal volume in au (real)
!!  nkpt  = total number of kpoints (integer)
!!  wkpt(nkpt) = weights of kpoints (real)
!!  nsymcrys = number of crystal symmetry operations(integer)
!!  symcrys(3,3,nsymcrys) = symmetry operations in cartisian coordinates(real)
!!  nstval = total number of valence states(integer)
!!  evalv(nstval,nspin,nkpt) = eigen value for each band in Ha(real)
!!  occv(nstval,nspin,nkpt) = occupation number
!!  efermi = Fermi energy in Ha(real)
!!  pmat(nstval,nstval,nkpt,3,nspin) = momentum matrix elements in cartesian coordinates(complex)
!!  v1,v2,v3 = desired component of the dielectric function(integer) 1=x,2=y,3=z
!!  nmesh = desired number of energy mesh points(integer)
!!  de = desired step in energy(real); nmesh*de=maximum energy for plotting
!!  sc = scissors shift in Ha(real)
!!  brod = broadening in Ha(real)
!!  tol = tolerance:how close to the singularity exact exact is calculated(real)
!!  fnam=root for filenames that will contain the output  :
!!   fnam1=trim(fnam)//'-ChiTotIm.out'
!!   fnam2=trim(fnam)//'-ChiTotRe.out'
!!   fnam3=trim(fnam)//'-ChiIm.out'
!!   fnam4=trim(fnam)//'-ChiRe.out'
!!   fnam5=trim(fnam)//'-ChiAbs.out'
!!  ncid=Netcdf id to save output data.
!!
!! OUTPUT
!!  Calculates the second harmonic generation susceptibility on a desired energy mesh and
!!  for desired direction of polarisation. The output is in files named
!!  ChiEOTot.out : Im\chi_{v1v2v3}(\omega,\omega,0) and Re\chi_{v1v2v3}(\omega,\omega,0)
!!  ChiEOIm.out  : contributions to the Im\chi_{v1v2v3}(\omega,\omega,0) from various terms
!!  ChiEORe.out  : contributions to Re\chi_{v1v2v3}(\omega,\omega,-0) from various terms
!!  ChiEOAbs.out : abs\chi_{v1v2v3}(\omega,\omega,0). The headers in these files contain
!!  information about the calculation.
!!
!!  NOTES:
!!    - The routine has been written using notations of Ref. 2
!!    - This routine does not symmetrize the tensor (up to now)
!!    - Sum over all the states and use occupation factors instead of looping only on resonant contributions
!!
!! PARENTS
!!      optic
!!
!! CHILDREN
!!      xmpi_max,xmpi_min,xmpi_split_work,xmpi_sum
!!
!! SOURCE

subroutine linelop(icomp,itemp,nspin,omega,nkpt,wkpt,nsymcrys,symcrys,nstval,evalv,occv,efermi, &
  pmat,v1,v2,v3,nmesh,de,sc,brod,tol,fnam,do_antiresonant,ncid,comm)

!Arguments ------------------------------------
integer, intent(in) :: icomp,itemp,nspin, ncid
real(dp), intent(in) :: omega
integer, intent(in) :: nkpt
real(dp), intent(in) :: wkpt(nkpt)
integer, intent(in) :: nsymcrys
real(dp), intent(in) :: symcrys(3,3,nsymcrys)
integer, intent(in) :: nstval
real(dp), intent(in) :: evalv(nstval,nkpt,nspin)
real(dp), intent(in) :: occv(nstval,nkpt,nspin)
real(dp), intent(in) :: efermi
complex(dpc), intent(in) :: pmat(nstval,nstval,nkpt,3,nspin)
integer, intent(in) :: v1
integer, intent(in) :: v2
integer, intent(in) :: v3
integer, intent(in) :: nmesh
integer, intent(in) :: comm
real(dp), intent(in) :: de
real(dp), intent(in) :: sc
real(dp), intent(in) :: brod
real(dp), intent(in) :: tol
character(len=*), intent(in) :: fnam
logical, intent(in) :: do_antiresonant

!Local variables -------------------------
!no_abirules
! present calculation related (user specific)
integer :: iw
integer :: i,j,k,lx,ly,lz
integer :: isp,isym,ik
integer :: ist1,istl,istn,istm
real(dp) :: ha2ev
real(dp) :: t1,t2,t3,tst
real(dp) :: ene,totre,totabs,totim
real(dp) :: el,en,em
real(dp) :: emin,emax,my_emin,my_emax
real(dp) :: const_esu,const_au,au2esu
real(dp) :: wmn,wnm,wln,wnl,wml,wlm
complex(dpc) :: idel,w,zi
character(len=fnlen) :: fnam1,fnam2,fnam3,fnam4,fnam5
! local allocatable arrays
real(dp), allocatable :: s(:,:)
real(dp), allocatable :: sym(:,:,:)
integer :: start4(4),count4(4)
! DBYG
 integer :: istp
 real(dp) :: ep, wmp, wpn
 real(dp), allocatable :: enk(:) ! (n) = \omega_n(k), with scissor included !
 real(dp) :: fn, fm, fl, fnm, fnl, fml, fln, fmn
 complex(dpc), allocatable :: delta(:,:,:) ! (m,n,a) = \Delta_{mn}^{a}
 complex(dpc), allocatable :: rmna(:,:,:) ! (m,n,a) = r_{mn}^{a}
 complex(dpc), allocatable :: rmnbc(:,:,:,:) ! (m,n,b,c) = r^b_{mn;c}(k)
 complex(dpc), allocatable :: roverw(:,:,:,:) ! (m,n,b,c) = [r^b_{mn}(k)/w_{mn(k)];c
 complex(dpc), allocatable :: chi(:) ! \chi_{II}^{abc}(-\omega,\omega,0)
 complex(dpc), allocatable :: eta(:) ! \eta_{II}^{abc}(-\omega,\omega,0)
 complex(dpc), allocatable :: sigma(:) ! \frac{i}{\omega} \sigma_{II}^{abc}(-\omega,\omega,0)
 complex(dpc), allocatable :: chi2tot(:)
 complex(dpc) :: num1, num2, den1, den2, term1, term2
 complex(dpc) :: chi1, chi1_1, chi1_2, chi2_1b, chi2_2b
 complex(dpc), allocatable :: chi2(:) ! Second term that depends on the frequency ! (omega)
 complex(dpc) :: eta1, eta2, eta2_1, eta2_2
 complex(dpc) :: sigma1, sigma1_1, sigma1_2, sigma2
 !Parallelism
 integer :: my_rank, nproc, master=0
 integer :: ierr
 integer :: my_k1, my_k2
 character(500) :: msg
 integer :: fout1,fout2,fout3,fout4,fout5

! *********************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

!calculate the constant
 zi=(0._dp,1._dp)
 idel=zi*brod
! Disable symmetries for now
 const_au=-2._dp/(omega*dble(nsymcrys))
 au2esu=5.8300348177d-8
 const_esu=const_au*au2esu
 ha2ev=13.60569172*2._dp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!5.8300348177d-8 : au2esu : bohr*c*10^4/4pi*2*ry2ev
!bohr: 5.2917ifc nlinopt.f907E-11
!c: 2.99792458   velocity of sound
!ry2ev: 13.60569172
!au2esu=(5.29177E-11*2.99792458*1.0E4)/(13.60569172*2)
!this const includes (e^3*hbar^3*hbar^3)/(vol*hbar^5*m_e^3)
!mass comes from converting P_mn to r_mn
!hbar^3 comes from converting all frequencies to energies in denominator
!hbar^3 comes from operator for momentum (hbar/i nabla)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!output file names
 fnam1=trim(fnam)//'-ChiEOTotIm.out'
 fnam2=trim(fnam)//'-ChiEOTotRe.out'
 fnam3=trim(fnam)//'-ChiEOIm.out'
 fnam4=trim(fnam)//'-ChiEORe.out'
 fnam5=trim(fnam)//'-ChiEOAbs.out'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fool proof:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!If there exists inversion symmetry exit with a mesg.
 tst=1.d-09
 do isym=1,nsymcrys
   t1=symcrys(1,1,isym)+1
   t2=symcrys(2,2,isym)+1
   t3=symcrys(3,3,isym)+1
!  test if diagonal elements are -1
   if (abs(t1).lt.tst.and.abs(t2).lt.tst.and.abs(t3).lt.tst) then
!    test if off-diagonal elements are zero
     if (abs(symcrys(1,2,isym)).lt.tst.and.abs(symcrys(1,3,isym)).lt.tst &
     .and.abs(symcrys(2,1,isym)).lt.tst.and.abs(symcrys(2,3,isym)).lt.tst.and.  &
     abs(symcrys(3,1,isym)).lt.tst.and.abs(symcrys(3,2,isym)).lt.tst) then
       write(std_out,*) '-----------------------------------------'
       write(std_out,*) '    the crystal has inversion symmetry   '
       write(std_out,*) '    the LEO susceptibility is zero       '
       write(std_out,*) '-----------------------------------------'
       MSG_ERROR("Aborting now")
     end if
   end if
 end do
!check polarisation
 if (v1.le.0.or.v2.le.0.or.v3.le.0.or.v1.gt.3.or.v2.gt.3.or.v3.gt.3) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in linelop:                        '
   write(std_out,*) '    the polarisation directions incorrect    '
   write(std_out,*) '    1=x,  2=y  and 3=z                       '
   write(std_out,*) '---------------------------------------------'
   MSG_ERROR("Aborting now")
 end if
!number of energy mesh points
 if (nmesh.le.0) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in linelop:                        '
   write(std_out,*) '    number of energy mesh points incorrect   '
   write(std_out,*) '    number has to be integer greater than 0  '
   write(std_out,*) '    nmesh*de = max energy for calculation    '
   write(std_out,*) '---------------------------------------------'
   MSG_ERROR("Aborting now")
 end if
!step in energy
 if (de.le.0._dp) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in linelop:                        '
   write(std_out,*) '    energy step is incorrect                 '
   write(std_out,*) '    number has to real greater than 0.0      '
   write(std_out,*) '    nmesh*de = max energy for calculation    '
   write(std_out,*) '---------------------------------------------'
   MSG_ERROR("Aborting now")
 end if
!broadening
 if (brod.gt.0.009) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    ATTENTION: broadening is quite high      '
   write(std_out,*) '    ideally should be less than 0.005        '
   write(std_out,*) '---------------------------------------------'
 else if (brod.gt.0.015) then
   write(std_out,*) '----------------------------------------'
   write(std_out,*) '    ATTENTION: broadening is too high   '
   write(std_out,*) '    ideally should be less than 0.005   '
   write(std_out,*) '----------------------------------------'
 end if
!tolerance
 if (tol.gt.0.006) then
   write(std_out,*) '----------------------------------------'
   write(std_out,*) '    ATTENTION: tolerance is too high    '
   write(std_out,*) '    ideally should be less than 0.004   '
   write(std_out,*) '----------------------------------------'
 end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!fool proof ends
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!allocate local arrays
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 ABI_MALLOC(enk,(nstval))
 ABI_MALLOC(delta,(nstval,nstval,3))
 ABI_MALLOC(rmnbc,(nstval,nstval,3,3))
 ABI_MALLOC(roverw,(nstval,nstval,3,3))
 ABI_MALLOC(rmna,(nstval,nstval,3))
 ABI_MALLOC(chi,(nmesh))
 ABI_MALLOC(eta,(nmesh))
 ABI_MALLOC(sigma,(nmesh))
 ABI_MALLOC(chi2,(nmesh))
 ABI_MALLOC(sym,(3,3,3))
 ABI_MALLOC(s,(3,3))

!generate the symmetrizing tensor
 sym(:,:,:)=0._dp
 do isym=1,nsymcrys
   s(:,:)=symcrys(:,:,isym)
   do i=1,3
     do j=1,3
       do k=1,3
         sym(i,j,k)=sym(i,j,k)+(s(i,v1)*s(j,v2)*s(k,v3))
       end do
     end do
   end do
 end do

!initialise
 delta(:,:,:)=0._dp
 rmnbc(:,:,:,:)=0._dp
 chi(:)=0._dp
 chi2(:) = 0._dp
 eta(:)=0._dp
 sigma(:)=0._dp
 my_emin=HUGE(0._dp)
 my_emax=-HUGE(0._dp)

 ! Split work
 call xmpi_split_work(nkpt,comm,my_k1,my_k2)

! loop over kpts
 do ik=my_k1,my_k2
   write(std_out,*) "P-",my_rank,": ",ik,'of',nkpt
!  loop over spins
   do isp=1,nspin
!    Calculate the scissor corrected energies and the energy window
     do ist1=1,nstval
       en = evalv(ist1,ik,isp)
       my_emin=min(my_emin,en)
       my_emax=max(my_emax,en)
       if(en > efermi) then
         en = en + sc
       end if
       enk(ist1) = en
     end do

!    calculate \Delta_nm and r_mn^a
     do istn=1,nstval
       en = enk(istn)
       do istm=1,nstval
         em = enk(istm)
         wmn = em - en
         delta(istn,istm,1:3)=pmat(istn,istn,ik,1:3,isp)-pmat(istm,istm,ik,1:3,isp)
         if(abs(wmn) < tol) then
           rmna(istm,istn,1:3) = 0._dp
         else
           rmna(istm,istn,1:3)=-zi*pmat(istm,istn,ik,1:3,isp)/wmn
         end if
       end do
     end do
!    calculate \r^b_mn;c
     do istm=1,nstval
       em = enk(istm)
       do istn=1,nstval
         en = enk(istn)
         wmn = em - en
         if(abs(wmn) > tol) then
           do ly = 1,3
             do lz = 1,3
               num1 = (rmna(istm,istn,ly)*delta(istm,istn,lz))+(rmna(istm,istn,lz)*delta(istm,istn,ly))
               den1 = wmn
               term1 = num1/den1
               term2 = 0._dp
               do istp=1,nstval
                 ep = enk(istp)
                 wmp = em - ep
                 wpn = ep - en
                 num2 = (wmp*rmna(istm,istp,ly)*rmna(istp,istn,lz))-(wpn*rmna(istm,istp,lz)*rmna(istp,istn,ly))
                 den2 = wmn
                 term2 = term2 + (num2/den2)
               end do
               rmnbc(istm,istn,ly,lz) = -term1-(zi*term2)
               roverw(istm,istn,ly,lz) = (rmnbc(istm,istn,ly,lz)/wmn) - (rmna(istm,istn,ly)/(wmn**2))*delta(istm,istn,lz)
             end do
           end do
         end if
       end do
     end do

!    initialise the factors
!    start the calculation
     do istn=1,nstval
       en=enk(istn)
       if (do_antiresonant .and. en .ge. efermi) then
         cycle
       end if
       fn=occv(istn,ik,isp)
       do istm=1,nstval
         em=enk(istm)
         if (do_antiresonant .and. em .le. efermi) then
           cycle
         end if
         wmn=em-en
         wnm=-wmn
         fm = occv(istm,ik,isp)
         fnm = fn - fm
         fmn = fm - fn
         eta1 = 0._dp
         eta2_1 = 0._dp
         eta2_2 = 0._dp
         sigma1_1 = 0._dp
         sigma1_2 = 0._dp
         sigma2 = 0._dp
         if(abs(wmn) > tol) then
           do lx = 1,3
             do ly = 1,3
               do lz = 1,3
                 eta1 = eta1 + sym(lx,ly,lz)*(fnm*rmna(istn,istm,lx)*(roverw(istm,istn,lz,ly)))
                 eta2_1 = eta2_1 + sym(lx,ly,lz)*(fnm*(rmna(istn,istm,lx)*rmnbc(istm,istn,ly,lz)))
                 eta2_2 = eta2_2 + sym(lx,ly,lz)*(fnm*(rmnbc(istn,istm,lx,lz)*rmna(istm,istn,ly)))
                 sigma1_1 = sigma1_1 + sym(lx,ly,lz)*(fnm*delta(istn,istm,lx)*rmna(istn,istm,ly)*rmna(istm,istn,lz))/(wmn**2)
                 sigma1_2 = sigma1_2 + sym(lx,ly,lz)*(fnm*delta(istn,istm,lx)*rmna(istn,istm,lz)*rmna(istm,istn,ly))/(wmn**2)
                 sigma2 = sigma2 + sym(lx,ly,lz)*(fnm*rmnbc(istn,istm,lz,lx)*rmna(istm,istn,ly))/wmn
               end do
             end do
           end do
         end if
         chi1_1 = 0._dp
         chi1_2 = 0._dp
         chi2_1b = 0._dp
         chi2_2b = 0._dp
         chi2(:) = 0._dp
!        Three band terms
         do istl=1,nstval
           el=enk(istl)
           fl = occv(istl,ik,isp)
           wlm = el-em
           wln = el-en
           wnl = en-el
           wml = em-el
           fnl = fn-fl
           fln = fl-fn
           fml = fm-fl
           do lx = 1,3
             do ly = 1,3
               do lz = 1,3
                 if(abs(wlm) > tol) then
                   chi1_1 = chi1_1 + sym(lx,ly,lz)*(fnm*rmna(istn,istm,lx)*rmna(istm,istl,lz)*rmna(istl,istn,ly))/(wlm)
                   chi2_1b = chi2_1b + sym(lx,ly,lz)*(fnm*rmna(istn,istl,lx)*rmna(istl,istm,lz)*rmna(istm,istn,ly))/(wlm)
                 end if
                 if(abs(wln) > tol) then
                   chi1_2 = chi1_2 + sym(lx,ly,lz)*(fnm*rmna(istn,istm,lx)*rmna(istm,istl,ly)*rmna(istl,istn,lz))/(wln)
                   chi2_2b = chi2_2b + sym(lx,ly,lz)*(fmn*rmna(istl,istm,lx)*rmna(istm,istn,ly)*rmna(istn,istl,lz))/(wnl)
                 end if
               end do
             end do
           end do
         end do

         sigma1 = 0.5_dp*(sigma1_1-sigma1_2)
         eta2 = 0.5_dp*(eta2_1-eta2_2)
         chi1 = chi1_1 + chi1_2
!
!        calculate over the desired energy mesh and sum over k-points
         do iw=1,nmesh
           w=(iw-1)*de+idel
           ! Better way to compute it
           chi(iw) = chi(iw) + 0.5_dp*wkpt(ik)*((chi1/(wmn-w)) + ((chi2_1b+chi2_2b)/(wmn-w)))*const_esu
           eta(iw) = eta(iw) + 0.5_dp*zi*wkpt(ik)*((eta1/(wmn-w)) + (eta2/((wmn-w)**2)))*const_esu
           sigma(iw) = sigma(iw) + 0.5_dp*zi*wkpt(ik)*((sigma1/(wmn-w))- (sigma2/(wmn-w)))*const_esu
         end do
       end do ! istn and istm
     end do
   end do ! spins
 end do ! k-points

 call xmpi_sum(chi,comm,ierr)
 call xmpi_sum(eta,comm,ierr)
 call xmpi_sum(sigma,comm,ierr)
 call xmpi_min(my_emin,emin,comm,ierr)
 call xmpi_max(my_emax,emax,comm,ierr)

 if (my_rank == master) then

   if (ncid /= nctk_noid) then
     start4 = [1, 1, icomp, itemp]
     count4 = [2, nmesh, 1, 1]
     ABI_MALLOC(chi2tot, (nmesh))
     chi2tot = chi + eta + sigma
#ifdef HAVE_NETCDF
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "leo_chi"), c2r(chi), start=start4, count=count4))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "leo_eta"), c2r(eta), start=start4, count=count4))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "leo_sigma"), c2r(sigma), start=start4, count=count4))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "leo_chi2tot"), c2r(chi2tot), start=start4, count=count4))
#endif
     ABI_FREE(chi2tot)
   end if

  ! write output in SI units and esu (esu to SI(m/v)=(value_esu)*(4xpi)/30000)
   if (open_file(fnam1,msg,newunit=fout1,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam2,msg,newunit=fout2,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam3,msg,newunit=fout3,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam4,msg,newunit=fout4,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam5,msg,newunit=fout5,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   ! write headers
   write(fout1, '(a,3i3)' ) ' #calculated the component:',v1,v2,v3
   write(fout1, '(a,es16.6)' ) ' #tolerance:',tol
   write(fout1, '(a,es16.6,a)' ) ' #broadening:',brod,'Ha'
   write(fout1, '(a,es16.6,a)' ) ' #scissors shift:',sc,'Ha'
   write(fout1, '(a,es16.6,a,es16.6,a)' ) ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout1, '(a)' )' # Energy      Tot-Im Chi(-w,w,0)  Tot-Im Chi(-w,w,0)'
   write(fout1, '(a)' )' # eV          *10^-7 esu        *10^-12 m/V SI units '
   write(fout1, '(a)' )' # '

   write(fout2, '(a,3i3)' ) ' #calculated the component:',v1,v2,v3
   write(fout2, '(a,es16.6)') ' #tolerance:',tol
   write(fout2, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout2, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout2, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout2, '(a)')' # Energy      Tot-Re Chi(-w,w,0)  Tot-Re Chi(-w,w,0)'
   write(fout2, '(a)')' # eV          *10^-7 esu        *10^-12 m/V SI units '
   write(fout2, '(a)')' # '

   write(fout3, '(a,3i3)') ' #calculated the component:',v1,v2,v3
   write(fout3, '(a,es16.6)') ' #tolerance:',tol
   write(fout3, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout3, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout3, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout3, '(a)')' # Energy(eV) Chi(w) Eta(w) Sigma(w)'
   write(fout3, '(a)')' # in esu'
   write(fout3, '(a)')' # '

   write(fout4, '(a,3i3)') ' #calculated the component:',v1,v2,v3
   write(fout4, '(a,es16.6)') ' #tolerance:',tol
   write(fout4, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout4, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout4, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout4, '(a)')' # Energy(eV) Chi(w) Eta(w) Sigma(w)'
   write(fout4, '(a)')' # in esu'
   write(fout4, '(a)')' # '

   write(fout5, '(a,3i3)') ' #calculated the component:',v1,v2,v3
   write(fout5, '(a,es16.6)') ' #tolerance:',tol
   write(fout5, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout5, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout5, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout5, '(a)')' # Energy(eV)  |TotChi(-w,w,0)|   |Tot Chi(-w,w,0)|'
   write(fout5, '(a)')' # eV          *10^-7 esu        *10^-12 m/V SI units '
   write(fout5, '(a)')' # '

   totim=0._dp
   totre=0._dp
   totabs=0._dp
   do iw=2,nmesh
     ene=(iw-1)*de
     ene=ene*ha2ev
     totim=aimag(chi(iw)+eta(iw)+sigma(iw))/1.d-7
     write(fout1,'(f15.6,2es15.6)') ene,totim,totim*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totim=0._dp
     totre=dble(chi(iw)+eta(iw)+sigma(iw))/1.d-7
     write(fout2,'(f15.6,2es15.6)') ene,totre,totre*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totre=0._dp
     write(fout3,'(f15.6,3es15.6)') ene,aimag(chi(iw))/1.d-7,      &
     aimag(eta(iw))/1.d-7,aimag(sigma(iw))/1.d-7
     write(fout4,'(f15.6,3es15.6)') ene,dble(chi(iw))/1.d-7,       &
     dble(eta(iw))/1.d-7,dble(sigma(iw))/1.d-7
     totabs=abs(chi(iw)+eta(iw)+sigma(iw))/1.d-7
     write(fout5,'(f15.6,2es15.6)') ene,totabs,totabs*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totabs=0._dp
   end do

   close(fout1)
   close(fout2)
   close(fout3)
   close(fout4)
   close(fout5)
   ! print information
   write(std_out,*) ' '
   write(std_out,*) 'information about calculation just performed:'
   write(std_out,*) ' '
   write(std_out,*) 'calculated the component:',v1,v2,v3 ,'of LEO susceptibility'
   write(std_out,*) 'tolerance:',tol
   if (tol.gt.0.008) write(std_out,*) 'ATTENTION: tolerance is too high'
   write(std_out,*) 'broadening:',brod,'Hartree'
   if (brod.gt.0.009) then
     write(std_out,*) ' '
     write(std_out,*) 'ATTENTION: broadening is quite high'
     write(std_out,*) ' '
   else if (brod.gt.0.015) then
     write(std_out,*) ' '
     write(std_out,*) 'ATTENTION: broadening is too high'
     write(std_out,*) ' '
   end if
   write(std_out,*) 'scissors shift:',sc,'Hartree'
   write(std_out,*) 'energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Hartree'

 end if

 ! deallocate local arrays
 ABI_FREE(enk)
 ABI_FREE(delta)
 ABI_FREE(rmnbc)
 ABI_FREE(roverw)
 ABI_FREE(rmna)
 ABI_FREE(chi)
 ABI_FREE(chi2)
 ABI_FREE(eta)
 ABI_FREE(sigma)
 ABI_FREE(s)
 ABI_FREE(sym)

end subroutine linelop
!!***

!----------------------------------------------------------------------

!!****f* m_optic_tools/nonlinopt
!! NAME
!! nonlinopt
!!
!! FUNCTION
!! Compute the frequency dependent nonlinear electro-optic susceptibility for semiconductors
!!
!! INPUTS
!!  icomp=Sequential index associated to computed tensor components (used for netcdf output)
!!  itemp=Temperature index (used for netcdf output)
!!  nspin = number of spins(integer)
!!  omega = crystal volume in au (real)
!!  nkpt  = total number of kpoints (integer)
!!  wkpt(nkpt) = weights of kpoints (real)
!!  nsymcrys = number of crystal symmetry operations(integer)
!!  symcrys(3,3,nsymcrys) = symmetry operations in cartisian coordinates(real)
!!  nstval = total number of valence states(integer)
!!  evalv(nstval,nspin,nkpt) = eigen value for each band in Ha(real)
!!  occv(nstval,nspin,nkpt) = occupation number
!!  efermi = Fermi energy in Ha(real)
!!  pmat(nstval,nstval,nkpt,3,nspin) = momentum matrix elements in cartesian coordinates(complex)
!!  v1,v2,v3 = desired component of the dielectric function(integer) 1=x,2=y,3=z
!!  nmesh = desired number of energy mesh points(integer)
!!  de = desired step in energy(real); nmesh*de=maximum energy for plotting
!!  sc = scissors shift in Ha(real)
!!  brod = broadening in Ha(real)
!!  tol = tolerance:how close to the singularity exact exact is calculated(real)
!!  fnam=root for filenames that will contain the output  :
!!   fnam1=trim(fnam)//'-ChiTotIm.out'
!!   fnam2=trim(fnam)//'-ChiTotRe.out'
!!   fnam3=trim(fnam)//'-ChiIm.out'
!!   fnam4=trim(fnam)//'-ChiRe.out'
!!   fnam5=trim(fnam)//'-ChiAbs.out'
!!
!! OUTPUT
!!  Calculates the nonlinear electro-optical susceptibility on a desired energy mesh and
!!  for desired direction of polarisation. The output is in files named
!!  ChiEOTot.out : Im\chi_{v1v2v3}(\omega,\omega,0) and Re\chi_{v1v2v3}(\omega,\omega,0)
!!  ChiEOIm.out  : contributions to the Im\chi_{v1v2v3}(\omega,\omega,0) from various terms
!!  ChiEORe.out  : contributions to Re\chi_{v1v2v3}(\omega,\omega,-0) from various terms
!!  ChiEOAbs.out : abs\chi_{v1v2v3}(\omega,\omega,0). The headers in these files contain
!!  information about the calculation.
!!  ncid=Netcdf id to save output data.
!!
!! COMMENTS
!!    - The routine has been written using notations of Ref. 2
!!    - This routine does not symmetrize the tensor (up to now)
!!    - Sum over all the states and use occupation factors instead of looping only on resonant contributions
!!
!! PARENTS
!!      optic
!!
!! CHILDREN
!!      xmpi_max,xmpi_min,xmpi_split_work,xmpi_sum
!!
!! SOURCE

subroutine nonlinopt(icomp,itemp,nspin,omega,nkpt,wkpt,nsymcrys,symcrys,nstval,evalv,occv,efermi, &
  pmat,v1,v2,v3,nmesh,de,sc,brod,tol,fnam,do_antiresonant,ncid,comm)

!Arguments ------------------------------------
integer, intent(in) :: icomp,itemp,nspin, ncid
real(dp), intent(in) :: omega
integer, intent(in) :: nkpt
real(dp), intent(in) :: wkpt(nkpt)
integer, intent(in) :: nsymcrys
real(dp), intent(in) :: symcrys(3,3,nsymcrys)
integer, intent(in) :: nstval
real(dp), intent(in) :: evalv(nstval,nkpt,nspin)
real(dp), intent(in) :: occv(nstval,nkpt,nspin)
real(dp), intent(in) :: efermi
complex(dpc), intent(in) :: pmat(nstval,nstval,nkpt,3,nspin)
integer, intent(in) :: v1
integer, intent(in) :: v2
integer, intent(in) :: v3
integer, intent(in) :: nmesh
integer, intent(in) :: comm
real(dp), intent(in) :: de
real(dp), intent(in) :: sc
real(dp), intent(in) :: brod
real(dp), intent(in) :: tol
character(len=*), intent(in) :: fnam
logical, intent(in) :: do_antiresonant

!Local variables -------------------------
integer :: iw
integer :: i,j,k,lx,ly,lz
integer :: isp,isym,ik
integer :: ist1,istl,istn,istm
real(dp) :: ha2ev
real(dp) :: t1,t2,t3,tst
real(dp) :: ene,totre,totabs,totim
real(dp) :: el,en,em
real(dp) :: emin,emax, my_emin,my_emax
real(dp) :: const_esu,const_au,au2esu
real(dp) :: wmn,wnm,wln,wnl,wml,wlm
complex(dpc) :: idel,w,zi
character(len=fnlen) :: fnam1,fnam2,fnam3,fnam4,fnam5,fnam6,fnam7
! local allocatable arrays
 integer :: start4(4),count4(4)
 real(dp) :: s(3,3),sym(3,3,3)
 integer :: istp
 real(dp) :: ep, wmp, wpn
 real(dp), allocatable :: enk(:) ! (n) = \omega_n(k), with scissor included !
 real(dp) :: fn, fm, fl, fnm, fnl, fml, fln, flm
 complex(dpc), allocatable :: delta(:,:,:) ! (m,n,a) = \Delta_{mn}^{a}
 complex(dpc), allocatable :: rmna(:,:,:) ! (m,n,a) = r_{mn}^{a}
 complex(dpc), allocatable :: rmnbc(:,:,:,:) ! (m,n,b,c) = r^b_{mn;c}(k)
 complex(dpc), allocatable :: roverw(:,:,:,:) ! (m,n,b,c) = [r^b_{mn}(k)/w_{mn(k)];c
 complex(dpc), allocatable :: chiw(:), chi2w(:) ! \chi_{II}^{abc}(-\omega,\omega,0)
 complex(dpc), allocatable :: etaw(:), eta2w(:) ! \eta_{II}^{abc}(-\omega,\omega,0)
 complex(dpc), allocatable :: sigmaw(:) ! \frac{i}{\omega} \sigma_{II}^{abc}(-\omega,\omega,0)
 complex(dpc) :: num1, num2, den1, den2, term1, term2
 complex(dpc) :: chi1, chi2_1, chi2_2
 complex(dpc), allocatable :: chi2(:) ! Second term that depends on the frequency ! (omega)
 complex(dpc), allocatable :: eta1(:) ! Second term that depends on the frequency ! (omega)
 complex(dpc), allocatable :: chi2tot(:)
 complex(dpc) :: eta1_1, eta1_2, eta2_1, eta2_2
 complex(dpc) :: sigma2_1, sigma1
 complex(dpc), allocatable :: symrmn(:,:,:) ! (m,l,n) = 1/2*(rml^b rln^c+rml^c rln^b)
 complex(dpc) :: symrmnl(3,3), symrlmn(3,3), symrmln(3,3)
!Parallelism
 integer :: my_rank, nproc
 integer,parameter :: master=0
 integer :: ierr
 integer :: my_k1, my_k2
 character(500) :: msg
 integer :: fout1,fout2,fout3,fout4,fout5,fout6,fout7

! *********************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

!calculate the constant
 zi=(0._dp,1._dp)
 idel=zi*brod
 const_au=-2._dp/(omega*dble(nsymcrys))
 !const_au=-2._dp/(omega)
 au2esu=5.8300348177d-8
 const_esu=const_au*au2esu
 ha2ev=13.60569172*2._dp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!5.8300348177d-8 : au2esu : bohr*c*10^4/4pi*2*ry2ev
!bohr: 5.2917ifc nlinopt.f907E-11
!c: 2.99792458   velocity of sound
!ry2ev: 13.60569172
!au2esu=(5.29177E-11*2.99792458*1.0E4)/(13.60569172*2)
!this const includes (e^3*hbar^3*hbar^3)/(vol*hbar^5*m_e^3)
!mass comes from converting P_mn to r_mn
!hbar^3 comes from converting all frequencies to energies in denominator
!hbar^3 comes from operator for momentum (hbar/i nabla)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!output file names
 fnam1=trim(fnam)//'-ChiSHGTotIm.out'
 fnam2=trim(fnam)//'-ChiSHGTotRe.out'
 fnam3=trim(fnam)//'-ChiSHGIm.out'
 fnam4=trim(fnam)//'-ChiSHGRe.out'
 fnam5=trim(fnam)//'-ChiSHGAbs.out'
 fnam6=trim(fnam)//'-ChiSHGImDec.out'
 fnam7=trim(fnam)//'-ChiSHGReDec.out'

!If there exists inversion symmetry exit with a mesg.
 tst=1.d-09
 do isym=1,nsymcrys
   t1=symcrys(1,1,isym)+1
   t2=symcrys(2,2,isym)+1
   t3=symcrys(3,3,isym)+1
!  test if diagonal elements are -1
   if (abs(t1).lt.tst.and.abs(t2).lt.tst.and.abs(t3).lt.tst) then
!    test if off-diagonal elements are zero
     if (abs(symcrys(1,2,isym)).lt.tst.and.abs(symcrys(1,3,isym)).lt.tst &
     .and.abs(symcrys(2,1,isym)).lt.tst.and.abs(symcrys(2,3,isym)).lt.tst.and.  &
     abs(symcrys(3,1,isym)).lt.tst.and.abs(symcrys(3,2,isym)).lt.tst) then
       write(std_out,*) '-----------------------------------------'
       write(std_out,*) '    the crystal has inversion symmetry   '
       write(std_out,*) '    the nl electro-optical susceptibility'
       write(std_out,*) '    is zero                              '
       write(std_out,*) '-----------------------------------------'
       MSG_ERROR("Aborting now")
     end if
   end if
 end do

!check polarisation
 if (v1.le.0.or.v2.le.0.or.v3.le.0.or.v1.gt.3.or.v2.gt.3.or.v3.gt.3) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in nonlinopt:                        '
   write(std_out,*) '    the polarisation directions incorrect    '
   write(std_out,*) '    1=x,  2=y  and 3=z                       '
   write(std_out,*) '---------------------------------------------'
   MSG_ERROR("Aborting now")
 end if

!number of energy mesh points
 if (nmesh.le.0) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in nonlinopt:                        '
   write(std_out,*) '    number of energy mesh points incorrect   '
   write(std_out,*) '    number has to be integer greater than 0  '
   write(std_out,*) '    nmesh*de = max energy for calculation    '
   write(std_out,*) '---------------------------------------------'
   MSG_ERROR("Aborting now")
 end if

!step in energy
 if (de.le.0._dp) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in nonlinopt:                        '
   write(std_out,*) '    energy step is incorrect                 '
   write(std_out,*) '    number has to real greater than 0.0      '
   write(std_out,*) '    nmesh*de = max energy for calculation    '
   write(std_out,*) '---------------------------------------------'
   MSG_ERROR("Aborting now")
 end if

!broadening
 if (brod.gt.0.009) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    ATTENTION: broadening is quite high      '
   write(std_out,*) '    ideally should be less than 0.005        '
   write(std_out,*) '---------------------------------------------'
 else if (brod.gt.0.015) then
   write(std_out,*) '----------------------------------------'
   write(std_out,*) '    ATTENTION: broadening is too high   '
   write(std_out,*) '    ideally should be less than 0.005   '
   write(std_out,*) '----------------------------------------'
 end if

!tolerance
 if (tol.gt.0.006) then
   write(std_out,*) '----------------------------------------'
   write(std_out,*) '    ATTENTION: tolerance is too high    '
   write(std_out,*) '    ideally should be less than 0.004   '
   write(std_out,*) '----------------------------------------'
 end if

 ! allocate local arrays
 ABI_MALLOC(enk,(nstval))
 ABI_MALLOC(delta,(nstval,nstval,3))
 ABI_MALLOC(rmnbc,(nstval,nstval,3,3))
 ABI_MALLOC(roverw,(nstval,nstval,3,3))
 ABI_MALLOC(rmna,(nstval,nstval,3))
 ABI_MALLOC(chiw,(nmesh))
 ABI_MALLOC(etaw,(nmesh))
 ABI_MALLOC(chi2w,(nmesh))
 ABI_MALLOC(eta2w,(nmesh))
 ABI_MALLOC(sigmaw,(nmesh))
 ABI_MALLOC(chi2,(nmesh))
 ABI_MALLOC(eta1,(nmesh))
 ABI_MALLOC(symrmn,(nstval,nstval,nstval))

!generate the symmetrizing tensor
 sym(:,:,:)=0._dp
 do isym=1,nsymcrys
   s(:,:)=symcrys(:,:,isym)
   do i=1,3
     do j=1,3
       do k=1,3
         sym(i,j,k)=sym(i,j,k)+(s(i,v1)*s(j,v2)*s(k,v3))
       end do
     end do
   end do
 end do


!initialise
 delta(:,:,:)=0._dp
 rmnbc(:,:,:,:)=0._dp
 chiw(:)=0._dp
 chi2w(:)=0._dp
 chi2(:) = 0._dp
 etaw(:)=0._dp
 eta2w(:)=0._dp
 sigmaw(:)=0._dp
 my_emin=HUGE(0._dp)
 my_emax=-HUGE(0._dp)

 ! Split work
 call xmpi_split_work(nkpt,comm,my_k1,my_k2)

! loop over kpts
 do ik=my_k1,my_k2
   write(std_out,*) "P-",my_rank,": ",ik,'of',nkpt
   do isp=1,nspin
     ! Calculate the scissor corrected energies and the energy window
     do ist1=1,nstval
       en = evalv(ist1,ik,isp)
       my_emin=min(my_emin,en)
       my_emax=max(my_emax,en)
       if(en > efermi) then
         en = en + sc
       end if
       enk(ist1) = en
     end do

!    calculate \Delta_nm and r_mn^a
     do istn=1,nstval
       en = enk(istn)
       do istm=1,nstval
         em = enk(istm)
         wmn = em - en
         delta(istn,istm,1:3)=pmat(istn,istn,ik,1:3,isp)-pmat(istm,istm,ik,1:3,isp)
         if(abs(wmn) < tol) then
           rmna(istm,istn,1:3) = 0._dp
         else
           rmna(istm,istn,1:3)=pmat(istm,istn,ik,1:3,isp)/wmn
         end if
       end do
     end do
!    calculate \r^b_mn;c
     do istm=1,nstval
       em = enk(istm)
       do istn=1,nstval
         en = enk(istn)
         wmn = em - en
         if (abs(wmn) < tol) then ! Degenerate energies
           rmnbc(istm,istn,:,:) = 0.0
           roverw(istm,istn,:,:) = 0.0
           cycle
         end if
         do ly = 1,3
           do lz = 1,3
             num1 = (rmna(istm,istn,ly)*delta(istm,istn,lz))+(rmna(istm,istn,lz)*delta(istm,istn,ly))
             den1 = wmn
             term1 = num1/den1
             term2 = 0._dp
             do istp=1,nstval
               ep = enk(istp)
               wmp = em - ep
               wpn = ep - en
               num2 = (wmp*rmna(istm,istp,ly)*rmna(istp,istn,lz))-(wpn*rmna(istm,istp,lz)*rmna(istp,istn,ly))
               den2 = wmn
               term2 = term2 + (num2/den2)
             end do
             rmnbc(istm,istn,ly,lz) = -term1-(zi*term2)
             roverw(istm,istn,ly,lz) = (rmnbc(istm,istn,ly,lz)/wmn) - (rmna(istm,istn,ly)/(wmn**2))*delta(istm,istn,lz)
           end do
         end do
       end do
     end do

!    initialise the factors
!    start the calculation
     do istn=1,nstval
       en=enk(istn)
       fn=occv(istn,ik,isp)
       if(do_antiresonant .and. en .ge. efermi) then
         cycle
       end if
       do istm=1,nstval
         em=enk(istm)
         if (do_antiresonant .and. em .le. efermi) then
           cycle
         end if
         wmn=em-en
         wnm=-wmn
         fm = occv(istm,ik,isp)
         fnm = fn - fm
         if(abs(wmn) > tol) then
           chi1 = 0._dp
           chi2(:) = 0._dp
           chi2_1 = 0._dp
           chi2_2 = 0._dp
           eta1(:) = 0._dp
           eta1_1 = 0._dp
           eta1_2 = 0._dp
           eta2_1 = 0._dp
           eta2_2 = 0._dp
           sigma1 = 0._dp
           sigma2_1 = 0._dp
           ! Three band terms
           do istl=1,nstval
             el=enk(istl)
             fl = occv(istl,ik,isp)
             wlm = el-em
             wln = el-en
             wnl = -wln
             wml = em-el
             fnl = fn-fl
             fml = fm-fl
             flm = -fml
             fln = -fnl
             do ly = 1,3
               do lz = 1,3
                 symrmnl(ly,lz) = 0.5_dp*(rmna(istm,istn,ly)*rmna(istn,istl,lz)+rmna(istm,istn,lz)*rmna(istn,istl,ly))
                 symrlmn(ly,lz) = 0.5_dp*(rmna(istl,istm,ly)*rmna(istm,istn,lz)+rmna(istl,istm,lz)*rmna(istm,istn,ly))
                 symrmln(ly,lz) = 0.5_dp*(rmna(istm,istl,ly)*rmna(istl,istn,lz)+rmna(istm,istl,lz)*rmna(istl,istn,ly))
               end do
             end do

             do lx = 1,3
               do ly = 1,3
                 do lz = 1,3
                   sigma1 = sigma1 + sym(lx,ly,lz)*(wnl*rmna(istl,istm,lx)*symrmnl(ly,lz)-wlm*rmna(istn,istl,lx)*symrlmn(ly,lz))
                   eta2_2 = eta2_2 + sym(lx,ly,lz)*fnm*rmna(istn,istm,lx)*symrmln(ly,lz)*(wml-wln)
                   if(abs(wln-wml) > tol) then
                     chi1 = chi1 + sym(lx,ly,lz)*(rmna(istn,istm,lx)*symrmln(ly,lz))/(wln-wml)
                   end if
                   eta1_1 = eta1_1 + sym(lx,ly,lz)*wln*rmna(istn,istl,lx)*symrlmn(ly,lz)
                   eta1_2 = eta1_2 - sym(lx,ly,lz)*wml*rmna(istl,istm,lx)*symrmnl(ly,lz)
                   if(abs(wnl-wmn) > tol) then
                     chi2_1 = chi2_1 - sym(lx,ly,lz)*(fnm*rmna(istl,istm,lx)*symrmnl(ly,lz)/(wnl-wmn))
                   end if
                   if(abs(wmn-wlm) > tol) then
                     chi2_2 = chi2_2 - sym(lx,ly,lz)*(fnm*rmna(istn,istl,lx)*symrlmn(ly,lz)/(wmn-wlm))
                   end if
                 end do
               end do
             end do
           end do

           ! Two band terms
           eta2_1 = 0.0_dp
           sigma2_1 = 0.0_dp
           do lx = 1,3
             do ly = 1,3
               do lz = 1,3
                 eta2_1 = eta2_1 + sym(lx,ly,lz)*fnm*rmna(istn,istm,lx)*0.5_dp &
&                    *(delta(istm,istn,ly)*rmna(istm,istn,lz)+delta(istm,istn,lz)*rmna(istm,istn,ly))
                 ! Correct version (Sipe 1993)
                 sigma2_1 = sigma2_1 + sym(lx,ly,lz)*fnm*rmna(istn,istm,lx)*0.5_dp &
&                    *(rmna(istm,istn,ly)*delta(istn,istm,lz)+rmna(istm,istn,lz)*delta(istn,istm,ly))

                 ! Incorrect version (Hughes 1996)
                 !sigma2_1 = fnm*delta(istn,istm,v1)*0.5_dp*(rmna(istm,istn,v2)*rmna(istn,istm,v3)+rmna(istm,istn,v3)*rmna(istn,istm,v2))
               end do
             end do
           end do
!
!          calculate over the desired energy mesh and sum over k-points
           do iw=1,nmesh
             w=(iw-1)*de+idel
             chi2w(iw) = chi2w(iw) + zi*wkpt(ik)*((2.0_dp*fnm*chi1/(wmn-2.0_dp*w)))*const_esu ! Inter(2w) from chi
             chiw(iw) = chiw(iw) + zi*wkpt(ik)*((chi2_1+chi2_2)/(wmn-w))*const_esu ! Inter(w) from chi
             eta2w(iw) = eta2w(iw) + zi*wkpt(ik)*(8.0_dp*(eta2_1/((wmn**2)*(wmn-2.0_dp*w))) &
&                 + 2.0_dp*eta2_2/((wmn**2)*(wmn-2.0_dp*w)))*const_esu ! Intra(2w) from eta
             etaw(iw) = etaw(iw) + zi*wkpt(ik)*((eta1_1 + eta1_2)*fnm/((wmn**2)*(wmn-w)))*const_esu ! Intra(w) from eta
             sigmaw(iw) = sigmaw(iw) + 0.5_dp*zi*wkpt(ik)*(fnm*sigma1/((wmn**2)*(wmn-w)) &
&                 + (sigma2_1/((wmn**2)*(wmn-w))))*const_esu ! Intra(1w) from sigma
           end do
         end if
       end do ! end loop over istn and istm
     end do
   end do ! spins
 end do ! k-points

 ! Collect info among the nodes
 call xmpi_min(my_emin,emin,comm,ierr)
 call xmpi_max(my_emax,emax,comm,ierr)

 call xmpi_sum(chiw,comm,ierr)
 call xmpi_sum(etaw,comm,ierr)
 call xmpi_sum(chi2w,comm,ierr)
 call xmpi_sum(eta2w,comm,ierr)
 call xmpi_sum(sigmaw,comm,ierr)

 ! Master writes the output
 if (my_rank == master) then

   if (ncid /= nctk_noid) then
     start4 = [1, 1, icomp, itemp]
     count4 = [2, nmesh, 1, 1]
     ABI_MALLOC(chi2tot, (nmesh))
     chi2tot = chiw + chi2w + etaw + eta2w + sigmaw
#ifdef HAVE_NETCDF
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "leo2_chi2tot"), c2r(chi2tot), start=start4, count=count4))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "leo2_chiw"), c2r(chiw), start=start4, count=count4))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "leo2_etaw"), c2r(etaw), start=start4, count=count4))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "leo2_chi2w"), c2r(chi2w), start=start4, count=count4))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "leo2_eta2w"), c2r(eta2w), start=start4, count=count4))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "leo2_sigmaw"), c2r(sigmaw), start=start4, count=count4))
#endif
     ABI_FREE(chi2tot)
   end if

   ! write output in SI units and esu (esu to SI(m/v)=(value_esu)*(4xpi)/30000)
   if (open_file(fnam1,msg,newunit=fout1,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam2,msg,newunit=fout2,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam3,msg,newunit=fout3,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam4,msg,newunit=fout4,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam5,msg,newunit=fout5,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam6,msg,newunit=fout6,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(fnam7,msg,newunit=fout7,action='WRITE',form='FORMATTED') /= 0) then
     MSG_ERROR(msg)
   end if
  !!write headers
   write(fout1, '(a,3i3)' ) ' #calculated the component:',v1,v2,v3
   write(fout1, '(a,es16.6)' ) ' #tolerance:',tol
   write(fout1, '(a,es16.6,a)' ) ' #broadening:',brod,'Ha'
   write(fout1, '(a,es16.6,a)' ) ' #scissors shift:',sc,'Ha'
   write(fout1, '(a,es16.6,a,es16.6,a)' ) ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout1, '(a)' )' # Energy      Tot-Im Chi(-w,w,0)  Tot-Im Chi(-w,w,0)'
   write(fout1, '(a)' )' # eV          *10^-7 esu        *10^-12 m/V SI units '
   write(fout1, '(a)' )' # '

   write(fout2, '(a,3i3)' ) ' #calculated the component:',v1,v2,v3
   write(fout2, '(a,es16.6)') ' #tolerance:',tol
   write(fout2, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout2, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout2, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout2, '(a)')' # Energy      Tot-Re Chi(-w,w,0)  Tot-Re Chi(-w,w,0)'
   write(fout2, '(a)')' # eV          *10^-7 esu        *10^-12 m/V SI units '
   write(fout2, '(a)')' # '

   write(fout3, '(a,3i3)') ' #calculated the component:',v1,v2,v3
   write(fout3, '(a,es16.6)') ' #tolerance:',tol
   write(fout3, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout3, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout3, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout3, '(a)')' # Energy(eV) Inter(2w) inter(1w) intra(2w) intra(1w)'
   write(fout3, '(a)')' # in esu'
   write(fout3, '(a)')' # '

   write(fout4, '(a,3i3)') ' #calculated the component:',v1,v2,v3
   write(fout4, '(a,es16.6)') ' #tolerance:',tol
   write(fout4, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout4, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout4, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout4, '(a)')' # Energy(eV) Inter(2w) inter(1w) intra(2w) intra(1w)'
   write(fout4, '(a)')' # in esu'
   write(fout4, '(a)')' # '

   write(fout5, '(a,3i3)') ' #calculated the component:',v1,v2,v3
   write(fout5, '(a,es16.6)') ' #tolerance:',tol
   write(fout5, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout5, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout5, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout5, '(a)')' # Energy(eV)  |TotChi(-w,w,0)|   |Tot Chi(-w,w,0)|'
   write(fout5, '(a)')' # eV          *10^-7 esu        *10^-12 m/V SI units '
   write(fout5, '(a)')' # '

   write(fout6, '(a,3i3)') ' #calculated the component:',v1,v2,v3
   write(fout6, '(a,es16.6)') ' #tolerance:',tol
   write(fout6, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout6, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout6, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout6, '(a)')' # Energy(eV) Chi(w) Eta(w) Sigma(w)'
   write(fout6, '(a)')' # in esu'
   write(fout6, '(a)')' # '

   write(fout7, '(a,3i3)') ' #calculated the component:',v1,v2,v3
   write(fout7, '(a,es16.6)') ' #tolerance:',tol
   write(fout7, '(a,es16.6,a)') ' #broadening:',brod,'Ha'
   write(fout7, '(a,es16.6,a)') ' #scissors shift:',sc,'Ha'
   write(fout7, '(a,es16.6,a,es16.6,a)') ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout7, '(a)')' # Energy(eV) Chi(w) Eta(w) Sigma(w)'
   write(fout7, '(a)')' # in esu'
   write(fout7, '(a)')' # '

   totim=0._dp
   totre=0._dp
   totabs=0._dp
   do iw=2,nmesh
     ene=(iw-1)*de
     ene=ene*ha2ev

     totim=aimag(chiw(iw)+chi2w(iw)+etaw(iw)+eta2w(iw)+sigmaw(iw))/1.d-7
     write(fout1,'(f15.6,2es15.6)') ene,totim,totim*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totim=0._dp

     totre=dble(chiw(iw)+chi2w(iw)+eta2w(iw)+etaw(iw)+sigmaw(iw))/1.d-7
     write(fout2,'(f15.6,2es15.6)') ene,totre,totre*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totre=0._dp

     write(fout3,'(f15.6,4es15.6)') ene,aimag(chi2w(iw))/1.d-7,aimag(chiw(iw))/1.d-7,     &
     aimag(eta2w(iw))/1.d-7,aimag(etaw(iw))/1.d-7+aimag(sigmaw(iw))/1.d-7

     write(fout4,'(f15.6,4es15.6)') ene,dble(chi2w(iw))/1.d-7,aimag(chiw(iw))/1.d-7,       &
     dble(eta2w(iw))/1.d-7,dble(etaw(iw))/1.d-7+dble(sigmaw(iw))/1.d-7

     totabs=abs(chiw(iw)+chi2w(iw)+etaw(iw)+eta2w(iw)+sigmaw(iw))/1.d-7
     write(fout5,'(f15.6,2es15.6)') ene,totabs,totabs*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totabs=0._dp

     write(fout6,'(f15.6,4es15.6)') ene,aimag(chi2w(iw)+chiw(iw))/1.d-7,      &
     aimag(eta2w(iw)+etaw(iw))/1.d-7,aimag(sigmaw(iw))/1.d-7

     write(fout7,'(f15.6,4es15.6)') ene,dble(chi2w(iw)+chiw(iw))/1.d-7,       &
     dble(eta2w(iw)+etaw(iw))/1.d-7,dble(sigmaw(iw))/1.d-7
   end do

   close(fout1)
   close(fout2)
   close(fout3)
   close(fout4)
   close(fout5)
   close(fout6)
   close(fout7)

   ! print information
   write(std_out,*) ' '
   write(std_out,*) 'information about calculation just performed:'
   write(std_out,*) ' '
   write(std_out,*) 'calculated the component:',v1,v2,v3 ,'of the nonlinear electro-optical susceptibility'
   write(std_out,*) 'tolerance:',tol
   if (tol.gt.0.008) write(std_out,*) 'ATTENTION: tolerance is too high'
   write(std_out,*) 'broadening:',brod,'Hartree'
   if (brod.gt.0.009) then
     write(std_out,*) ' '
     write(std_out,*) 'ATTENTION: broadening is quite high'
     write(std_out,*) ' '
   else if (brod.gt.0.015) then
     write(std_out,*) ' '
     write(std_out,*) 'ATTENTION: broadening is too high'
     write(std_out,*) ' '
   end if
   write(std_out,*) 'scissors shift:',sc,'Hartree'
   write(std_out,*) 'energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Hartree'

 end if

 ! deallocate local arrays
 ABI_FREE(enk)
 ABI_FREE(delta)
 ABI_FREE(rmnbc)
 ABI_FREE(roverw)
 ABI_FREE(rmna)
 ABI_FREE(chiw)
 ABI_FREE(chi2w)
 ABI_FREE(chi2)
 ABI_FREE(etaw)
 ABI_FREE(eta1)
 ABI_FREE(symrmn)
 ABI_FREE(eta2w)
 ABI_FREE(sigmaw)

end subroutine nonlinopt
!!***

!----------------------------------------------------------------------

END MODULE m_optic_tools
