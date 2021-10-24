!!****m* ABINIT/m_optic_tools
!! NAME
!! m_optic_tools
!!
!! FUNCTION
!!  Helper functions used in the optic code
!!
!! COPYRIGHT
!! Copyright (C) 2002-2021 ABINIT group (SSharma,MVer,VRecoules,TD,YG, NAP)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! COMMENTS
!!
!!  Right now the routine sums over the k-points. In future linear tetrahedron method might be useful.
!!
!!  Reference articles:
!!
!!      1. S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl, Phys. Rev. B {\bf 67} 165332 2003 [[cite:Sharma2003]]
!!      2. J. L. P. Hughes and J. E. Sipe, Phys. Rev. B {\bf 53} 10 751 1996 [[cite:Hughes1996]]
!!      3. S. Sharma and C. Ambrosch-Draxl, Physica Scripta T 109 2004 [[cite:Sharma2004]]
!!      4. J. E. Sipe and Ed. Ghahramani, Phys. Rev. B {\bf 48} 11 705 1993 [[cite:Sipe1993]]
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_optic_tools

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
 use m_numeric_tools,   only : c2r
 use m_io_tools,        only : open_file
 use m_crystal,         only : crystal_t

 implicit none

 private

 public :: pmat2cart
 public :: pmat_renorm
 public :: linopt           ! Compute dielectric function for semiconductors
 public :: nlinopt          ! Second harmonic generation susceptibility for semiconductors
 public :: linelop          ! Linear electro-optic susceptibility for semiconductors
 public :: nonlinopt        ! nonlinear electro-optic susceptibility for semiconductors

contains
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

subroutine pmat2cart(eigen11, eigen12, eigen13, mband, nkpt, nsppol, pmat, rprimd)

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
         pmat(iband2,iband1,ikpt,:,isppol) =  &
          rprim(:,1)*cmplx(eigen11(1,iband2,iband1,ikpt,isppol),eigen11(2,iband2,iband1,ikpt,isppol),kind=dp) &
         +rprim(:,2)*cmplx(eigen12(1,iband2,iband1,ikpt,isppol),eigen12(2,iband2,iband1,ikpt,isppol),kind=dp) &
         +rprim(:,3)*cmplx(eigen13(1,iband2,iband1,ikpt,isppol),eigen13(2,iband2,iband1,ikpt,isppol),kind=dp)
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
!!  fermie = Fermi level
!!  sc = scissor shift for conduction bands
!!  eig = ground state eigenvalues
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

subroutine pmat_renorm(fermie, eig, mband, nkpt, nsppol, pmat, sc)

!Arguments -----------------------------------------------
!scalars
 integer, intent(in) :: nsppol
 integer, intent(in) :: nkpt
 integer, intent(in) :: mband
 real(dp), intent(in) :: fermie
 real(dp), intent(in) :: sc
!arrays
 real(dp), intent(in) :: eig(mband,nkpt,nsppol)
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
       e1 = eig(iband1,ikpt,isppol)
       if (e1 > fermie) cycle
       do iband2=1,mband ! conduction states
         e2 = eig(iband2,ikpt,isppol)
         if (e2 < fermie) cycle
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
!!  nband_sum=Number of bands included in the sum. Must be <= mband
!!  pmat(mband,mband,nkpt,3,nsppol)=momentum matrix elements in cartesian coordinates(complex)
!!  v1,v2=desired component of the dielectric function(integer) 1=x,2=y,3=z
!!  nmesh=desired number of energy mesh points(integer)
!!  de=desired step in energy(real); nmesh*de=maximum energy
!!  sc=scissors shift in Ha(real)
!!  brod=broadening in Ha(real)
!!  fnam=root for filename that will contain the output filename will be trim(fnam)//'-linopt.out'
!!  ncid=Netcdf id to save output data.
!!  prtlincompmatrixelements=if set to 1, the matrix elements are dumped in the _OPTIC.nc file for post processing.
!!
!! SIDE EFFECTS
!!  Dielectric function for semiconductors, on a desired energy mesh and for a desired
!!  direction of polarisation is written to file.
!!  The output is in a file named trim(fnam)//'-linopt.out' and contains
!!  Im(\epsilon_{v1v2}(\omega), Re(\epsilon_{v1v2}(\omega) and abs(\epsilon_{v1v2}(\omega).
!!
!!  If 'prtlincompmatrixelements' is set to 1, the matrix elements and other quantities used to build
!!  the chi tensor are stored in the _OPTIC.nc file as well. This includes the matrix elements,
!!  the occupations, the renormalized but unshifted eigenvalues and the kpts weights.
!!
!!  Comment:
!!    Right now the routine sums over the kpoints. In future linear tetrahedron method should be useful.
!!
!! PARENTS
!!      optic
!!
!! CHILDREN
!!      xmpi_max,xmpi_min,xmpi_split_work,xmpi_sum
!!
!! SOURCE

subroutine linopt(icomp, itemp, nband_sum, cryst, ks_ebands, EPBSt, pmat, &
  v1, v2, nmesh, de, sc, brod, fnam, ncid, prtlincompmatrixelements, comm)

!Arguments ------------------------------------
integer, intent(in) :: icomp,itemp,nband_sum, ncid
type(crystal_t), intent(in) :: cryst
type(ebands_t),intent(in) :: ks_ebands,EPBSt
complex(dpc), intent(in) :: pmat(ks_ebands%mband, ks_ebands%mband, ks_ebands%nkpt, 3, ks_ebands%nsppol)
integer, intent(in) :: v1, v2, nmesh
real(dp), intent(in) :: de, sc, brod
character(len=*), intent(in) :: fnam
integer, intent(in) :: comm
integer, intent(in) :: prtlincompmatrixelements

!Local variables -------------------------
integer,parameter :: master=0
integer :: isp,i,j,isym,lx,ly,ik,ist1,ist2,iw,nkpt
integer :: my_rank, nproc, my_k1, my_k2, ierr, fout1, mband, nsppol
#ifdef HAVE_NETCDF
integer :: ncerr
#endif
logical :: do_linewidth
real(dp) :: deltav1v2, ha2ev, tmpabs, renorm_factor,emin,emax
real(dp) :: ene,abs_eps,re_eps
complex(dpc) :: e1,e2,e12, e1_ep,e2_ep,e12_ep, b11,b12, ieta, w
character(len=fnlen) :: fnam1
character(len=500) :: msg
! allocatable arrays
real(dp) :: s(3,3),sym(3,3)
real(dp), allocatable :: im_refract(:),re_refract(:)
complex(dpc), allocatable :: chi(:,:), matrix_elements(:,:,:,:), renorm_eigs(:,:,:), eps(:)

! *********************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)
 nkpt = ks_ebands%nkpt
 nsppol = ks_ebands%nsppol
 mband = ks_ebands%mband
 ABI_CHECK(nband_sum <= mband, "nband_sum <= mband")

 if (my_rank == master) then
   ! check polarisation
   if (v1.le.0.or.v2.le.0.or.v1.gt.3.or.v2.gt.3) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    Error in linopt:                         '
     write(std_out,*) '    the polarisation directions incorrect    '
     write(std_out,*) '    1=x and 2=y and 3=z                      '
     write(std_out,*) '---------------------------------------------'
     ABI_ERROR("Aborting now")
   end if
   ! number of energy mesh points
   if (nmesh.le.0) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    Error in linopt:                         '
     write(std_out,*) '    number of energy mesh points incorrect   '
     write(std_out,*) '    number has to integer greater than 0     '
     write(std_out,*) '    nmesh*de = max energy for calculation    '
     write(std_out,*) '---------------------------------------------'
     ABI_ERROR("Aborting now")
   end if
   ! step in energy
   if (de.le.zero) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    Error in linopt:                         '
     write(std_out,*) '    energy step is incorrect                 '
     write(std_out,*) '    number has to real greater than 0.0      '
     write(std_out,*) '    nmesh*de = max energy for calculation    '
     write(std_out,*) '---------------------------------------------'
     ABI_ERROR("Aborting now")
   end if
   ! broadening
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
   ! fermi energy
   if(ks_ebands%fermie<-1.0d4) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    ATTENTION: Fermi energy seems extremely  '
     write(std_out,*) '    low                                      '
     write(std_out,*) '---------------------------------------------'
   end if
   ! scissors operator
   if (sc.lt.zero) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    Error in linopt:                         '
     write(std_out,*) '    scissors shift is incorrect              '
     write(std_out,*) '    number has to be greater than 0.0      '
     write(std_out,*) '---------------------------------------------'
     ABI_ERROR("Aborting now")
   end if
 end if

 do_linewidth = allocated(EPBSt%linewidth)
! TODO: activate this, and remove do_linewidth - always add it in even if 0.
! if (.not. allocated(EPBSt%linewidth)) then
!   ABI_MALLOC(EPBSt%linewidth, (1, mband, my_k2-my_k1+1, nsppol))
!   EPBSt%linewidth = zero
! end if

 ABI_MALLOC(chi, (nmesh, nsppol))
 ABI_MALLOC(eps, (nmesh))
 ABI_MALLOC(im_refract, (nmesh))
 ABI_MALLOC(re_refract, (nmesh))
 ieta=(zero, 1._dp)*brod
 renorm_factor=1._dp/(cryst%ucvol*dble(cryst%nsym))
 ha2ev=13.60569172*2._dp

 ! output file names
 fnam1=trim(fnam)//'-linopt.out'

 ! construct symmetrisation tensor
 sym = zero
 do isym=1,cryst%nsym
   s(:,:)=cryst%symrel_cart(:,:,isym)
   do i=1,3
     do j=1,3
       sym(i,j)=sym(i,j)+s(i,v1)*s(j,v2)
     end do
   end do
 end do

 ! calculate the energy window
 emin=zero
 emax=zero
 do ik=1,nkpt
   do isp=1,nsppol
     do ist1=1,nband_sum
       emin=min(emin,EPBSt%eig(ist1,ik,isp))
       emax=max(emax,EPBSt%eig(ist1,ik,isp))
     end do
   end do
 end do

 ! Split work
 call xmpi_split_work(nkpt,comm,my_k1,my_k2)
 ! if we print matrix elements, allocate full arrays for each process
 ! this is not optimized memory-wise since we could just allocate what is needed
 ! however we would need to write all data using mpi-io.
 if (prtlincompmatrixelements == 1) then
   ABI_CALLOC(matrix_elements, (mband, mband, nkpt, nsppol))
   ABI_CALLOC(renorm_eigs, (mband, nkpt, nsppol))
 endif

 ! start calculating linear optical response
 chi(:,:)=zero
 do isp=1,nsppol
   do ik=my_k1,my_k2
     write(std_out,*) "P-",my_rank,": ",ik,'of',nkpt
     do ist1=1,nband_sum
       e1=ks_ebands%eig(ist1,ik,isp)
       e1_ep=EPBSt%eig(ist1,ik,isp)
       ! TODO: unless memory is a real issue, should set lifetimes to 0 and do this sum systematically
       ! instead of putting an if statement in a loop! See above
       if(do_linewidth) then
         e1_ep = e1_ep + EPBSt%linewidth(1,ist1,ik,isp)*(0.0_dp,1.0_dp)
       end if
       do ist2=1,nband_sum
         e2=ks_ebands%eig(ist2,ik,isp)
         e2_ep=EPBSt%eig(ist2,ik,isp)
         if(do_linewidth) then
           e2_ep = e2_ep - EPBSt%linewidth(1,ist2,ik,isp)*(0.0_dp,1.0_dp)
         end if
         if (ist1.ne.ist2) then
           ! scissors correction of momentum matrix
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
           ! e12=e1-e2-sc
           b11=zero
           ! symmetrization of momentum matrix
           do lx=1,3
             do ly=1,3
               b11=b11+(sym(lx,ly)*pmat(ist1,ist2,ik,lx,isp)* &
               conjg(pmat(ist1,ist2,ik,ly,isp)))
             end do
           end do
           b12=b11*renorm_factor*(1._dp/(e12**2))
           ! store data for printing if necessary
           if (prtlincompmatrixelements == 1) then
             matrix_elements(ist1,ist2,ik,isp) = b12
             renorm_eigs(ist1,ik,isp) = e1_ep
             renorm_eigs(ist2,ik,isp) = e2_ep
           endif
           ! calculate on the desired energy grid
           do iw=2,nmesh
             w=(iw-1)*de+ieta
             chi(iw,isp)=chi(iw,isp)+(ks_ebands%wtk(ik)*(ks_ebands%occ(ist1,ik,isp)-ks_ebands%occ(ist2,ik,isp))* &
             (b12/(-e12_ep-w)))
           end do ! frequencies
         end if
       end do  ! states 2
     end do  ! states 1
   end do ! k points
 end do ! spin

 call xmpi_sum(chi,comm,ierr)
 if (prtlincompmatrixelements == 1) then
   ! gather all data to main process in order to write them using a single process
   ! in the netcdf file. This could be avoided by doing mpiio.
   call xmpi_sum(matrix_elements,comm,ierr)
   call xmpi_sum(renorm_eigs,comm,ierr)
 endif

 ! calculate epsilon
 eps(1) = zero
 deltav1v2=zero; if (v1 == v2) deltav1v2=one
 do iw=2,nmesh
   eps(iw)=deltav1v2+4._dp*pi*sum(chi(iw,:))
 end do

 if (my_rank == master) then
   !  open the output files
   if (open_file(fnam1,msg,newunit=fout1,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   ! write output
   write(fout1, '(a,2i3,a)' )' #calculated the component:',v1,v2,'  of dielectric function'
   write(std_out,*) 'calculated the component:',v1,v2,'  of dielectric function'
   write(fout1, '(a,2es16.6)' ) ' #broadening:', real(ieta),aimag(ieta)
   write(std_out,*) ' with broadening:',ieta
   write(fout1, '(a,es16.6)' ) ' #scissors shift:',sc
   write(std_out,*) 'and scissors shift:',sc
   write(fout1, '(a,es16.6,a,es16.6,a)' ) ' #energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(std_out,*) 'energy window:',(emax-emin)*ha2ev,'eV',(emax-emin),'Ha'
   write(fout1,*)
   if(nsppol==1)write(fout1, '(a)' ) ' # Energy(eV)         Im(eps(w))'
   if(nsppol==2)write(fout1, '(a)' ) ' # Energy(eV)         Im(eps(w))         Spin up       Spin down '
   do iw=2,nmesh
     ene=(iw-1)*de*ha2ev
     if(nsppol==1)write(fout1, '(2es16.6)' ) ene,aimag(eps(iw))
     if(nsppol==2)write(fout1, '(4es16.6)' ) ene,aimag(eps(iw)),4._dp*pi*aimag(chi(iw,1)),4._dp*pi*aimag(chi(iw,2))
   end do
   write(fout1,*)
   write(fout1,*)
   if(nsppol==1)write(fout1, '(a)' ) ' # Energy(eV)         Re(eps(w))'
   if(nsppol==2)write(fout1, '(a)' ) ' # Energy(eV)         Re(eps(w))         Spin up       Spin down    +delta(diag) '
   do iw=2,nmesh
     ene=(iw-1)*de*ha2ev
     if(nsppol==1)write(fout1, '(2es16.6)' ) ene,dble(eps(iw))
     if(nsppol==2)write(fout1, '(5es16.6)' ) ene,dble(eps(iw)),4._dp*pi*dble(chi(iw,1)),4._dp*pi*dble(chi(iw,2)),deltav1v2
   end do
   write(fout1,*)
   write(fout1,*)
   write(fout1, '(a)' )' # Energy(eV)         abs(eps(w))'
   do iw=2,nmesh
     ene=(iw-1)*de*ha2ev
     abs_eps=abs(eps(iw))
     re_eps=dble(eps(iw))
     write(fout1, '(2es16.6)' ) ene,abs_eps
     re_refract(iw)=sqrt(half*(abs_eps+re_eps))
     im_refract(iw)=sqrt(half*(abs_eps-re_eps))
   end do
   write(fout1,*)
   write(fout1,*)
   write(fout1, '(a)' )' # Energy(eV)         Im(refractive index(w)) aka kappa'
   do iw=2,nmesh
     ene=(iw-1)*de*ha2ev
     write(fout1, '(2es16.6)' ) ene,im_refract(iw)
   end do
   write(fout1,*)
   write(fout1,*)
   write(fout1, '(a)' )' # Energy(eV)         Re(refractive index(w)) aka n'
   do iw=2,nmesh
     ene=(iw-1)*de*ha2ev
     write(fout1, '(2es16.6)' ) ene,re_refract(iw)
   end do
   write(fout1,*)
   write(fout1,*)
   write(fout1, '(a)' )' # Energy(eV)         Reflectivity(w) from vacuum, at normal incidence'
   do iw=2,nmesh
     ene=(iw-1)*de*ha2ev
     write(fout1, '(2es16.6)' ) ene, ((re_refract(iw)-one)**2+im_refract(iw)**2)/((re_refract(iw)+one)**2+im_refract(iw)**2)
   end do
   write(fout1,*)
   write(fout1,*)
   write(fout1, '(a)' )' # Energy(eV)         absorption coeff (in 10^6 m-1) = omega Im(eps) / c n(eps)'
   do iw=2,nmesh
     ene=(iw-1)*de
     tmpabs=zero
     if ( re_refract(iw) > tol10 ) then
       tmpabs = aimag(eps(iw))*ene / re_refract(iw) / Sp_Lt / Bohr_meter * 1.0d-6
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
   if (prtlincompmatrixelements == 1) then
     ! write matrix elements and other quantities used to build the chi tensor.
     write(std_out, '(a)') 'Writing linopt matrix elements in _OPTIC.nc file.'
     ncerr = nf90_put_var(ncid, nctk_idname(ncid, "linopt_matrix_elements"), c2r(matrix_elements),&
                          start=[1, 1, 1, 1, 1, icomp, itemp])
     NCF_CHECK(ncerr)
     ncerr = nf90_put_var(ncid, nctk_idname(ncid, "linopt_renorm_eigs"), c2r(renorm_eigs), start=[1, 1, 1, 1])
     NCF_CHECK(ncerr)

     ! write occupations and kpt weights
     ncerr = nf90_put_var(ncid, nctk_idname(ncid, "linopt_occupations"), ks_ebands%occ, start=[1, 1, 1])
     NCF_CHECK(ncerr)
     ncerr = nf90_put_var(ncid, nctk_idname(ncid, "linopt_wkpts"), ks_ebands%wtk, start=[1])
     NCF_CHECK(ncerr)
     write(std_out, '(a)') 'Writing linopt matrix elements done.'
   endif

#endif
 end if ! rank == master

 ABI_FREE(chi)
 ABI_FREE(eps)
 ABI_FREE(im_refract)
 ABI_FREE(re_refract)

 ABI_SFREE(matrix_elements)
 ABI_SFREE(renorm_eigs)

end subroutine linopt
!!***

!----------------------------------------------------------------------

!!****f* m_optic_tools/nlinopt
!! NAME
!! nlinopt
!!
!! FUNCTION
!! Compute second harmonic generation susceptibility for semiconductors
!!
!! INPUTS
!!  icomp=Sequential index associated to computed tensor components (used for netcdf output)
!!  itemp=Temperature index (used for netcdf output)
!!  nband_sum=Number of bands included in the sum. Must be <= mband
!!  fermie = Fermi energy in Ha(real)
!!  pmat(mband,mband,nkpt,3,nsppol) = momentum matrix elements in cartesian coordinates(complex)
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

subroutine nlinopt(icomp, itemp, nband_sum, cryst, ks_ebands, pmat, &
                   v1, v2, v3, nmesh, de, sc, brod, tol, fnam, ncid, comm)

!Arguments ------------------------------------
integer, intent(in) :: icomp, itemp, nband_sum, ncid
type(crystal_t),intent(in) :: cryst
type(ebands_t),intent(in) :: ks_ebands
complex(dpc), intent(in) :: pmat(ks_ebands%mband, ks_ebands%mband, ks_ebands%nkpt, 3, ks_ebands%nsppol)
integer, intent(in) :: v1, v2, v3, nmesh, comm
real(dp), intent(in) :: de, sc, brod, tol
character(len=*), intent(in) :: fnam

!Local variables -------------------------
integer,parameter :: master=0
integer :: iw, mband,i,j,k,lx,ly,lz
integer :: isp,isym,ik,ist1,ist2,istl,istn,istm
integer :: my_rank, nproc, my_k1, my_k2, ierr
integer :: fout1,fout2,fout3,fout4,fout5,fout6,fout7
real(dp) :: f1,f2,f3, ha2ev
real(dp) :: t1,t2,t3
real(dp) :: ene,totre,totabs,totim
real(dp) :: e1,e2,el,en,em,emin,emax,my_emin,my_emax
real(dp) :: const_esu,const_au,au2esu,wmn,wnm,wln,wnl,wml,wlm
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
complex(dpc), allocatable :: px(:,:,:,:,:), py(:,:,:,:,:), pz(:,:,:,:,:)
complex(dpc), allocatable :: delta(:,:,:), inter2w(:), inter1w(:)
complex(dpc), allocatable :: intra2w(:), intra1w(:), intra1wS(:),chi2tot(:)

! *********************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)
 mband = ks_ebands%mband

!calculate the constant
 zi=(0._dp,1._dp)
 idel=zi*brod
 const_au=-2._dp/(cryst%ucvol*dble(cryst%nsym))
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
   ! If there exists inversion symmetry exit with a message.
   if (cryst%idx_spatial_inversion() /= 0) then
     write(std_out,*) '-------------------------------------------'
     write(std_out,*) '    The crystal has inversion symmetry     '
     write(std_out,*) '    The SHG susceptibility is zero         '
     write(std_out,*) '    Action : set num_nonlin_comp to zero   '
     write(std_out,*) '-------------------------------------------'
     ABI_ERROR("Aborting now")
   end if
   ! check polarisation
   if (v1.le.0.or.v2.le.0.or.v3.le.0.or.v1.gt.3.or.v2.gt.3.or.v3.gt.3) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    Error in nlinopt:                        '
     write(std_out,*) '    Incorrect polarisation directions        '
     write(std_out,*) '    1=x,  2=y  and 3=z                       '
     write(std_out,*) '    Action : check your input file,          '
     write(std_out,*) '    use only 1, 2 or 3 to define directions  '
     write(std_out,*) '---------------------------------------------'
     ABI_ERROR("Aborting now")
   end if
   !number of energy mesh points
   if (nmesh.le.0) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    Error in nlinopt:                        '
     write(std_out,*) '    number of energy mesh points incorrect   '
     write(std_out,*) '    number has to be integer greater than 0  '
     write(std_out,*) '    nmesh*de = max energy for calculation    '
     write(std_out,*) '---------------------------------------------'
     ABI_ERROR("Aborting now")
   end if
   !step in energy
   if (de.le.zero) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    Error in nlinopt:                        '
     write(std_out,*) '    energy step is incorrect                 '
     write(std_out,*) '    number has to real greater than 0.0      '
     write(std_out,*) '    nmesh*de = max energy for calculation    '
     write(std_out,*) '---------------------------------------------'
     ABI_ERROR("Aborting now")
   end if
   !broadening
   if (brod.gt.0.009) then
     write(std_out,*) '---------------------------------------------'
     write(std_out,*) '    WARNING : broadening is quite high       '
     write(std_out,*) '    ideally should be less than 0.005        '
     write(std_out,*) '---------------------------------------------'
   else if (brod.gt.0.015) then
     write(std_out,*) '----------------------------------------'
     write(std_out,*) '    WARNING : broadening is too high    '
     write(std_out,*) '    ideally should be less than 0.005   '
     write(std_out,*) '----------------------------------------'
   end if
   !tolerance
   if (tol.gt.0.006) then
     write(std_out,*) '----------------------------------------'
     write(std_out,*) '    WARNING : tolerance is too high     '
     write(std_out,*) '    ideally should be less than 0.004   '
     write(std_out,*) '----------------------------------------'
   end if
 end if

 !allocate local arrays
 ABI_MALLOC(px, (mband, mband, 3, 3, 3))
 ABI_MALLOC(py, (mband, mband, 3, 3, 3))
 ABI_MALLOC(pz,(mband,mband, 3, 3, 3))
 ABI_MALLOC(inter2w, (nmesh))
 ABI_MALLOC(inter1w, (nmesh))
 ABI_MALLOC(intra2w, (nmesh))
 ABI_MALLOC(intra1w, (nmesh))
 ABI_MALLOC(intra1wS, (nmesh))
 ABI_MALLOC(delta, (mband, mband, 3))

 !generate the symmetrizing tensor
 sym = zero
 do isym=1,cryst%nsym
   s(:,:)=cryst%symrel_cart(:,:,isym)
   do i=1,3
     do j=1,3
       do k=1,3
         sym(i,j,k)=sym(i,j,k)+(s(i,v1)*s(j,v2)*s(k,v3))
       end do
     end do
   end do
 end do
 ! Disable symmetries for now
 !sym(:,:,:) = zero
 !sym(v1,v2,v3) = nsym

 ! Split work
 call xmpi_split_work(ks_ebands%nkpt, comm, my_k1, my_k2)

 ! initialise
 inter2w(:)=zero
 inter1w(:)=zero
 intra2w(:)=zero
 intra1w(:)=zero
 intra1wS(:)=zero
 delta(:,:,:)=zero

 my_emin=HUGE(zero)
 my_emax=-HUGE(zero)

 ! loop over kpts
 do ik=my_k1,my_k2
   write(std_out,*) "P-",my_rank,": ",ik,'of',ks_ebands%nkpt
   ! loop over spins
   do isp=1,ks_ebands%nsppol
     !  loop over states
     do ist1=1,mband
       e1 = ks_ebands%eig(ist1,ik,isp)
       if (e1.lt.ks_ebands%fermie) then   ! ist1 is a valence state
         do ist2=1,mband
           e2 = ks_ebands%eig(ist2,ik,isp)
           if (e2.gt.ks_ebands%fermie) then ! ist2 is a conduction state
             ! symmetrize the momentum matrix elements
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
             end do ! end loop over states
           end if
         end do
       end if
     end do

     ! calculate the energy window and \Delta_nm
     do ist1=1,mband
       my_emin=min(my_emin, ks_ebands%eig(ist1,ik,isp))
       my_emax=max(my_emax, ks_ebands%eig(ist1,ik,isp))
       do ist2=1,mband
         delta(ist1,ist2,1:3)=pmat(ist1,ist1,ik,1:3,isp)-pmat(ist2,ist2,ik,1:3,isp)
       end do
     end do
     ! initialise the factors
     ! factors are named according to the Ref. article 2.
     b111=zero
     b121=zero
     b131=zero
     b112=zero
     b122=zero
     b132=zero
     b113=zero
     b123=zero
     b133=zero
     b211=zero
     b221=zero
     b212=zero
     b222=zero
     b213=zero
     b223=zero
     b231=zero
     b241=zero
     b242=zero
     b243=zero
     b311=zero
     b312=zero
     b313=zero
     b331=zero
     ! start the calculation
     do istn=1,mband
       en=ks_ebands%eig(istn,ik,isp)
       if (en.lt.ks_ebands%fermie) then  ! istn is a valence state
         do istm=1,mband
           em=ks_ebands%eig(istm,ik,isp)
           if (em.gt.ks_ebands%fermie) then  ! istm is a conduction state
             em = em + sc ! Should add the scissor to conduction energies
             wmn=em-en
             wnm=-wmn
             ! calculate the matrix elements for two band intraband term
             mat2w_tra=zero
             mat1w3_tra=zero
             do lx=1,3
               do ly=1,3
                 do lz=1,3
                   mat2w_tra=mat2w_tra+px(istn,istm,lx,ly,lz)*pmat(istm,istn,ik,lz,isp)    &
                   *delta(istm,istn,ly)
                   mat1w3_tra=mat1w3_tra+px(istn,istm,lx,ly,lz)*pmat(istm,istn,ik,ly,isp)  &
                   *delta(istm,istn,lz)
                   ! NOTE:: lx to ly m to n in pmat matrices respectively
                   ! Changes are made so that this (b3) term is according to paper
                   ! [[cite:Sipe1993]] (Ref. 4) rather than [[cite:Hughes1996]] (Ref 2) in which this term is incorrect
                 end do
               end do
             end do
             b331=mat1w3_tra/wnm
             b11=zero
             b12_13=zero
             b24=zero
             b31_32=zero
             b21_22=zero

             b231=8._dp*mat2w_tra/wmn
             b331=mat1w3_tra/(wnm)
             ! istl < istn
             do istl=1,istn-1  ! istl is a valence state below istn
               el=ks_ebands%eig(istl,ik,isp)
               wln=el-en  ! do not add sc to the valence bands!
               wml=em-el
               wnl=-wln
               wlm=-wml
               ! calculate the matrix elements for three band terms
               mat2w=zero
               mat1w1=zero
               mat1w2=zero
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
               b221=zero
               b211=mat1w1/wml
               b241=-mat2w/wml
               b311=mat1w2/wlm

               if (abs(wln).gt.tol) then
                 b111=b111/wln
                 b121=b121/wln
                 b131=b131/wln
                 b221=mat1w2/wln
                 b241=b241+(mat2w/wln)
                 b311=b311+(mat1w1/wln)
               else
                 b111=zero
                 b121=zero
                 b131=zero
                 b221=zero
               end if
               t1=wln-wnm
               if (abs(t1).gt.tol) then
                 b131=b131/t1
               else
                 b131=zero
               end if
               b11=b11-2._dp*b111
               b12_13=b12_13+b121+b131
               b21_22=b21_22-b211+b221
               b24=b24+2._dp*b241
               b31_32=b31_32+b311
             end do ! istl

             ! istn < istl < istm
             do istl=istn+1,istm-1
               el=ks_ebands%eig(istl,ik,isp)
               ! calculate the matrix elements for three band terms
               mat2w=zero
               mat1w1=zero
               mat1w2=zero
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
               if (el.lt.ks_ebands%fermie) then
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
               b112=zero
               b122=mat1w1*(1._dp/(wnm+wlm))
               b132=mat1w2*(1._dp/(wnm+wnl))
               b242=zero
               b222=zero
               b212=zero
               b312=zero
               if (abs(wnl).gt.tol) then
                 b112=mat2w/wln
                 b122=b122/wnl
                 b132=b132/wnl
                 b242=mat2w/wln
                 b222=mat1w2/wln
                 b312=mat1w1/wln
               else
                 b122=zero
                 b132=zero
               end if
               if (abs(wlm).gt.tol) then
                 b112=b112/wml
                 b122=b122/wlm
                 b132=b132/wlm
                 b242=b242-(mat2w/wml)
                 b212=mat1w1/wml
                 b312=b312+(mat1w2/wlm)
               else
                 b112=zero
                 b122=zero
                 b132=zero
                 b212=zero
               end if
               t1=wlm-wnl
               if (abs(t1).gt.tol) then
                 b112=b112/t1
               else
                 b112=zero
               end if
               b11=b11+2._dp*b112
               b12_13=b12_13-b122+b132
               b24=b24+2._dp*b242
               b21_22=b21_22-b212+b222
               b31_32=b31_32+b312
             end do ! istl

             ! istl > istm    !
             do istl=istm+1,mband
               el=ks_ebands%eig(istl,ik,isp)+sc
               wln=el-en
               wnl=-wln
               wml=em-el
               wlm=-wml
               ! calculate the matrix elements for three band terms
               mat2w=zero
               mat1w1=zero
               mat1w2=zero
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

               b113=mat2w*(1._dp/(wnl+wml))*(1._dp/wnl)
               b123=mat1w1*(1._dp/wnl)
               b133=mat1w2*(1._dp/wnl)*(1._dp/(wnl+wnm))
               b243=mat2w/wln
               b223=mat1w2/wln
               b213=zero
               b313=-1._dp*mat1w1/wnl
               if (abs(wml).gt.tol) then
                 b113=b113/wml
                 b123=b123/wml
                 b133=b133/wml
                 b243=b243-(mat2w/wml)
                 b213=mat1w1/wml
                 b313=b313+(mat1w2/wlm)
               else
                 b113=zero
                 b123=zero
                 b133=zero
               end if
               t1=wnm-wml
               if (abs(t1).gt.tol) then
                 b123=b123/t1
               else
                 b123=zero
               end if
               b11=b11+2._dp*b113
               b12_13=b12_13+b123-b133
               b21_22=b21_22-b213+b223
               b24=b24+2._dp*b243
               b31_32=b31_32+b313
             end do ! istl

             b11=b11*zi*(1._dp/wnm)*const_esu
             b12_13=b12_13*zi*(1._dp/wnm)*const_esu
             b24=(b24+b231)*zi*(1._dp/(wnm**3))*const_esu
             b21_22=(b21_22)*zi*(1._dp/(wnm**3))*const_esu
             b31_32=(b31_32-b331)*zi*(1._dp/(wmn**3))*const_esu*0.5_dp
             ! calculate over the desired energy mesh and sum over k-points
             do iw=1,nmesh
               w=(iw-1)*de+idel
               inter2w(iw)=inter2w(iw)+(ks_ebands%wtk(ik)*(b11/(wmn-2._dp*w))) ! Inter(2w) from chi
               inter1w(iw)=inter1w(iw)+(ks_ebands%wtk(ik)*(b12_13/(wmn-w))) ! Inter(1w) from chi
               intra2w(iw)=intra2w(iw)+(ks_ebands%wtk(ik)*(b24/(wmn-2._dp*w))) ! Intra(2w) from eta
               intra1w(iw)=intra1w(iw)+(ks_ebands%wtk(ik)*((b21_22)/(wmn-w))) ! Intra(1w) from eta
               intra1wS(iw)=intra1wS(iw)+(ks_ebands%wtk(ik)*((b31_32)/(wmn-w))) ! Intra(1w) from sigma
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
     ABI_ERROR(msg)
   end if
   if (open_file(fnam2,msg,newunit=fout2,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   if (open_file(fnam3,msg,newunit=fout3,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   if (open_file(fnam4,msg,newunit=fout4,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   if (open_file(fnam5,msg,newunit=fout5,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   if (open_file(fnam6,msg,newunit=fout6,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   if (open_file(fnam7,msg,newunit=fout7,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   ! write headers
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

   totim=zero
   totre=zero
   totabs=zero
   do iw=2,nmesh
     ene=(iw-1)*de
     ene=ene*ha2ev

     totim=aimag(inter2w(iw)+inter1w(iw)+intra2w(iw)+intra1w(iw)+intra1wS(iw))/1.d-7
     write(fout1,'(f15.6,2es15.6)') ene,totim,totim*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totim=zero

     totre=dble(inter2w(iw)+inter1w(iw)+intra2w(iw)+intra1w(iw)+intra1wS(iw))/1.d-7
     write(fout2,'(f15.6,2es15.6)') ene,totre,totre*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totre=zero

     write(fout3,'(f15.6,4es15.6)') ene,aimag(inter2w(iw))/1.d-7,      &
     aimag(inter1w(iw))/1.d-7,aimag(intra2w(iw))/1.d-7, aimag(intra1w(iw)+intra1wS(iw))/1.d-7

     write(fout4,'(f15.6,4es15.6)') ene,dble(inter2w(iw))/1.d-7,       &
     dble(inter1w(iw))/1.d-7,dble(intra2w(iw))/1.d-7,dble(intra1w(iw)+intra1wS(iw))/1.d-7

     totabs=abs(inter2w(iw)+inter1w(iw)+intra2w(iw)+intra1w(iw)+intra1wS(iw))/1.d-7
     write(fout5,'(f15.6,2es15.6)') ene,totabs,totabs*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totabs=zero

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
   ! print information
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
 ABI_FREE(px)
 ABI_FREE(py)
 ABI_FREE(pz)
 ABI_FREE(inter2w)
 ABI_FREE(inter1w)
 ABI_FREE(intra2w)
 ABI_FREE(intra1w)
 ABI_FREE(intra1wS)
 ABI_FREE(delta)

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
!!  nband_sum=Number of bands included in the sum. Must be <= mband
!!  pmat(mband,mband,nkpt,3,nsppol) = momentum matrix elements in cartesian coordinates(complex)
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

subroutine linelop(icomp, itemp, nband_sum, cryst, ks_ebands, &
                   pmat,v1,v2,v3,nmesh,de,sc,brod,tol,fnam,do_antiresonant,ncid,comm)

!Arguments ------------------------------------
 integer, intent(in) :: icomp, itemp, nband_sum, ncid
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ks_ebands
 complex(dpc), intent(in) :: pmat(ks_ebands%mband, ks_ebands%mband, ks_ebands%nkpt, 3, ks_ebands%nsppol)
 integer, intent(in) :: v1, v2, v3
 integer, intent(in) :: nmesh
 integer, intent(in) :: comm
 real(dp), intent(in) :: de
 real(dp), intent(in) :: sc
 real(dp), intent(in) :: brod
 real(dp), intent(in) :: tol
 character(len=*), intent(in) :: fnam
 logical, intent(in) :: do_antiresonant

!Local variables -------------------------
 integer,parameter :: master = 0
 integer :: iw
 integer :: i,j,k,lx,ly,lz
 integer :: isp,isym,ik
 integer :: ist1,istl,istn,istm, mband
 real(dp) :: ha2ev
 real(dp) :: t1,t2,t3
 real(dp) :: ene,totre,totabs,totim
 real(dp) :: el,en,em
 real(dp) :: emin,emax,my_emin,my_emax
 real(dp) :: const_esu,const_au,au2esu
 real(dp) :: wmn,wnm,wln,wnl,wml,wlm
 complex(dpc) :: idel,w,zi
 character(len=fnlen) :: fnam1,fnam2,fnam3,fnam4,fnam5
! local allocatable arrays
 real(dp), allocatable :: s(:,:), sym(:,:,:)
 integer :: start4(4),count4(4)
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
 integer :: my_rank, nproc, ierr, my_k1, my_k2
 integer :: fout1,fout2,fout3,fout4,fout5
 character(500) :: msg

! *********************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

!calculate the constant
 zi=(0._dp,1._dp)
 idel=zi*brod
! Disable symmetries for now
 const_au=-2._dp/(cryst%ucvol*dble(cryst%nsym))
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

 ! If there exists inversion symmetry exit with a mesg.
 if (cryst%idx_spatial_inversion() /= 0) then
   write(std_out,*) '-----------------------------------------'
   write(std_out,*) '    the crystal has inversion symmetry   '
   write(std_out,*) '    the LEO susceptibility is zero       '
   write(std_out,*) '-----------------------------------------'
   ABI_ERROR("Aborting now")
 end if

 ! check polarisation
 if (v1.le.0.or.v2.le.0.or.v3.le.0.or.v1.gt.3.or.v2.gt.3.or.v3.gt.3) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in linelop:                        '
   write(std_out,*) '    the polarisation directions incorrect    '
   write(std_out,*) '    1=x,  2=y  and 3=z                       '
   write(std_out,*) '---------------------------------------------'
   ABI_ERROR("Aborting now")
 end if

 ! number of energy mesh points
 if (nmesh.le.0) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in linelop:                        '
   write(std_out,*) '    number of energy mesh points incorrect   '
   write(std_out,*) '    number has to be integer greater than 0  '
   write(std_out,*) '    nmesh*de = max energy for calculation    '
   write(std_out,*) '---------------------------------------------'
   ABI_ERROR("Aborting now")
 end if

 ! step in energy
 if (de.le.zero) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in linelop:                        '
   write(std_out,*) '    energy step is incorrect                 '
   write(std_out,*) '    number has to real greater than 0.0      '
   write(std_out,*) '    nmesh*de = max energy for calculation    '
   write(std_out,*) '---------------------------------------------'
   ABI_ERROR("Aborting now")
 end if

 ! broadening
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

 ! tolerance
 if (tol.gt.0.006) then
   write(std_out,*) '----------------------------------------'
   write(std_out,*) '    ATTENTION: tolerance is too high    '
   write(std_out,*) '    ideally should be less than 0.004   '
   write(std_out,*) '----------------------------------------'
 end if

 mband = ks_ebands%mband
 ABI_MALLOC(enk, (mband))
 ABI_MALLOC(delta, (mband, mband, 3))
 ABI_MALLOC(rmnbc,(mband,mband, 3, 3))
 ABI_MALLOC(roverw,(mband, mband, 3, 3))
 ABI_MALLOC(rmna, (mband, mband, 3))
 ABI_MALLOC(chi, (nmesh))
 ABI_MALLOC(eta, (nmesh))
 ABI_MALLOC(sigma, (nmesh))
 ABI_MALLOC(chi2, (nmesh))
 ABI_MALLOC(sym, (3, 3, 3))
 ABI_MALLOC(s, (3, 3))

 ABI_CHECK(nband_sum <= mband, "nband_sum <= mband")

 ! generate the symmetrizing tensor
 sym(:,:,:)=zero
 do isym=1,cryst%nsym
   s(:,:)=cryst%symrel_cart(:,:,isym)
   do i=1,3
     do j=1,3
       do k=1,3
         sym(i,j,k)=sym(i,j,k)+(s(i,v1)*s(j,v2)*s(k,v3))
       end do
     end do
   end do
 end do

 ! initialise
 delta(:,:,:)=zero
 rmnbc(:,:,:,:)=zero
 chi(:)=zero
 chi2(:) = zero
 eta(:)=zero
 sigma(:)=zero
 my_emin=HUGE(zero)
 my_emax=-HUGE(zero)

 ! Split work
 call xmpi_split_work(ks_ebands%nkpt,comm,my_k1,my_k2)

 ! loop over kpts
 do ik=my_k1,my_k2
   write(std_out,*) "P-",my_rank,": ",ik,'of',ks_ebands%nkpt
   do isp=1,ks_ebands%nsppol
     ! Calculate the scissor corrected energies and the energy window
     do ist1=1,nband_sum
       en = ks_ebands%eig(ist1,ik,isp)
       my_emin=min(my_emin,en)
       my_emax=max(my_emax,en)
       if(en > ks_ebands%fermie) then
         en = en + sc
       end if
       enk(ist1) = en
     end do

     ! calculate \Delta_nm and r_mn^a
     do istn=1,nband_sum
       en = enk(istn)
       do istm=1,nband_sum
         em = enk(istm)
         wmn = em - en
         delta(istn,istm,1:3)=pmat(istn,istn,ik,1:3,isp)-pmat(istm,istm,ik,1:3,isp)
         if(abs(wmn) < tol) then
           rmna(istm,istn,1:3) = zero
         else
           rmna(istm,istn,1:3)=-zi*pmat(istm,istn,ik,1:3,isp)/wmn
         end if
       end do
     end do

     ! calculate \r^b_mn;c
     do istm=1,nband_sum
       em = enk(istm)
       do istn=1,nband_sum
         en = enk(istn)
         wmn = em - en
         if(abs(wmn) > tol) then
           do ly = 1,3
             do lz = 1,3
               num1 = (rmna(istm,istn,ly)*delta(istm,istn,lz))+(rmna(istm,istn,lz)*delta(istm,istn,ly))
               den1 = wmn
               term1 = num1/den1
               term2 = zero
               do istp=1,nband_sum
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

     ! initialise the factors
     ! start the calculation
     do istn=1,nband_sum
       en=enk(istn)
       if (do_antiresonant .and. en .ge. ks_ebands%fermie) then
         cycle
       end if
       fn=ks_ebands%occ(istn,ik,isp)
       do istm=1,nband_sum
         em=enk(istm)
         if (do_antiresonant .and. em .le. ks_ebands%fermie) then
           cycle
         end if
         wmn=em-en
         wnm=-wmn
         fm = ks_ebands%occ(istm,ik,isp)
         fnm = fn - fm
         fmn = fm - fn
         eta1 = zero
         eta2_1 = zero
         eta2_2 = zero
         sigma1_1 = zero
         sigma1_2 = zero
         sigma2 = zero
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
         chi1_1 = zero
         chi1_2 = zero
         chi2_1b = zero
         chi2_2b = zero
         chi2(:) = zero
         ! Three band terms
         do istl=1,nband_sum
           el=enk(istl)
           fl = ks_ebands%occ(istl,ik,isp)
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

         !  calculate over the desired energy mesh and sum over k-points
         do iw=1,nmesh
           w=(iw-1)*de+idel
           ! Better way to compute it
           chi(iw) = chi(iw) + 0.5_dp*ks_ebands%wtk(ik)*((chi1/(wmn-w)) + ((chi2_1b+chi2_2b)/(wmn-w)))*const_esu
           eta(iw) = eta(iw) + 0.5_dp*zi*ks_ebands%wtk(ik)*((eta1/(wmn-w)) + (eta2/((wmn-w)**2)))*const_esu
           sigma(iw) = sigma(iw) + 0.5_dp*zi*ks_ebands%wtk(ik)*((sigma1/(wmn-w))- (sigma2/(wmn-w)))*const_esu
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
     ABI_ERROR(msg)
   end if
   if (open_file(fnam2,msg,newunit=fout2,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   if (open_file(fnam3,msg,newunit=fout3,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   if (open_file(fnam4,msg,newunit=fout4,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   if (open_file(fnam5,msg,newunit=fout5,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
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

   totim=zero
   totre=zero
   totabs=zero
   do iw=2,nmesh
     ene=(iw-1)*de
     ene=ene*ha2ev
     totim=aimag(chi(iw)+eta(iw)+sigma(iw))/1.d-7
     write(fout1,'(f15.6,2es15.6)') ene,totim,totim*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totim=zero
     totre=dble(chi(iw)+eta(iw)+sigma(iw))/1.d-7
     write(fout2,'(f15.6,2es15.6)') ene,totre,totre*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totre=zero
     write(fout3,'(f15.6,3es15.6)') ene,aimag(chi(iw))/1.d-7,      &
     aimag(eta(iw))/1.d-7,aimag(sigma(iw))/1.d-7
     write(fout4,'(f15.6,3es15.6)') ene,dble(chi(iw))/1.d-7,       &
     dble(eta(iw))/1.d-7,dble(sigma(iw))/1.d-7
     totabs=abs(chi(iw)+eta(iw)+sigma(iw))/1.d-7
     write(fout5,'(f15.6,2es15.6)') ene,totabs,totabs*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totabs=zero
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
!!  nband_sum=Number of bands included in the sum. Must be <= mband
!!  pmat(mband,mband,nkpt,3,nsppol) = momentum matrix elements in cartesian coordinates(complex)
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

subroutine nonlinopt(icomp, itemp, nband_sum, cryst, ks_ebands, &
                      pmat, v1, v2, v3, nmesh, de, sc, brod, tol, fnam, do_antiresonant, ncid, comm)

!Arguments ------------------------------------
integer, intent(in) :: icomp, itemp, nband_sum, ncid
type(crystal_t),intent(in) :: cryst
type(ebands_t),intent(in) :: ks_ebands
complex(dpc), intent(in) :: pmat(ks_ebands%mband, ks_ebands%mband, ks_ebands%nkpt, 3, ks_ebands%nsppol)
integer, intent(in) :: v1, v2, v3
integer, intent(in) :: nmesh
integer, intent(in) :: comm
real(dp), intent(in) :: de, sc, brod, tol
character(len=*), intent(in) :: fnam
logical, intent(in) :: do_antiresonant

!Local variables -------------------------
integer :: iw,i,j,k,lx,ly,lz,mband
integer :: isp,isym,ik,ist1,istl,istn,istm
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
 real(dp) :: ep, wmp, wpn, wtk
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
 integer,parameter :: master = 0
 integer :: ierr
 integer :: my_k1, my_k2
 character(500) :: msg
 integer :: fout1,fout2,fout3,fout4,fout5,fout6,fout7

! *********************************************************************

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

!calculate the constant
 zi=(0._dp,1._dp)
 idel=zi*brod
 const_au=-2._dp/(cryst%ucvol*dble(cryst%nsym))
 !const_au=-2._dp/(cryst%ucvol)
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

 ! If there exists inversion symmetry exit with a mesg.
 if (cryst%idx_spatial_inversion() /= 0) then
   write(std_out,*) '-----------------------------------------'
   write(std_out,*) '    the crystal has inversion symmetry   '
   write(std_out,*) '    the nl electro-optical susceptibility'
   write(std_out,*) '    is zero                              '
   write(std_out,*) '-----------------------------------------'
   ABI_ERROR("Aborting now")
 end if

 ! check polarisation
 if (v1.le.0.or.v2.le.0.or.v3.le.0.or.v1.gt.3.or.v2.gt.3.or.v3.gt.3) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in nonlinopt:                        '
   write(std_out,*) '    the polarisation directions incorrect    '
   write(std_out,*) '    1=x,  2=y  and 3=z                       '
   write(std_out,*) '---------------------------------------------'
   ABI_ERROR("Aborting now")
 end if

 ! number of energy mesh points
 if (nmesh.le.0) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in nonlinopt:                        '
   write(std_out,*) '    number of energy mesh points incorrect   '
   write(std_out,*) '    number has to be integer greater than 0  '
   write(std_out,*) '    nmesh*de = max energy for calculation    '
   write(std_out,*) '---------------------------------------------'
   ABI_ERROR("Aborting now")
 end if

 ! step in energy
 if (de.le.zero) then
   write(std_out,*) '---------------------------------------------'
   write(std_out,*) '    Error in nonlinopt:                        '
   write(std_out,*) '    energy step is incorrect                 '
   write(std_out,*) '    number has to real greater than 0.0      '
   write(std_out,*) '    nmesh*de = max energy for calculation    '
   write(std_out,*) '---------------------------------------------'
   ABI_ERROR("Aborting now")
 end if

 ! broadening
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

 ! tolerance
 if (tol.gt.0.006) then
   write(std_out,*) '----------------------------------------'
   write(std_out,*) '    ATTENTION: tolerance is too high    '
   write(std_out,*) '    ideally should be less than 0.004   '
   write(std_out,*) '----------------------------------------'
 end if

 ! allocate local arrays
 mband = ks_ebands%mband
 ABI_CHECK(nband_sum <= mband, "nband_sum <= mband")
 ABI_MALLOC(enk, (mband))
 ABI_MALLOC(delta, (mband, mband, 3))
 ABI_MALLOC(rmnbc, (mband, mband, 3, 3))
 ABI_MALLOC(roverw, (mband, mband, 3, 3))
 ABI_MALLOC(rmna, (mband, mband, 3))
 ABI_MALLOC(chiw, (nmesh))
 ABI_MALLOC(etaw, (nmesh))
 ABI_MALLOC(chi2w, (nmesh))
 ABI_MALLOC(eta2w, (nmesh))
 ABI_MALLOC(sigmaw, (nmesh))
 ABI_MALLOC(chi2, (nmesh))
 ABI_MALLOC(eta1, (nmesh))
 ABI_MALLOC(symrmn, (mband, mband, mband))

 ! generate the symmetrizing tensor
 sym = zero
 do isym=1,cryst%nsym
   s(:,:)=cryst%symrel_cart(:,:,isym)
   do i=1,3
     do j=1,3
       do k=1,3
         sym(i,j,k)=sym(i,j,k)+(s(i,v1)*s(j,v2)*s(k,v3))
       end do
     end do
   end do
 end do

 ! initialise
 delta = zero
 rmnbc = zero
 chiw = zero
 chi2w = zero
 chi2 = zero
 etaw = zero
 eta2w = zero
 sigmaw = zero
 my_emin=HUGE(one)
 my_emax=-HUGE(one)

 ! Split work
 call xmpi_split_work(ks_ebands%nkpt, comm, my_k1, my_k2)

! loop over kpts
 do ik=my_k1,my_k2
   write(std_out,*) "P-",my_rank,": ",ik,'of ', ks_ebands%nkpt
   do isp=1,ks_ebands%nsppol
     ! Calculate the scissor corrected energies and the energy window
     do ist1=1,nband_sum
       en = ks_ebands%eig(ist1,ik,isp)
       my_emin=min(my_emin,en)
       my_emax=max(my_emax,en)
       if(en > ks_ebands%fermie) then
         en = en + sc
       end if
       enk(ist1) = en
     end do

     ! calculate \Delta_nm and r_mn^a
     do istn=1,nband_sum
       en = enk(istn)
       do istm=1,nband_sum
         em = enk(istm)
         wmn = em - en
         delta(istn,istm,1:3)=pmat(istn,istn,ik,1:3,isp)-pmat(istm,istm,ik,1:3,isp)
         if(abs(wmn) < tol) then
           rmna(istm,istn,1:3) = zero
         else
           rmna(istm,istn,1:3)=pmat(istm,istn,ik,1:3,isp)/wmn
         end if
       end do
     end do

     ! calculate \r^b_mn;c
     do istm=1,nband_sum
       em = enk(istm)
       do istn=1,nband_sum
         en = enk(istn)
         wmn = em - en
         if (abs(wmn) < tol) then ! Degenerate energies
           rmnbc(istm,istn,:,:) = zero
           roverw(istm,istn,:,:) = zero
           cycle
         end if
         do ly = 1,3
           do lz = 1,3
             num1 = (rmna(istm,istn,ly)*delta(istm,istn,lz))+(rmna(istm,istn,lz)*delta(istm,istn,ly))
             den1 = wmn
             term1 = num1/den1
             term2 = zero
             do istp=1,nband_sum
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

     ! initialise the factors
     ! start the calculation
     do istn=1,nband_sum
       en=enk(istn)
       fn=ks_ebands%occ(istn,ik,isp)
       if(do_antiresonant .and. en .ge. ks_ebands%fermie) then
         cycle
       end if
       do istm=1,nband_sum
         em=enk(istm)
         if (do_antiresonant .and. em .le. ks_ebands%fermie) then
           cycle
         end if
         wmn=em-en
         wnm=-wmn
         fm = ks_ebands%occ(istm,ik,isp)
         fnm = fn - fm
         if(abs(wmn) > tol) then
           chi1 = zero
           chi2(:) = zero
           chi2_1 = zero
           chi2_2 = zero
           eta1(:) = zero
           eta1_1 = zero
           eta1_2 = zero
           eta2_1 = zero
           eta2_2 = zero
           sigma1 = zero
           sigma2_1 = zero
           ! Three band terms
           do istl=1,nband_sum
             el=enk(istl)
             fl = ks_ebands%occ(istl,ik,isp)
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
           eta2_1 = zero
           sigma2_1 = zero
           do lx = 1,3
             do ly = 1,3
               do lz = 1,3
                 eta2_1 = eta2_1 + sym(lx,ly,lz)*fnm*rmna(istn,istm,lx)*0.5_dp &
                    *(delta(istm,istn,ly)*rmna(istm,istn,lz)+delta(istm,istn,lz)*rmna(istm,istn,ly))
                 ! Correct version (Sipe 1993)
                 sigma2_1 = sigma2_1 + sym(lx,ly,lz)*fnm*rmna(istn,istm,lx)*0.5_dp &
                    *(rmna(istm,istn,ly)*delta(istn,istm,lz)+rmna(istm,istn,lz)*delta(istn,istm,ly))

                 ! Incorrect version (Hughes 1996)
                 !sigma2_1 = fnm*delta(istn,istm,v1)*0.5_dp*(rmna(istm,istn,v2)*rmna(istn,istm,v3)+rmna(istm,istn,v3)*rmna(istn,istm,v2))
               end do
             end do
           end do

           ! calculate over the desired energy mesh and sum over k-points
           wtk = ks_ebands%wtk(ik)
           do iw=1,nmesh
             w=(iw-1)*de+idel
             chi2w(iw) = chi2w(iw) + zi*wtk*((2.0_dp*fnm*chi1/(wmn-2.0_dp*w)))*const_esu ! Inter(2w) from chi
             chiw(iw) = chiw(iw) + zi*wtk*((chi2_1+chi2_2)/(wmn-w))*const_esu ! Inter(w) from chi
             eta2w(iw) = eta2w(iw) + zi*wtk*(8.0_dp*(eta2_1/((wmn**2)*(wmn-2.0_dp*w))) &
                 + 2.0_dp*eta2_2/((wmn**2)*(wmn-2.0_dp*w)))*const_esu ! Intra(2w) from eta
             etaw(iw) = etaw(iw) + zi*wtk*((eta1_1 + eta1_2)*fnm/((wmn**2)*(wmn-w)))*const_esu ! Intra(w) from eta
             sigmaw(iw) = sigmaw(iw) + 0.5_dp*zi*wtk*(fnm*sigma1/((wmn**2)*(wmn-w)) &
                 + (sigma2_1/((wmn**2)*(wmn-w))))*const_esu ! Intra(1w) from sigma
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
     ABI_ERROR(msg)
   end if
   if (open_file(fnam2,msg,newunit=fout2,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   if (open_file(fnam3,msg,newunit=fout3,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   if (open_file(fnam4,msg,newunit=fout4,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   if (open_file(fnam5,msg,newunit=fout5,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   if (open_file(fnam6,msg,newunit=fout6,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
   end if
   if (open_file(fnam7,msg,newunit=fout7,action='WRITE',form='FORMATTED') /= 0) then
     ABI_ERROR(msg)
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

   totim=zero
   totre=zero
   totabs=zero
   do iw=2,nmesh
     ene=(iw-1)*de
     ene=ene*ha2ev

     totim=aimag(chiw(iw)+chi2w(iw)+etaw(iw)+eta2w(iw)+sigmaw(iw))/1.d-7
     write(fout1,'(f15.6,2es15.6)') ene,totim,totim*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totim=zero

     totre=dble(chiw(iw)+chi2w(iw)+eta2w(iw)+etaw(iw)+sigmaw(iw))/1.d-7
     write(fout2,'(f15.6,2es15.6)') ene,totre,totre*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totre=zero

     write(fout3,'(f15.6,4es15.6)') ene,aimag(chi2w(iw))/1.d-7,aimag(chiw(iw))/1.d-7,     &
     aimag(eta2w(iw))/1.d-7,aimag(etaw(iw))/1.d-7+aimag(sigmaw(iw))/1.d-7

     write(fout4,'(f15.6,4es15.6)') ene,dble(chi2w(iw))/1.d-7,aimag(chiw(iw))/1.d-7,       &
     dble(eta2w(iw))/1.d-7,dble(etaw(iw))/1.d-7+dble(sigmaw(iw))/1.d-7

     totabs=abs(chiw(iw)+chi2w(iw)+etaw(iw)+eta2w(iw)+sigmaw(iw))/1.d-7
     write(fout5,'(f15.6,2es15.6)') ene,totabs,totabs*4._dp*pi*(1._dp/30000._dp)*(1._dp/1.d-5)
     totabs=zero

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

end module m_optic_tools
