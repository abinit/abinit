!{\src2tex{textfont=tt}}
!!****f* ABINIT/initaim
!! NAME
!! initaim
!!
!! FUNCTION
!! Initialization for the 3D interpolation for the AIM code:
!!  - this procedure reads the charge density of the electrons of valence on
!!    the equidistant 3D grid (*_DEN output file of ABINIT) and the core charge
!!    density of electrons from *.fc files (fhi package)
!!  - the Cholesky decomposition  of the general matrix for
!!    the computation of the 1D spline coeficients in each direction is done.
!!    Warning - the procedure is modified to use periodic boundary conditions
!!    already during the decomposition
!!  - the second derivations of valence density in three directions are computed
!!    and stored in the real space grid of the density for interpolation.
!!  - the core density is stored separately in the radial grid together with th
!!    second radial derivation
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! aim_dtset= the structured entity containing all input variables
!!
!! OUTPUT
!! znucl_batom= the nuclear charge of the Bader atom
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  thie routine works on the data contained in the aim_fields and aim_prom modules
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!      drvaim
!!
!! CHILDREN
!!      bschg1,hdr_bcast,hdr_echo,hdr_fort_read,hdr_free,hdr_ncread,inspln
!!      lubksb,ludcmp,metric,xmpi_bcast
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine initaim(aim_dtset,znucl_batom)

 use defs_basis
 use defs_abitypes
 use defs_aimfields
 use defs_aimprom
 use defs_parameters
 use m_profiling_abi
 use m_errors
 use m_xmpi
#ifdef HAVE_TRIO_NETCDF
 use netcdf
#endif
 use m_hdr

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'initaim'
 use interfaces_28_numeric_noabirule
 use interfaces_41_geometry
 use interfaces_63_bader, except_this_one => initaim
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(aim_dataset_type),intent(in) :: aim_dtset

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: fform0,id,ierr,ii,info,jj,kk,kod,mm,ndtmax,nn,nsa,nsb,nsc,nsym,me,nproc,npsp
 integer :: unth,comm
#ifdef HAVE_TRIO_NETCDF
 integer :: den_id
#endif
 real(dp) :: ss,ucvol,znucl_batom
 type(hdr_type) :: hdr
!arrays
 integer :: ipiv(3)
 integer,allocatable :: symrel(:,:,:)
 real(dp) :: aa(3),bb(3),gmet(3,3),gprimd(3,3),rmet(3,3),yy(3,3)
 real(dp),allocatable :: tnons(:,:),znucl(:),zionpsp(:)
 real(dp),pointer :: ptc(:),ptd(:),ptf(:),ptp(:),ptsd(:)

! *********************************************************************

!DEBUG
!write(std_out,*) ' initaim : enter '
!ENDDEBUG

 comm = xmpi_world
 me=xmpi_comm_rank(comm)
 nproc=xmpi_comm_size(comm)

 slc=0    ! code for follow

!The use of the "hdr" routines is much better for the future
!maintenance of the code. Indeed, the content of the header
!will continue to change from time to time, and the associated
!changes should be done in one place only.

!Read ABINIT header ----------------------------------------------------------
 if(me==master)then
   if (aim_iomode == IO_MODE_ETSF) then
     call hdr_ncread(hdr, untad, fform0)
   else
     call hdr_fort_read(hdr, untad, fform0)
   end if
 end if
 ABI_CHECK(fform0 /= 0, "hdr_read returned fform == 0")
 call hdr_bcast(hdr,master,me,comm)

!Echo part of the header
 call hdr_echo(hdr, fform0, 4, unit=std_out)
 call hdr_echo(hdr, fform0, 4, unit=untout)

 natom=hdr%natom
 ngfft(1:3)=hdr%ngfft(:)
 nsym=hdr%nsym
 npsp=hdr%npsp
 ntypat=hdr%ntypat
 rprimd(:,:)=hdr%rprimd(:,:)

 ABI_ALLOCATE(zionpsp,(npsp))
 ABI_ALLOCATE(znucl,(ntypat))
 ABI_ALLOCATE(typat,(natom))
 ABI_ALLOCATE(xred,(3,natom))
 ABI_ALLOCATE(symrel,(3,3,nsym))
 ABI_ALLOCATE(tnons,(3,nsym))
 ABI_ALLOCATE(xatm,(3,natom))

 symrel(:,:,:)=hdr%symrel(:,:,:)
 typat(:)=hdr%typat(:)
 tnons(:,:)=hdr%tnons(:,:)
 znucl(:)=hdr%znucltypat(:)
 zionpsp(:)=hdr%zionpsp(:)
 xred(:,:)=hdr%xred(:,:)

!This is to deallocate records of hdr
 call hdr_free(hdr)

!-------------------------------------------------------------------------------

 ABI_ALLOCATE(dvl,(ngfft(1),ngfft(2),ngfft(3)))

 if(me==master)then
   if (aim_iomode == IO_MODE_ETSF) then
#ifdef HAVE_TRIO_NETCDF
     ! netcdf array has shape [cplex, n1, n2, n3, nspden]), here we read only the total density.
     NCF_CHECK(nf90_inq_varid(untad, "density", den_id))
     NCF_CHECK(nf90_get_var(untad, den_id, dvl, start=[1,1,1,1], count=[1, ngfft(1), ngfft(2), ngfft(3), 1]))
#endif
   else
     read(untad,iostat=nn) dvl(1:ngfft(1),1:ngfft(2),1:ngfft(3))
     ABI_CHECK(nn==0,"error of reading !")
   end if
 end if
 call xmpi_bcast(dvl, master, comm, ierr)

 write(std_out,*)ch10,' initaim : the valence density has been read' ,ch10

!INITIALISATION OF SOME IMPORTANT FIELDS

!Only interpolation is computed (inside vgh_rho) in reduced
!coordinates. In all other routines the cart. coordinates (CC) are used.

!transformation of the atom positions to CC
 do ii=1,natom
   xatm(:,ii)=xred(:,ii)
   call bschg1(xatm(:,ii),1)
 end do

!Generation of the neighbouring cells + transf to CC
 nn=0
 nsa=aim_dtset%nsa ; nsb=aim_dtset%nsb ; nsc=aim_dtset%nsc
 do ii=-nsa,nsa
   do jj=-nsb,nsb
     do kk=-nsc,nsc
       nn=nn+1
       atp(1,nn)=ii*1._dp
       atp(2,nn)=jj*1._dp
       atp(3,nn)=kk*1._dp
       call bschg1(atp(:,nn),1)
     end do
   end do
 end do
 nnpos=nn

!DEBUG
!write(std_out,*)' initaim : nnpos=',nnpos
!ENDDEBUG

 batcell=nsa*(2*nsb+1)*(2*nsc+1)+(2*nsc+1)*nsb+nsc+1
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 maxatdst=min(maxatdst, nsa*sqrt(rmet(1,1)), nsb*sqrt(rmet(2,2)), nsc*sqrt(rmet(3,3)) )
 if (maxcpdst > maxatdst) maxcpdst=0.75*maxatdst


!RPRIM ITS INVERSE AND TRANSPOSE

 do ii=1,3
   do jj=1,3
     yy(ii,jj)=rprimd(ii,jj)
   end do
 end do
 call ludcmp(yy,3,3,ipiv,id,info)
 ABI_CHECK(info==0,'Error inverting rprimd')

 do  ii=1,3
   do jj=1,3
     ivrprim(ii,jj)=0._dp
   end do
   ivrprim(ii,ii)=1._dp
 end do
 do ii=1,3
   call lubksb(yy,3,3,ipiv,ivrprim(:,ii))
 end do
 do ii=1,3
   do jj=1,3
     trivrp(ii,jj)=ivrprim(jj,ii)
   end do
 end do

 write(std_out,'(" INVERSE OF RPRIMD: ",/,3F16.8,/,3F16.8,/,3F16.8,/)') &
& ((ivrprim(ii,jj), jj=1,3), ii=1,3)
 write(untout,'(" INVERSE OF RPRIMD: ",/,3F16.8,/,3F16.8,/,3F16.8,/)') &
& ((ivrprim(ii,jj), jj=1,3), ii=1,3)

 write(std_out,*) "ATOMS (index,at.number,Zionic,position(xcart.))"
 write(std_out,*) "======================================="
 do ii=1,natom
   jj=typat(ii)
   write(std_out,'(I4,2F10.6,3F16.8)') ii, znucl(jj), zionpsp(jj), (xatm(kk,ii),kk=1,3)
 end do
 write(untout,*) "ATOMS (index,at.number,Zionic,position(xcart.))"
 write(untout,*) "======================================="
 do ii=1,natom
   jj=typat(ii)
   write(untout,'(I4,2F10.6,3F16.8)') ii, znucl(jj), zionpsp(jj), (xatm(kk,ii),kk=1,3)
 end do

!STEPS IN REAL SPACE GRID (REDUCED)
 do ii=1,3
   dix(ii)=1._dp/ngfft(ii)
 end do

!READING OF THE CORE DENSITY
 write(std_out,*)ch10,' initaim : will read the core densities' ,ch10

 ABI_ALLOCATE(ndat,(ntypat))
 ABI_ALLOCATE(rminl,(natom))
 ndtmax=0
 if(me==master)then
   do ii=1,ntypat
     unth=unt+ii
!    DEBUG
!    write(std_out,*)' read from unit ',unth
!    call flush(std_out)
!    stop
!    ENDDEBUG
     read(unth,*) ndat(ii),ss
     if (ndat(ii)>ndtmax) ndtmax=ndat(ii)
   end do
 end if
 call xmpi_bcast(ndat,master,comm,ierr)
 call xmpi_bcast(ndtmax,master,comm,ierr)
 call xmpi_bcast(ss,master,comm,ierr)

!FIELDS FOR STORING CORE DENSITY

 ABI_ALLOCATE(rrad,(ndtmax,ntypat))
 ABI_ALLOCATE(crho,(ndtmax,ntypat))
 ABI_ALLOCATE(sp2,(ndtmax,ntypat))
 ABI_ALLOCATE(sp3,(ndtmax,ntypat))
 ABI_ALLOCATE(sp4,(ndtmax,ntypat))
 ABI_ALLOCATE(corlim,(ntypat))

 sp2(:,:)=zero
 sp3(:,:)=zero
 sp4(:,:)=zero

!Reading of the core densities
 corlim(:)=0
 kod=0
 if(me==master)then
   do ii=1,ntypat
     unth=unt+ii
     do jj=1,ndat(ii)
       read(unth,*) rrad(jj,ii),crho(jj,ii),sp2(jj,ii),sp3(jj,ii)
       crho(jj,ii)=crho(jj,ii)/4._dp/pi
       if ((crho(jj,ii) < aim_rhocormin).and.(corlim(ii)==0)) corlim(ii)=jj
       sp2(jj,ii)=sp2(jj,ii)/4._dp/pi
       sp3(jj,ii)=sp3(jj,ii)/4._dp/pi   ! ATENTION!!! in sp3 is just second derivation
     end do
     do jj=1,ndat(ii)-1
       sp4(jj,ii)=(sp3(jj+1,ii)-sp3(jj,ii))/(6._dp*(rrad(jj+1,ii)-rrad(jj,ii)))
     end do
     if (corlim(ii)==0) corlim(ii)=ndat(ii)
   end do
 end if
 call xmpi_bcast(rrad,master,comm,ierr)
 call xmpi_bcast(crho,master,comm,ierr)
 call xmpi_bcast(sp2,master,comm,ierr)
 call xmpi_bcast(sp3,master,comm,ierr)
 call xmpi_bcast(sp4,master,comm,ierr)
 call xmpi_bcast(corlim,master,comm,ierr)

 write(std_out,*)ch10,' initaim : the core densities have been read' ,ch10

!CORRECTION OF THE CORE DENSITY NORMALISATION
 crho(:,:)=1.0003*crho(:,:)
 sp2(:,:)=1.0003*sp2(:,:)
 sp3(:,:)=1.0003*sp3(:,:)
 sp4(:,:)=1.0003*sp4(:,:)

!FIELDS FOR INTERPOLATIONS OF THE VALENCE DENSITY

 ABI_ALLOCATE(dig1,(ngfft(1)))
 ABI_ALLOCATE(dig2,(ngfft(2)))
 ABI_ALLOCATE(dig3,(ngfft(3)))
 ABI_ALLOCATE(llg1,(ngfft(1)))
 ABI_ALLOCATE(llg2,(ngfft(2)))
 ABI_ALLOCATE(llg3,(ngfft(3)))
 ABI_ALLOCATE(cdig1,(ngfft(1)-1))
 ABI_ALLOCATE(cdig2,(ngfft(2)-1))
 ABI_ALLOCATE(cdig3,(ngfft(3)-1))
 ABI_ALLOCATE(ddx,(ngfft(1),ngfft(2),ngfft(3)))
 ABI_ALLOCATE(ddy,(ngfft(1),ngfft(2),ngfft(3)))
 ABI_ALLOCATE(ddz,(ngfft(1),ngfft(2),ngfft(3)))

!DECOMPOSITION OF THE MATRIX FOR THE DETERMINATION OF COEFFICIENTS
!FOR CUBIC SPLINE INTERPOLATION (using the periodic boundary conditions)

!MAIN DIAGONAL (aa) AND SECONDARY DIAGONAL (bb) MATRIX ELEMENTS

 nmax=ngfft(1)
 do ii=2,3
   if (ngfft(ii) > nmax) nmax=ngfft(ii)
 end do
 nullify(ptf,ptsd)
 nullify(ptd,ptc,ptp)
 aa(:)=2.0*dix(:)**2/3.0
 bb(:)=dix(:)**2/6.0

 do ii=1,3
   if(ii==1) then
     ptd=>dig1;ptc=>cdig1;ptp=>llg1
   elseif (ii==2) then
     ptd=>dig2;ptc=>cdig2;ptp=>llg2
   else
     ptd=>dig3;ptc=>cdig3;ptp=>llg3
   end if
   ptd(1)=sqrt(aa(ii))
   ptc(1)=bb(ii)/ptd(1)
   ptp(1)=ptc(1)
   do jj=2,ngfft(ii)-1
     ptd(jj)=aa(ii)-ptc(jj-1)**2
     if(ptd(jj)<zero) then 
       MSG_ERROR('Matrix is not positive definite !')
     end if
     ptd(jj)=sqrt(ptd(jj))
     if (jj==ngfft(ii)-1) then
       ptc(jj)=(bb(ii)-ptp(jj-1)*ptc(jj-1))/ptd(jj)
       ptp(jj)=ptc(jj)
       exit
     end if
     ptc(jj)=bb(ii)/ptd(jj)
     ptp(jj)=-ptp(jj-1)*ptc(jj-1)/ptd(jj)
   end do
   ss=0._dp
   do jj=1,ngfft(ii)-1
     ss=ss+ptp(jj)**2
   end do
   ss=aa(ii)-ss
   if(ss<zero) then 
     MSG_ERROR('Matrix is not positive definite !')
   end if
   ptd(ngfft(ii))=sqrt(ss)
   ptp(ngfft(ii))=ptd(ngfft(ii))


!  INITIALISATION OF THE SECOND DERIVATIVE FIELDS

   nn=ii+1
   if (nn>3) nn=nn-3
   mm=ii+2
   if (mm>3) mm=mm-3
   do jj=1,ngfft(nn)
     do kk=1,ngfft(mm)
!      The calcul of the second derivations on the grid
       call inspln(ii,jj,kk)
     end do
   end do
   nullify(ptd,ptc,ptp)
 end do
 nullify(ptd,ptc,ptp)

 znucl_batom=znucl(typat(aim_dtset%batom))

 ABI_DEALLOCATE(znucl)
 ABI_DEALLOCATE(zionpsp)
 ABI_DEALLOCATE(symrel)
 ABI_DEALLOCATE(tnons)

!the pointers are obsolete - to remove later

end subroutine initaim
!!***
