!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_xc
!! NAME
!! defs_xc
!!
!! FUNCTION
!! This module contains definitions for the XC routines in the Src_3xc
!! directory and for the routines who call the XC routines, especially when
!! optional arguments are defined.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
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


module defs_xc

 use m_profiling_abi

 implicit none
!!***

 interface

!!****f* defs_xc/drivexc
!! NAME
!! drivexc
!!
!! FUNCTION
!! routine interface
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      drivexc_main
!!
!! CHILDREN
!!
!! SOURCE
  subroutine drivexc(exc,ixc,npts,nspden,order,rho_updn,vxc,ndvxc,ngr2,nvxcdgr,   & !Mandatory arguments
&  dvxc,d2vxc,exexch,grho2_updn,vxcgr)    !Optional arguments
   use defs_basis
   integer,intent(in) :: ixc,npts,nspden,order
   integer,intent(in) :: ndvxc,ngr2,nvxcdgr
   integer,intent(in),optional ::exexch
   real(dp),intent(in) :: rho_updn(npts,nspden)
   real(dp),intent(in), optional :: grho2_updn(npts,ngr2)
   real(dp),intent(out) :: exc(npts),vxc(npts,nspden)
   real(dp),intent(out), optional :: d2vxc(npts),dvxc(npts,ndvxc),vxcgr(npts,nvxcdgr)
  end subroutine drivexc
!!***
 end interface

 interface
!!****f* defs_xc/rhohxc
!! NAME
!! rhohxc
!!
!! FUNCTION
!! routine interface
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      calc_vhxc_me,energy,hybrid_corr,m_kxc,nonlinear,nres2vres,odamix,prcref
!!      prcref_PMA,respfn,rhotov,scfcv,setvtr,xchybrid_ncpp_cc
!!
!! CHILDREN
!!
!! SOURCE
  subroutine rhohxc(dtset,enxc,gsqcut,izero,kxc,mpi_enreg,nfft,ngfft,&
&  nhat,nhatdim,nhatgr,nhatgrdim,nkxc,nspden,n3xccc,option,rhog,rhor,rprimd,&
&  strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d,k3xc)
   use defs_basis
   use defs_abitypes
!scalars
   integer,intent(in) :: izero,n3xccc,nfft,nhatdim,nhatgrdim,nkxc,nspden,option,usexcnhat
   real(dp),intent(in) :: gsqcut
   real(dp),intent(out) :: enxc,vxcavg
   type(MPI_type),intent(inout) :: mpi_enreg
   type(dataset_type),intent(in) :: dtset
!arrays
   integer,intent(in) :: ngfft(18)
   real(dp),intent(in) :: nhat(nfft,nspden*nhatdim),nhatgr(nfft,nspden,3*nhatgrdim)
   real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,nspden),rprimd(3,3)
   real(dp),intent(in) :: xccc3d(n3xccc)
   real(dp),intent(out) :: kxc(nfft,nkxc),strsxc(6),vhartr(nfft),vxc(nfft,nspden)
   real(dp),intent(out),optional :: k3xc(1:nfft)
  end subroutine rhohxc
!!***
 end interface

 interface
!!****f* defs_xc/xctetr
!! NAME
!! xctetr
!!
!! FUNCTION
!! routine interface
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!
!! SOURCE
  subroutine xctetr(exc,npt,order,rhor,rspts,vxc,& !Mandatory arguments
&                   d2vxc,dvxc)                    !Optional arguments
   use defs_basis
   integer,intent(in) :: npt,order
   real(dp),intent(in) :: rhor(npt),rspts(npt)
   real(dp),intent(out) :: exc(npt),vxc(npt)
   real(dp),intent(out), optional :: d2vxc(npt),dvxc(npt)
  end subroutine xctetr
!!***
 end interface

 interface
!!****f* defs_xc/xcpbe
!! NAME
!! xcpbe
!!
!! FUNCTION
!! routine interface
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!
!! SOURCE
  subroutine xcpbe(exci,npts,nspden,option,order,rho_updn,vxci,ndvxci,ngr2,& !Mandatory Arguments
&                  d2vxci,dvxcdgr,dvxci,exexch,grho2_updn)                     !Optional Arguments
   use defs_basis
   integer,intent(in) :: npts,nspden,option,order
   integer,intent(in) :: ndvxci,ngr2
   integer,intent(in),optional :: exexch
   real(dp),intent(in) :: rho_updn(npts,nspden)
   real(dp),intent(in), optional :: grho2_updn(npts,ngr2)
   real(dp),intent(out) :: exci(npts)
   real(dp),intent(out), optional :: d2vxci(npts),dvxcdgr(npts,3),dvxci(npts,ndvxci)
   real(dp),intent(out) :: vxci(npts,nspden)
  end subroutine xcpbe
!!***
 end interface

 interface
!!****f* defs_xc/xcspol
!! NAME
!! xcspol
!!
!! FUNCTION
!! routine interface
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!
!! SOURCE
  subroutine xcspol(exc,npts,nspden,order,rspts,vxc,zeta,ndvxc,& !Mandatory arguments
&                   dvxc)                                        !Optional arguments
   use defs_basis
   integer,intent(in) :: npts,nspden,order
   integer,intent(in) :: ndvxc
   real(dp),intent(in) :: rspts(npts),zeta(npts)
   real(dp),intent(out) :: exc(npts),vxc(npts,nspden)
   real(dp),intent(out), optional :: dvxc(npts,ndvxc)
  end subroutine xcspol
!!***
 end interface

 interface
!!****f* defs_xc/xcpzca
!! NAME
!! xcpzca
!!
!! FUNCTION
!! routine interface
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!
!! SOURCE
  subroutine xcpzca(exc,npt,order,rhor,rspts,vxc,&  !Mandatory arguments
&                   dvxc)                           !Optional arguments
   use defs_basis
   integer,intent(in) :: npt,order
   real(dp),intent(in) :: rhor(npt),rspts(npt)
   real(dp),intent(out) :: exc(npt),vxc(npt)
   real(dp),intent(out), optional :: dvxc(npt)
  end subroutine xcpzca
!!***
 end interface

 interface
!!****f* defs_xc/xcwign
!! NAME
!! xcwign
!!
!! FUNCTION
!! routine interface
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!
!! SOURCE
  subroutine xcwign(exc,npt,order,rspts,vxc,& !Mandatory arguments
&                   dvxc)                          !Optional arguments
   use defs_basis
   integer,intent(in) :: npt,order
   real(dp),intent(in) :: rspts(npt)
   real(dp),intent(out) :: exc(npt),vxc(npt)
   real(dp),intent(out), optional :: dvxc(npt)
  end subroutine xcwign
!!***
 end interface

 interface
!!****f* defs_xc/xchelu
!! NAME
!! xchelu
!!
!! FUNCTION
!! routine interface
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!
!! SOURCE
  subroutine xchelu(exc,npt,order,rspts,vxc,dvxc)  ! dvxc is optional
   use defs_basis
   integer,intent(in) :: npt,order
   real(dp),intent(in) :: rspts(npt)
   real(dp),intent(out) :: exc(npt),vxc(npt)
   real(dp),intent(out), optional :: dvxc(npt)
  end subroutine xchelu
!!***
 end interface

 interface
!!****f* defs_xc/xcxalp
!! NAME
!! xcxalp
!!
!! FUNCTION
!! routine interface
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      drivexc
!!
!! CHILDREN
!!
!! SOURCE
  subroutine xcxalp(exc,npt,order,rspts,vxc,dvxc)  ! dvxc is optional
   use defs_basis
   integer,intent(in) :: npt,order
   real(dp),intent(in) :: rspts(npt)
   real(dp),intent(out) :: exc(npt),vxc(npt)
   real(dp),intent(out),optional :: dvxc(npt)
  end subroutine xcxalp
!!***
 end interface


end module defs_xc
!!***
