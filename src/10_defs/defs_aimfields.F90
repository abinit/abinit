!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_aimfields
!! NAME
!! defs_aimfields
!!
!! FUNCTION
!! Default declarations and 2 spline subroutines for the aim.f utility.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (PCasek,FF,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! WARNING
!! This file does not follow the ABINIT coding rules (yet)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_aimfields

 use defs_basis

 implicit none

integer, save :: ngfft(3),nmax
integer, allocatable, save :: ndat(:)
real(dp), allocatable, target, save :: dig1(:),llg1(:),&
& dig2(:),llg2(:),dig3(:),llg3(:)
real(dp), allocatable, target, save :: cdig1(:),cdig2(:),cdig3(:)
real(dp), save :: dix(3)
real(dp), allocatable, target, save :: dvl(:,:,:),&
& ddx(:,:,:),ddy(:,:,:),ddz(:,:,:),rval(:,:,:)
real(dp), allocatable, save :: rrad(:,:),crho(:,:),sp2(:,:),sp3(:,:),&
& sp4(:,:),pdd(:),pd(:)
end module defs_aimfields
!!***
