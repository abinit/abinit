!!****m* ABINIT/defs_PSolver
!! NAME
!! defs_PSolver
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! for the bigDFT Poisson Solver
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_PSolver

 implicit none
!!***

#if defined HAVE_BIGDFT
 interface
!!****m* defs_PSolver/PSolver
!! NAME
!! PSolver
!!
!! FUNCTION
!!
!! PARENTS
!!      mklocl_realspace
!!
!! CHILDREN
!!
!! SOURCE
   subroutine PSolver(geocode,datacode,iproc,nproc,n01,n02,n03,xc,hx,hy,hz,&
                      rhopot,karray,pot_ion,eh,exc,vxc,offset,sumpion,nspin)
    use module_base
    use module_types
    use module_xc
    use yaml_output
    use Poisson_Solver, except_dp => dp, except_gp => gp, except_wp => wp
    implicit none
    character(len=1),intent(in) :: geocode,datacode
    logical,intent(in) :: sumpion
    integer, intent(in) :: iproc,nproc,n01,n02,n03,nspin
    type(xc_info),intent(in) :: xc
    real(gp),intent(in) :: hx,hy,hz
    real(dp),intent(in) :: offset
    real(dp), dimension(*),intent(in) :: karray
    real(gp), intent(out) :: eh,exc,vxc
    real(dp), dimension(*),intent(inout) :: rhopot
    real(wp), dimension(*),intent(inout) :: pot_ion
   end subroutine PSolver
 end interface
!!***
#endif

end module defs_PSolver
!!***

