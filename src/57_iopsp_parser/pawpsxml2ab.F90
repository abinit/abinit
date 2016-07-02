!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawpsxml2ab
!! NAME
!! pawpsxml2ab
!!
!! FUNCTION
!!  From a XML format pseudopotential file which has already been read in,
!!  convert to abinit internal datastructures.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2016 ABINIT group (FJ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! psxml  = pseudopotential data structure
!! option = 1 for inpsphead call
!!          2 for pspatm call
!! OUTPUT
!! pspheads data structure is filled
!!
!! PARENTS
!!      inpspheads
!!
!! CHILDREN
!!      pawpsp_read_header_xml,pawpsp_read_pawheader
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine pawpsxml2ab( psxml, pspheads,option )

 use defs_basis
 use defs_datatypes
 use m_profiling_abi
 use m_errors
 use m_pawpsp,only: pawpsp_read_header_xml,pawpsp_read_pawheader
 use m_pawxmlps, only : paw_setup_t

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawpsxml2ab'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(paw_setup_t),intent(in) :: psxml
 type(pspheader_type),intent(inout) :: pspheads !vz_i
 integer ,intent(in) ::option

!Local variables-------------------------------
!scalars
 integer :: ii,il,lloc,mesh_size_ae,mesh_size_proj
 real(dp) :: r2well
! character(len=100) :: xclibxc
! character(len=500) :: message
!arrays

! *********************************************************************

 if(option==1.or.option==2) then
   call pawpsp_read_header_xml(lloc,pspheads%lmax,pspheads%pspcod,&
&   pspheads%pspxc,psxml,r2well,pspheads%zionpsp,pspheads%znuclpsp)

   call pawpsp_read_pawheader(pspheads%pawheader%basis_size,&
&   pspheads%lmax,pspheads%pawheader%lmn_size,&
&   pspheads%pawheader%l_size,pspheads%pawheader%mesh_size,&
&   pspheads%pawheader%pawver,psxml,pspheads%pawheader%rpaw,&
&   pspheads%pawheader%rshp,pspheads%pawheader%shape_type)

   pspheads%nproj=0
   do il=0,pspheads%lmax
     do ii=1,pspheads%pawheader%basis_size
       if(psxml%valence_states%state(ii)%ll==il) pspheads%nproj(il)=pspheads%nproj(il)+1
     end do
   end do
   pspheads%nprojso=0
   pspheads%pspdat=27061961
   pspheads%pspso=0
   pspheads%xccc=1
   pspheads%title=psxml%atom%symbol
!
 end if

 if (option==2) then
   write(std_out,*) "paw_setup version= ",psxml%version
   write(std_out,*)"atom symbol= ",psxml%atom%symbol,"Z= ",psxml%atom%znucl,&
&   "core= ",psxml%atom%zion,"valence= ",psxml%atom%zval
   write(std_out,*) "xc_functional  xc_type=",psxml%xc_functional%functionaltype,&
&   "xc_name= ",psxml%xc_functional%name
   write(std_out,*)"generator  type=",psxml%generator%gen,"name= ",psxml%generator%name
   write(std_out,*)"PAW_radius rpaw=",psxml%rpaw
   write(std_out,*)"valence_states"
   do ii=1,psxml%valence_states%nval
     write(std_out,*)"state n=",psxml%valence_states%state(ii)%nn,&
&     "l= ",psxml%valence_states%state(ii)%ll,&
&     "f= ",psxml%valence_states%state(ii)%ff,&
&     "rc= ",psxml%valence_states%state(ii)%rc,&
&     "e= ",psxml%valence_states%state(ii)%ee,&
&     "id= ",psxml%valence_states%state(ii)%id
   end do
   do ii=1,psxml%ngrid
     write(std_out,*)"radial_grid  eq= ",psxml%radial_grid(ii)%eq,&
&     "a= ",psxml%radial_grid(ii)%aa,&
&     "n= ",psxml%radial_grid(ii)%nn,&
&     "d= ",psxml%radial_grid(ii)%dd,&
&     "b= ",psxml%radial_grid(ii)%bb,&
&     "istart= ",psxml%radial_grid(ii)%istart,&
&     "iend= ",psxml%radial_grid(ii)%iend,&
&     "id= ",psxml%radial_grid(ii)%id
   end do
   write(std_out,*)"shape_function  type= ",psxml%shape_function%gtype,&
&   "rc= ",psxml%shape_function%rc,&
&   "lamb",psxml%shape_function%lamb
   if(psxml%ae_core_density%tread) then
     write(std_out,*)"ae_core_density grid=  ",psxml%ae_core_density%grid
     write(std_out,*)psxml%ae_core_density%data
   end if
   if(psxml%pseudo_core_density%tread) then
     write(std_out,*)"pseudo_core_density grid= ",psxml%pseudo_core_density%grid
     write(std_out,*)psxml%pseudo_core_density%data
   end if
   if(psxml%pseudo_valence_density%tread) then
     write(std_out,*)"pseudo_valence_densit grid= ",psxml%pseudo_valence_density%grid
     write(std_out,*)psxml%pseudo_valence_density%data
   end if
   if(psxml%zero_potential%tread) then
     write(std_out,*)"zero_potential grid= ",psxml%zero_potential%grid
     write(std_out,*)psxml%zero_potential%data
   end if
   if(psxml%ae_core_kinetic_energy_density%tread) then
     write(std_out,*)"ae_core_kinetic_energy_density grid= ",psxml%ae_core_kinetic_energy_density%grid
     write(std_out,*)psxml%ae_core_kinetic_energy_density%data
   end if
   if(psxml%pseudo_core_kinetic_energy_density%tread) then
     write(std_out,*)"pseudo_core_kinetic_energy_density grid= ",psxml%pseudo_core_kinetic_energy_density%grid
     write(std_out,*)psxml%pseudo_core_kinetic_energy_density%data
   end if
   if(psxml%kresse_joubert_local_ionic_potential%tread) then
     write(std_out,*)"kresse_joubert_local_ionic_potential grid =",psxml%kresse_joubert_local_ionic_potential%grid
     write(std_out,*)psxml%kresse_joubert_local_ionic_potential%data
   end if
   if(psxml%blochl_local_ionic_potential%tread) then
     write(std_out,*)"blochl_local_ionic_potential grid= ",psxml%blochl_local_ionic_potential%grid
     write(std_out,*)psxml%blochl_local_ionic_potential%data
   end if
   do ii=1,psxml%ngrid
     if(trim(psxml%ae_partial_wave(1)%grid)==trim(psxml%radial_grid(ii)%id)) then
       mesh_size_ae=psxml%radial_grid(ii)%iend-psxml%radial_grid(ii)%istart+1
     end if
   end do
   do ii=1,psxml%ngrid
     if(trim(psxml%projector_function(1)%grid)==trim(psxml%radial_grid(ii)%id)) then
       mesh_size_proj=psxml%radial_grid(ii)%iend-psxml%radial_grid(ii)%istart+1
     end if
   end do
   do ii=1,psxml%valence_states%nval
     write(std_out,*)"ae_partial_wave state= ",psxml%ae_partial_wave(ii)%state,&
&     "grid= ",psxml%ae_partial_wave(ii)%grid
     write(std_out,*)psxml%ae_partial_wave(ii)%data(1:mesh_size_ae)
     write(std_out,*)"pseudo_partial_wave state= ",psxml%pseudo_partial_wave(ii)%state,&
&     "grid= ",psxml%pseudo_partial_wave(ii)%grid
     write(std_out,*)psxml%pseudo_partial_wave(ii)%data(1:mesh_size_ae)
     write(std_out,*)"projector_function state= ",psxml%projector_function(ii)%state,&
&     "grid= ",psxml%projector_function(ii)%grid
     write(std_out,*)psxml%projector_function(ii)%data(1:mesh_size_proj)
   end do
   write(std_out,*)"kinetic_energy_differences"
   write(std_out,*)psxml%kinetic_energy_differences%data
 end if


end subroutine pawpsxml2ab
!!***
