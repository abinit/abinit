!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_bse_interp
!! NAME
!!  setup_bse_interp
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2009-2018 ABINIT group (Y. Gillet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! ngfft_gw(18)=Information about 3D FFT for density and potentials, see ~abinit/doc/variables/vargs.htm#ngfft
!! acell(3)=Length scales of primitive translations (bohr)
!! rprim(3,3)=Dimensionless real space primitive translations.
!! Dtset<dataset_type>=All input variables for this dataset.
!!  Some of them might be redefined here TODO
!! Dtfil=filenames and unit numbers used in abinit. fnameabi_wfkfile is changed is Fortran file is not 
!! found but a netcdf version with similar name is available.
!!
!! OUTPUT
!! Cryst<crystal_structure>=Info on the crystalline Structure.
!! Kmesh<BZ_mesh_type>=Structure defining the k-sampling for the wavefunctions.
!! Qmesh<BZ_mesh_type>=Structure defining the q-sampling for the symmetrized inverse dielectric matrix.
!! Gsph_x<gsphere_t=Data type gathering info on the G-sphere for wave functions and e^{-1},
!! KS_BSt<Bandstructure_type>=The KS band structure (energies, occupancies, k-weights...)
!! Vcp<vcoul_t>=Structure gathering information on the Coulomb interaction in reciprocal space,
!!   including a possible cutoff in real space.
!! ngfft_osc(18)=Contain all needed information about the 3D FFT for the oscillator matrix elements.
!!   See ~abinit/doc/variables/vargs.htm#ngfft
!! Bsp<excparam>=Basic parameters defining the Bethe-Salpeter run. Completely initialed in output.
!! Hdr_wfk<Hdr_type>=The header of the WFK file.
!! Hdr_bse<Hdr_type>=Local header initialized from the parameters used for the Bethe-Salpeter calculation.
!! w_file=File name used to construct W. Set to ABI_NOFILE if no external file is used.
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      apply_scissor,double_grid_init,ebands_copy,ebands_init,ebands_print
!!      ebands_report_gap,ebands_update_occ,find_qmesh,gsph_extend,gsph_init
!!      init_transitions,kmesh_init,kmesh_print,make_mesh,print_gsphere
!!      vcoul_init,wfk_read_eigenvalues,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine setup_bse_interp(Dtset,Dtfil,BSp,Cryst,Kmesh,&
& Kmesh_dense,Qmesh_dense,KS_BSt_dense,QP_bst_dense,Gsph_x,Gsph_c,Vcp_dense,Hdr_wfk_dense,ngfftf,grid,comm)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_bs_defs
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_nctk
 use m_hdr

 use m_gwdefs,        only : GW_Q0_DEFAULT
 use m_io_tools,      only : file_exists
 use m_crystal,       only : crystal_t
 use m_bz_mesh,       only : kmesh_t, kmesh_init, kmesh_print, find_qmesh, make_mesh
 use m_double_grid,   only : double_grid_t, double_grid_init
 use m_ebands,        only : ebands_init, ebands_print, ebands_copy, ebands_free, ebands_update_occ, &
                             apply_scissor, ebands_report_gap
 use m_vcoul,         only : vcoul_t, vcoul_init
 use m_gsphere,       only : gsphere_t, gsph_init, print_gsphere, gsph_extend
 use m_wfk,           only : wfk_read_eigenvalues

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'setup_bse_interp'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(dataset_type),intent(in) :: Dtset
 type(datafiles_type),intent(inout) :: Dtfil
 type(excparam),intent(inout) :: Bsp
 type(hdr_type),intent(out) :: Hdr_wfk_dense
 type(crystal_t),intent(in) :: Cryst
 type(kmesh_t),intent(in) :: Kmesh
 type(kmesh_t),intent(out) :: Kmesh_dense,Qmesh_dense
 type(ebands_t),intent(out) :: KS_BSt_dense,QP_Bst_dense
 type(double_grid_t),intent(out) :: grid
 type(vcoul_t),intent(out) :: Vcp_dense
 type(gsphere_t),intent(out) :: Gsph_x,Gsph_c
!arrays
 integer,intent(in) :: ngfftf(18)

!Local variables ------------------------------
!scalars
 integer,parameter :: pertcase0=0,master=0
 integer :: bantot_dense,ib,ibtot,ik_ibz,isppol,jj
 integer :: nbnds_kss_dense
 integer :: spin,hexc_size
 integer :: my_rank
 integer :: it
 integer :: nprocs
 integer :: is1,is2,is3,is4
 real(dp) :: nelect_hdr_dense
 logical,parameter :: remove_inv=.FALSE.
 character(len=500) :: msg
 character(len=fnlen) :: wfk_fname_dense
 integer :: nqlwl
!arrays
 integer,allocatable :: npwarr(:)
 real(dp),allocatable :: shiftk(:,:)
 real(dp),allocatable :: doccde(:),eigen(:),occfact(:)
 real(dp),pointer :: energies_p_dense(:,:,:)
 complex(dpc),allocatable :: gw_energy(:,:,:)
 integer,allocatable :: nbands_temp(:)
 integer :: kptrlatt_dense(3,3)
 real(dp),allocatable :: qlwl(:,:)
 real(dp) :: minmax_tene(2)

!************************************************************************

 DBG_ENTER("COLL")

 kptrlatt_dense = zero

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 SELECT CASE(BSp%interp_mode)
 CASE (1,2,3,4)
   nbnds_kss_dense = -1
   wfk_fname_dense = Dtfil%fnameabi_wfkfine
   call wrtout(std_out,"BSE Interpolation: will read energies from: "//trim(wfk_fname_dense),"COLL")

   if (nctk_try_fort_or_ncfile(wfk_fname_dense, msg) /= 0) then
     MSG_ERROR(msg)
   end if

   Dtfil%fnameabi_wfkfine = wfk_fname_dense

   call wfk_read_eigenvalues(wfk_fname_dense,energies_p_dense,Hdr_wfk_dense,comm)
   nbnds_kss_dense = MAXVAL(Hdr_wfk_dense%nband)
 CASE DEFAULT
   MSG_ERROR("Not yet implemented")
 END SELECT
 
 nelect_hdr_dense = Hdr_wfk_dense%nelect

 if (ABS(Dtset%nelect-nelect_hdr_dense)>tol6) then
   write(msg,'(2(a,f8.2))')&
&   "File contains ", nelect_hdr_dense," electrons but nelect initialized from input is ",Dtset%nelect
   MSG_ERROR(msg)
 end if

 ! Setup of the k-point list and symmetry tables in the  BZ 
 SELECT CASE(BSp%interp_mode)
 CASE (1,2,3,4)
   if(Dtset%chksymbreak == 0) then
     ABI_MALLOC(shiftk,(3,Dtset%nshiftk))
     kptrlatt_dense(:,1) = BSp%interp_kmult(1)*Dtset%kptrlatt(:,1)
     kptrlatt_dense(:,2) = BSp%interp_kmult(2)*Dtset%kptrlatt(:,2)
     kptrlatt_dense(:,3) = BSp%interp_kmult(3)*Dtset%kptrlatt(:,3)
     do jj = 1,Dtset%nshiftk
       shiftk(:,jj) = Bsp%interp_kmult(:)*Dtset%shiftk(:,jj)
     end do
     call make_mesh(Kmesh_dense,Cryst,Dtset%kptopt,kptrlatt_dense,Dtset%nshiftk,shiftk,break_symmetry=.TRUE.)
     ABI_FREE(shiftk)
   else
     !Initialize Kmesh with no wrapping inside ]-0.5;0.5]
     call kmesh_init(Kmesh_dense,Cryst,Hdr_wfk_dense%nkpt,Hdr_wfk_dense%kptns,Dtset%kptopt)
   end if
 CASE DEFAULT
   MSG_ERROR("Not yet implemented")
 END SELECT

 ! Init Qmesh
 call find_qmesh(Qmesh_dense,Cryst,Kmesh_dense)

 call gsph_init(Gsph_c,Cryst,0,ecut=Dtset%ecuteps)

 call double_grid_init(Kmesh,Kmesh_dense,Dtset%kptrlatt,BSp%interp_kmult,grid)

 BSp%nkibz_interp = Kmesh_dense%nibz  !We might allow for a smaller number of points....

 call kmesh_print(Kmesh_dense,"Interpolated K-mesh for the wavefunctions",std_out,Dtset%prtvol,"COLL")
 call kmesh_print(Kmesh_dense,"Interpolated K-mesh for the wavefunctions",ab_out, 0,           "COLL")

 if (nbnds_kss_dense < Dtset%nband(1)) then
   write(msg,'(2(a,i0),3a,i0)')&
&    'Interpolated WFK file contains only ', nbnds_kss_dense,' levels instead of ',Dtset%nband(1),' required;',ch10,&
&    'The calculation will be done with nbands= ',nbnds_kss_dense
   MSG_WARNING(msg)
   MSG_ERROR("Not supported yet !")
 end if

 ABI_MALLOC(nbands_temp,(Hdr_wfk_dense%nkpt*Hdr_wfk_dense%nsppol))
 do isppol=1,Hdr_wfk_dense%nsppol
   do ik_ibz=1,Hdr_wfk_dense%nkpt
     nbands_temp(ik_ibz+(isppol-1)*Hdr_wfk_dense%nkpt) = Dtset%nband(1)
   end do
 end do

 call gsph_extend(Gsph_c,Cryst,Dtset%ecutwfn,Gsph_x)
 call print_gsphere(Gsph_x,unit=std_out,prtvol=Dtset%prtvol)

 nqlwl=1
 ABI_MALLOC(qlwl,(3,nqlwl))
 qlwl(:,nqlwl)= GW_Q0_DEFAULT

 ! Compute Coulomb term on the largest G-sphere.
 if (Gsph_x%ng > Gsph_c%ng ) then
   call vcoul_init(Vcp_dense,Gsph_x,Cryst,Qmesh_dense,Kmesh_dense,Dtset%rcut,Dtset%icutcoul,&
&    Dtset%vcutgeo,Dtset%ecutsigx,Gsph_x%ng,nqlwl,qlwl,ngfftf,comm)
 else
   call vcoul_init(Vcp_dense,Gsph_c,Cryst,Qmesh_dense,Kmesh_dense,Dtset%rcut,Dtset%icutcoul,&
&    Dtset%vcutgeo,Dtset%ecutsigx,Gsph_c%ng,nqlwl,qlwl,ngfftf,comm)
 end if

 ABI_FREE(qlwl)

 bantot_dense=SUM(Hdr_wfk_dense%nband(1:Hdr_wfk_dense%nkpt*Hdr_wfk_dense%nsppol))
 ABI_ALLOCATE(doccde,(bantot_dense))
 ABI_ALLOCATE(eigen,(bantot_dense))
 ABI_ALLOCATE(occfact,(bantot_dense))
 doccde=zero; eigen=zero; occfact=zero

 jj=0; ibtot=0
 do isppol=1,Hdr_wfk_dense%nsppol
   do ik_ibz=1,Hdr_wfk_dense%nkpt
     do ib=1,Hdr_wfk_dense%nband(ik_ibz+(isppol-1)*Hdr_wfk_dense%nkpt)
       ibtot=ibtot+1
       if (ib<=BSP%nbnds) then
         jj=jj+1
         occfact(jj)=Hdr_wfk_dense%occ(ibtot)
         eigen  (jj)=energies_p_dense(ib,ik_ibz,isppol)
       end if
     end do
   end do
 end do

 ABI_FREE(energies_p_dense)

 ABI_MALLOC(npwarr,(kmesh_dense%nibz))
 npwarr=BSP%npwwfn

 call ebands_init(bantot_dense,KS_BSt_dense,Dtset%nelect,doccde,eigen,Hdr_wfk_dense%istwfk,Kmesh_dense%ibz,nbands_temp,&
&  Kmesh_dense%nibz,npwarr,Hdr_wfk_dense%nsppol,Hdr_wfk_dense%nspinor,Hdr_wfk_dense%tphysel,Hdr_wfk_dense%tsmear,&
&  Hdr_wfk_dense%occopt,occfact,Kmesh_dense%wt,&
&  hdr_wfk_dense%charge, hdr_wfk_dense%kptopt, hdr_wfk_dense%kptrlatt_orig, hdr_wfk_dense%nshiftk_orig, &
&  hdr_wfk_dense%shiftk_orig, hdr_wfk_dense%kptrlatt, hdr_wfk_dense%nshiftk, hdr_wfk_dense%shiftk)

 ABI_DEALLOCATE(doccde)
 ABI_DEALLOCATE(eigen)
 ABI_DEALLOCATE(npwarr)

 ABI_FREE(nbands_temp)

 ABI_FREE(occfact)

 !TODO Occupancies are zero if NSCF. One should calculate the occupancies from the energies when
 ! the occupation scheme for semiconductors is used.
 call ebands_update_occ(KS_BSt_dense,Dtset%spinmagntarget,prtvol=Dtset%prtvol)

 call ebands_print(KS_BSt_dense,"Interpolated band structure read from the WFK file",unit=std_out,prtvol=Dtset%prtvol)

 call ebands_report_gap(KS_BSt_dense,header="Interpolated KS band structure",unit=std_out,mode_paral="COLL")

 BSp%nkbz_interp = Kmesh_dense%nbz

 call ebands_copy(KS_BSt_dense,QP_bst_dense)

 SELECT CASE (Bsp%calc_type)
 CASE (BSE_HTYPE_RPA_KS)
   if (ABS(BSp%mbpt_sciss)>tol6) then
     write(msg,'(a,f8.2,a)')' Applying a scissors operator energy= ',BSp%mbpt_sciss*Ha_eV," [eV] on top of the KS energies."
     call wrtout(std_out,msg,"COLL")
     call apply_scissor(QP_BSt_dense,BSp%mbpt_sciss)
   else
     write(msg,'(a,f8.2,a)')' Using KS energies since mbpt_sciss= ',BSp%mbpt_sciss*Ha_eV," [eV]."
     call wrtout(std_out,msg,"COLL")
   end if
   !
 CASE (BSE_HTYPE_RPA_QPENE) ! Read _GW files with the corrections TODO here I should introduce variable getgw
   MSG_ERROR("Not yet implemented with interpolation !")
 CASE (BSE_HTYPE_RPA_QP)
   MSG_ERROR("Not implemented error!")
 CASE DEFAULT
   write(msg,'(a,i0)')"Unknown value for Bsp%calc_type= ",Bsp%calc_type
   MSG_ERROR(msg)
 END SELECT

 call ebands_report_gap(QP_BSt_dense,header=" Interpolated QP band structure",unit=std_out,mode_paral="COLL")

 ! Transitions are ALWAYS ordered in c-v-k mode with k being the slowest index.
 ! FIXME: linewidths not coded.
 ABI_ALLOCATE(gw_energy,(BSp%nbnds,Kmesh_dense%nibz,Dtset%nsppol))
 gw_energy = QP_BSt_dense%eig

 ABI_ALLOCATE(Bsp%nreh_interp,(Hdr_wfk_dense%nsppol))
 Bsp%nreh_interp=zero

 call init_transitions(BSp%Trans_interp,BSp%lomo_spin,BSp%humo_spin,BSp%ircut,Bsp%uvcut,BSp%nkbz_interp,Bsp%nbnds,&
&  Bsp%nkibz_interp,Hdr_wfk_dense%nsppol,Hdr_wfk_dense%nspinor,gw_energy,QP_BSt_dense%occ,Kmesh_dense%tab,minmax_tene,&
&  Bsp%nreh_interp)

 ABI_DEALLOCATE(gw_energy)

 do spin=1,Dtset%nsppol
   write(msg,'(a,i2,a,i0)')" For spin: ",spin,' the number of resonant e-h transitions is: ',BSp%nreh_interp(spin)
   call wrtout(std_out,msg,"COLL")
 end do

 if (ANY(Bsp%nreh_interp/=Bsp%nreh_interp(1))) then
   write(msg,'(a,(i0))')" BSE code does not support different number of transitions for the two spin channels",Bsp%nreh
   MSG_ERROR(msg)
 end if
 !
 ! Create transition table vcks2t
 is1=BSp%lomo_min;is2=BSp%homo_max;is3=BSp%lumo_min;is4=BSp%humo_max
 ABI_ALLOCATE(Bsp%vcks2t_interp,(is1:is2,is3:is4,BSp%nkbz_interp,Dtset%nsppol))
 Bsp%vcks2t_interp = 0

 do spin=1,Dtset%nsppol
   do it=1,BSp%nreh_interp(spin)
     BSp%vcks2t_interp(BSp%Trans_interp(it,spin)%v,BSp%Trans_interp(it,spin)%c,&
& BSp%Trans_interp(it,spin)%k,spin) = it
   end do
 end do

 hexc_size = SUM(Bsp%nreh_interp); if (Bsp%use_coupling>0) hexc_size=2*hexc_size
 if (Bsp%nstates_interp<=0) then
   Bsp%nstates_interp=hexc_size
 else
   if (Bsp%nstates_interp>hexc_size) then
      Bsp%nstates_interp=hexc_size
      write(msg,'(2(a,i0),2a)')&
&      "Since the total size of excitonic Hamiltonian ",hexc_size," is smaller than Bsp%nstates ",Bsp%nstates_interp,ch10,&
&      "the number of excitonic states nstates has been modified"
     MSG_WARNING(msg)
   end if
 end if

 DBG_EXIT("COLL")

end subroutine setup_bse_interp
!!***
