!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_exc_spectra
!! NAME
!! m_exc_spectra
!!
!! FUNCTION
!!  Routines to compute the macroscopic dielectric function in the Bethe-Salpeter code.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2019 ABINIT and EXC groups (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida, M.Giantomassi, Y. Gillet)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_exc_spectra

 use defs_basis
 use defs_datatypes
 use m_bs_defs
 use m_abicore
 use iso_c_binding
 use m_xmpi
 use m_errors
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_ebands
 use m_hdr

 use m_io_tools,        only : open_file
 use m_fstrings,        only : toupper, strcat, sjoin, int2char4
 use m_numeric_tools,   only : simpson_int, simpson_cplx
 use m_hide_blas,       only : xdotu,xdotc
 use m_special_funcs,   only : gaussian
 use m_crystal,         only : crystal_t
 use m_bz_mesh,         only : kmesh_t
 use m_eprenorms,       only : eprenorms_t, renorm_bst
 use m_pawtab,          only : pawtab_type
 use m_paw_hr,          only : pawhur_t
 use m_wfd,             only : wfd_t
 !use m_bse_io,          only : exc_amplitude
 use m_wfd_optic,       only : calc_optical_mels

 implicit none

 private

 public :: build_spectra           ! Driver routine for the computation of optical spectra.
 public :: exc_write_data          ! This routine drives the writing of the files produced by the Bethe-Salpeter code.
 public :: exc_eps_rpa             ! Build epsilon within RPA and GW.
 public :: mdfs_ncwrite            ! Writes the MDF.nc file with the final results.
 public :: exc_write_tensor        ! Write of complex dielectric tensor
 !public :: exc_eps_resonant       ! Build the macroscopic dielectric function with excitonic effects.
!!***

contains

!!****f* m_exc_spectra/build_spectra
!! NAME
!!  build_spectra
!!
!! FUNCTION
!!  Driver routine for the computation of optical spectra.
!!
!! INPUTS
!!  usepaw=1 for PAW calculations, 0 otherwise.
!!  drude_plsmf=Drude plasma frequency.
!!  Bsp<excparam>=Data type gathering the paramenters used for the Bethe-Salpeter calculation.
!!    inclvkb=If different from 0, [Vnl,r] is included in the calculation of the matrix elements of the velocity operator.
!!  BS_files<excfiles>=filenames used in the Bethe-Salpeter part.
!!  Kmesh<kmesh_t>=the k-point sampling for the wave functions.
!!  Cryst<crystal_t>=Structure defining the crystalline structure.
!!  KS_BSt=The KS energies.
!!  QP_BSt=The QP energies.
!!  Psps <pseudopotential_type>=variables related to pseudopotentials.
!!  Pawtab(Cryst%ntypat*usepaw)<pawtab_type>=PAW tabulated starting data
!!  Hur(Cryst%natom*usepaw)<pawhur_t>=Only for PAW and LDA+U, quantities used to evaluate the commutator [H_u,r].
!!  Wfd<wfd_t>=Handler for the wavefunctions.
!!    nsppol=Number of independent spin polarizations.
!!    nspinor=Number of spinorial components.
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  No output. The routine calls specialized routines where the computation and the output of the spectra is done.
!!
!! PARENTS
!!      m_exc_diago
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine build_spectra(BSp,BS_files,Cryst,Kmesh,KS_BSt,QP_BSt,Psps,Pawtab,Wfd,Hur,drude_plsmf,comm,Epren)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 real(dp),intent(in) :: drude_plsmf
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 type(pseudopotential_type),intent(in) :: Psps
 type(kmesh_t),intent(in) :: Kmesh
 type(crystal_t),intent(in) :: Cryst
 type(ebands_t),intent(in) :: KS_BSt,QP_BSt
 type(wfd_t),intent(inout) :: Wfd
 type(eprenorms_t),optional,intent(in) :: Epren
!arrays
 type(pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*Wfd%usepaw)
 type(pawhur_t),intent(in) :: Hur(Cryst%natom*Wfd%usepaw)

!Local variables ------------------------------
!scalars
 integer :: my_rank,master,iq,io,nsppol,lomo_min,max_band,ncid
 integer :: itemp,ntemp
 logical :: do_ep_renorm
 real(dp) :: omegaev
 complex(dpc) :: ks_avg,gw_avg,exc_avg
 character(len=4) :: ts
 character(len=fnlen) :: path,prefix
 character(len=fnlen) :: filbseig, ost_fname
 !character(len=500) :: msg
 type(ebands_t) :: EPBSt, EP_QPBSt
!arrays
 real(dp),allocatable :: dos_exc(:),dos_gw(:),dos_ks(:)
 complex(dpc),allocatable :: eps_rpanlf(:,:),eps_gwnlf(:,:)
 complex(dpc),allocatable :: eps_exc(:,:),opt_cvk(:,:,:,:,:)

!************************************************************************

 my_rank = Wfd%my_rank
 master  = Wfd%master
 nsppol  = Wfd%nsppol

 do_ep_renorm = .False.
 ntemp = 1
 if (BSp%do_ep_renorm .and. PRESENT(Epren)) then
   do_ep_renorm = .True.
   ntemp = Epren%ntemp
 end if

 !
 ! =====================================================
 ! === Calculate fcv(k)=<c k s|e^{-iqr}|v k s> in BZ ===
 ! =====================================================
 lomo_min=Bsp%lomo_min; max_band=Bsp%nbnds
 ABI_MALLOC(opt_cvk,(lomo_min:max_band,lomo_min:max_band,BSp%nkbz,nsppol,BSp%nq))

 do iq=1,BSp%nq
   call calc_optical_mels(Wfd,Kmesh,KS_BSt,Cryst,Psps,Pawtab,Hur,BSp%inclvkb,Bsp%lomo_spin,lomo_min,max_band,&
&                         BSp%nkbz,BSp%q(:,iq),opt_cvk(:,:,:,:,iq))
 end do
 !
 ! ============================
 ! ==== Make EPS EXCITONIC ====
 ! ============================
 if (my_rank==master) then ! Only master works.

   ABI_MALLOC(eps_exc,(BSp%nomega,BSp%nq))
   ABI_MALLOC(dos_exc,(BSp%nomega))

   ABI_MALLOC(eps_rpanlf,(BSp%nomega,BSp%nq))
   ABI_MALLOC(dos_ks,(BSp%nomega))

   ABI_MALLOC(eps_gwnlf ,(BSp%nomega,BSp%nq))
   ABI_MALLOC(dos_gw,(BSp%nomega))

   do itemp = 1, ntemp

     call int2char4(itemp,ts)

     if(do_ep_renorm) then
       prefix = TRIM("_T") // ts
     else
       prefix = ""
     end if

     ost_fname = strcat(BS_files%out_basename,prefix,"_EXC_OST")

     !TODO for RPA
     call ebands_copy(KS_BST, EPBSt)
     call ebands_copy(QP_BST, EP_QPBSt)

     if (BS_files%in_eig /= BSE_NOFILE) then
       filbseig = strcat(BS_files%in_eig,prefix)
     else
       filbseig = strcat(BS_files%out_eig,prefix)
     end if

     if(do_ep_renorm) then
       ! No scissor with KSBST
       call renorm_bst(Epren, EPBSt, Cryst, itemp, do_lifetime=.TRUE.,do_check=.TRUE.)

       call renorm_bst(Epren, EP_QPBSt, Cryst, itemp, do_lifetime=.TRUE.,do_check=.FALSE.)
     end if


     if (BSp%use_coupling==0) then
       call exc_eps_resonant(BSp,filbseig,ost_fname,lomo_min,max_band,BSp%nkbz,nsppol,opt_cvk,&
&        Cryst%ucvol,BSp%nomega,BSp%omega,eps_exc,dos_exc,elph_lifetime=do_ep_renorm)
     else
       call exc_eps_coupling(Bsp,BS_files,lomo_min,max_band,BSp%nkbz,nsppol,opt_cvk,&
&        Cryst%ucvol,BSp%nomega,BSp%omega,eps_exc,dos_exc)
     end if
     !
     ! =======================================================
     ! === Make EPS RPA and GW without local-field effects ===
     ! =======================================================
     call wrtout(std_out," Calculating RPA NLF and QP NLF epsilon","COLL")

     call exc_eps_rpa(BSp%nbnds,BSp%lomo_spin,Bsp%lomo_min,BSp%homo_spin,Kmesh,EPBSt,BSp%nq,nsppol,opt_cvk,&
&      Cryst%ucvol,BSp%broad,BSp%nomega,BSp%omega,eps_rpanlf,dos_ks)

     call exc_eps_rpa(BSp%nbnds,BSp%lomo_spin,Bsp%lomo_min,BSp%homo_spin,Kmesh,EP_QPBSt,BSp%nq,nsppol,opt_cvk,&
&      Cryst%ucvol,Bsp%broad,BSp%nomega,BSp%omega,eps_gwnlf,dos_gw)
     !
     ! =========================
     ! === Write out Epsilon ===
     ! =========================
     !this is just for the automatic tests, It will be removed when fldiff
     !will be able to compare two optical spectral
     write(ab_out,*)" "
     write(ab_out,*)"Macroscopic dielectric function:"
     write(ab_out,*)"omega [eV] <KS_RPA_nlf>  <GW_RPA_nlf>  <BSE> "
     do io=1,MIN(10,BSp%nomega)
       omegaev = REAL(BSp%omega(io))*Ha_eV
       ks_avg  = SUM( eps_rpanlf(io,:)) / Bsp%nq
       gw_avg  = SUM( eps_gwnlf (io,:)) / Bsp%nq
       exc_avg = SUM( eps_exc   (io,:)) / Bsp%nq
       write(ab_out,'(7f9.4)')omegaev,ks_avg,gw_avg,exc_avg
     end do
     write(ab_out,*)" "

     !
     ! Master node writes final results on file.
     call exc_write_data(BSp,BS_files,"RPA_NLF_MDF",eps_rpanlf,prefix=prefix,dos=dos_ks)

     call exc_write_data(BSp,BS_files,"GW_NLF_MDF",eps_gwnlf,prefix=prefix,dos=dos_gw)

     call exc_write_data(BSp,BS_files,"EXC_MDF",eps_exc,prefix=prefix,dos=dos_exc)

     call wrtout(std_out," Checking Kramers Kronig on Excitonic Macroscopic Epsilon","COLL")
     call check_kramerskronig(BSp%nomega,REAL(BSp%omega),eps_exc(:,1))

     call wrtout(std_out," Checking Kramers Kronig on RPA NLF Macroscopic Epsilon","COLL")
     call check_kramerskronig(BSp%nomega,REAL(BSp%omega),eps_rpanlf(:,1))

     call wrtout(std_out," Checking Kramers Kronig on GW NLF Macroscopic Epsilon","COLL")
     call check_kramerskronig(BSp%nomega,REAL(BSp%omega),eps_gwnlf(:,1))

     call wrtout(std_out," Checking f-sum rule on Excitonic Macroscopic Epsilon","COLL")

     if (BSp%exchange_term>0) then
       MSG_COMMENT(' f-sum rule should be checked without LF')
     end if
     call check_fsumrule(BSp%nomega,REAL(BSp%omega),AIMAG(eps_exc(:,1)),drude_plsmf)

     call wrtout(std_out," Checking f-sum rule on RPA NLF Macroscopic Epsilon","COLL")
     call check_fsumrule(BSp%nomega,REAL(BSp%omega),AIMAG(eps_rpanlf(:,1)),drude_plsmf)

     call wrtout(std_out," Checking f-sum rule on GW NLF Macroscopic Epsilon","COLL")
     call check_fsumrule(BSp%nomega,REAL(BSp%omega),AIMAG(eps_gwnlf(:,1)),drude_plsmf)

#ifdef HAVE_NETCDF
     path = strcat(BS_files%out_basename, strcat(prefix,"_MDF.nc"))
     NCF_CHECK_MSG(nctk_open_create(ncid, path, xmpi_comm_self), sjoin("Creating MDF file:", path))
     NCF_CHECK(cryst%ncwrite(ncid))
     NCF_CHECK(ebands_ncwrite(QP_BSt, ncid))
     ! Write dielectric functions.
     call mdfs_ncwrite(ncid, Bsp, eps_exc,eps_rpanlf,eps_gwnlf)
     NCF_CHECK(nf90_close(ncid))
#else
     ABI_UNUSED(ncid)
#endif

     !TODO
     call ebands_free(EPBSt)
     call ebands_free(EP_QPBSt)

   end do

   ABI_FREE(eps_rpanlf)
   ABI_FREE(eps_gwnlf)
   ABI_FREE(eps_exc)
   ABI_FREE(dos_exc)
   ABI_FREE(dos_ks)
   ABI_FREE(dos_gw)
 end if ! my_rank==master

 ABI_FREE(opt_cvk)

 call xmpi_barrier(comm)

end subroutine build_spectra
!!***

!----------------------------------------------------------------------

!!****f* m_exc_spectra/exc_write_data
!! NAME
!!  exc_write_data
!!
!! FUNCTION
!!  This routine drives the writing of the files produced by the Bethe-Salpeter code.
!!
!! INPUTS
!! BSp<excparam>=Bethe-Salpeter Parameters.
!! what= "EXC_MDF"
!!       "RPA_NLF_MDF"
!!       "GW_NLF_MDF"
!! [dos(nomega)]
!!
!! OUTPUT
!!  Only writing.
!!
!! SIDE EFFECTS
!!  eps(BSp%nomega,BSp%nq) = Macroscopic dielectric function to be written.
!!
!! PARENTS
!!      m_exc_spectra,m_haydock
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine exc_write_data(BSp,BS_files,what,eps,prefix,dos)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: what
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
 character(len=*),optional,intent(in) :: prefix
!arrays
 real(dp),optional,intent(in) :: dos(BSp%nomega)
 complex(dpc),intent(in) :: eps(BSp%nomega,BSp%nq)

!Local variables ------------------------------
!scalars
 integer :: io,iq,funt
 real(dp) :: omegaev,step
 !real(dp),parameter :: SMALL=5.0d-99
!arrays
 real(dp) :: int_dos(BSp%nomega)
 real(dp) :: tmp_eps(2,BSp%nq)
 character(len=500) :: lf_type,block_type,wgg_type,frm,str_type,msg
 character(len=fnlen) :: fname

!************************************************************************

 if (PRESENT(prefix)) then
   fname = strcat(BS_files%out_basename,prefix,'_',toupper(what))
 else
   fname = strcat(BS_files%out_basename,'_',toupper(what))
 end if

 if (open_file(fname,msg,newunit=funt,form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 select case (toupper(what))
 case ("EXC_MDF")
   call wrtout(ab_out," Writing EXC Macroscopic dielectric function to file: "//trim(fname),"COLL")

   write(funt,'("# Macroscopic dielectric function obtained with the BS equation.")')

   lf_type = 'WITHOUT LOCAL FIELD EFFECTS'
   if (BSp%exchange_term>0) lf_type='LOCAL FIELD EFFECTS INCLUDED'
   call bsp_calctype2str(Bsp,str_type)
   write(funt,'("# ",a,"     " ,a)') TRIM(str_type), TRIM(lf_type)

   block_type = 'RESONANT-ONLY calculation'
   if (BSp%use_coupling>0) block_type = 'RESONANT+COUPLING calculation'
   write(funt,'("# ",a)') TRIM(block_type)

   if (BSp%use_coulomb_term) then
     wgg_type = "Coulomb term constructed with full W(G1,G2)"
     if ( BSp%use_diagonal_Wgg ) wgg_type = "Coulomb term constructed with diagonal approximation W(G1,G1)"
     write(funt,'("# ",a)') TRIM(wgg_type)
   end if

   write(funt,'(a,f7.4,a)')'# Scissor operator energy = ',BSp%mbpt_sciss*Ha_eV,' [eV]'

 case ("RPA_NLF_MDF")
   call wrtout(ab_out," Writing KS-RPA macroscopic dielectric function without local fields to file: "//trim(fname),"COLL")
   write(funt,'("# RPA macroscopic dielectric function without local fields")')

 case ("GW_NLF_MDF")
   call wrtout(ab_out," Writing GW-RPA macroscopic dielectric function without local fields to file: "//trim(fname),"COLL")

   write(funt,'("# GW Macroscopic dielectric function without local field effects ")')
   write(funt,'(a,f7.4,a)')'# Scissor operator energy = ',BSp%mbpt_sciss*Ha_eV,' [eV]'

 case default
   MSG_ERROR("Unknown value for what: "//trim(what))
 end select
 !
 ! Paramaters common to the different calculations.
 if (BSp%algorithm /= BSE_ALGO_HAYDOCK) then
   write(funt,'(a,i0)')"# nstates included in the diagonalization = ",BSp%nstates
 end if

 if (BSp%algorithm == BSE_ALGO_HAYDOCK) then
   write(funt,'(a,2f7.4)')'# Tolerance = ',BSp%haydock_tol
 end if

 write(funt,'(a,i0)')"# npweps  = ",BSp%npweps
 write(funt,'(a,i0)')"# npwwfn  = ",BSp%npwwfn
 write(funt,'(a,i0)')"# nbands  = ",BSp%nbnds
 write(funt,'(a,i0)')"# loband  = ",BSp%lomo_spin(1)
 if (Bsp%nsppol==2) write(funt,'(a,i0)')"# loband(spin=2) = ",BSp%lomo_spin(2)
 write(funt,'(a,i0)')"# nkibz   = ",BSp%nkibz
 write(funt,'(a,i0)')"# nkbz    = ",BSp%nkbz
 write(funt,'(a,f7.4,a)')'# Lorentzian broadening = ',BSp%broad*Ha_eV,' [eV]'
 !
 ! Write the list of q-points.
 write(funt,'(a)')"# List of q-points for the optical limit:"
 do iq=1,BSp%nq
   write(funt,'(a,3(f9.6,","),a)')'# q = ',BSp%q(:,iq),' [Reduced coords] '
 end do
 !
 ! Write spectra.
 if (.not.PRESENT(dos)) then
   write(funt,'(a)')"# omega [eV]    RE(eps(q=1)) IM(eps(q=1) RE(eps(q=2) ) ... "
   !write(frm,*)'(f7.3,',2*BSp%nq,'es12.4)'
   write(frm,*)'(f7.3,',2*BSp%nq,'(1x,f9.4))'
   do io=1,BSp%nomega
     omegaev = DBLE(BSp%omega(io))*Ha_eV
     tmp_eps(1,:) = REAL (eps(io,:))
     tmp_eps(2,:) = AIMAG(eps(io,:))
     !where (ABS(tmp_eps) < SMALL) ! this to improve the portability of the automatic tests.
     !  tmp_eps = zero
     !end where
     write(funt,frm) omegaev,(tmp_eps(:,iq), iq=1,BSp%nq)
   end do

 else
   write(funt,'(a)')"# omega [eV]    RE(eps(q=1)) IM(eps(q=1) RE(eps(q=2) ) ... DOS   IDOS"
   step = DBLE(BSp%omega(2) - BSp%omega(1))
   if ( ABS( step - DBLE((BSp%omega(BSp%nomega) - BSp%omega(BSp%nomega-1)))) > tol6 ) then
     MSG_WARNING("Frequency mesh must be linear for using simpson_int")
   end if
   call simpson_int(Bsp%nomega,step,dos,int_dos)
   !write(frm,*)'(f7.3,',2*BSp%nq,'es12.4,2es12.4)'
   write(frm,*)'(f7.3,',2*BSp%nq,'(1x,f9.4,1x,f9.4,1x,f9.4))'
   do io=1,BSp%nomega
     omegaev = DBLE(BSp%omega(io))*Ha_eV
     tmp_eps(1,:) = REAL (eps(io,:))
     tmp_eps(2,:) = AIMAG(eps(io,:))
     !where (ABS(tmp_eps) < SMALL) ! this to improve the portability of the automatic tests.
     !  tmp_eps = zero
     !end where
     !write(funt,frm) omegaev,(eps(io,iq), iq=1,BSp%nq), dos(io), int_dos(io)
     write(funt,frm) omegaev,(tmp_eps(:,iq), iq=1,BSp%nq), dos(io), int_dos(io)
   end do
 end if

 close(funt)

end subroutine exc_write_data
!!***

!----------------------------------------------------------------------

!!****f* m_exc_spectra/exc_eps_rpa
!! NAME
!!  exc_eps_rpa
!!
!! FUNCTION
!!  Build epsilon within RPA and GW.
!!
!! INPUTS
!! nkbz=Number of points in the BZ
!! nbnds=Number of bands
!! lomo_spin(nsppol)
!! lomo_min=Lowest occupied state
!! homo=Number of occupied states.
!! homo_spin(nsppol)
!! nsppol=Number of independent spin polarizations.
!! nomega=Number of frequencies
!! omega(nomega)=Frequency mesh.
!! ucvol=Unit cell volume.
!! broad=Broadening used for the DOS.
!! opt_cvk(nbnds,nbnds,nkbz)=Matrix elements <b k|e^{-iqr}|b" k> for a given q in the full BZ.
!!
!! OUTPUT
!!  eps_rpa(nomega)=RPA spectrum without local-field effects.
!!  dos(nomega)=The DOS.
!!
!! PARENTS
!!      m_exc_spectra,m_haydock
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine exc_eps_rpa(nbnds,lomo_spin,lomo_min,homo_spin,Kmesh,Bst,nq,nsppol,opt_cvk,ucvol,broad,nomega,omega,eps_rpa,dos)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbnds,lomo_min,nsppol,nomega,nq
 real(dp),intent(in) :: ucvol,broad
 type(kmesh_t),intent(in) :: Kmesh
 type(ebands_t),intent(in) :: BSt
!arrays
 integer,intent(in) :: lomo_spin(nsppol),homo_spin(nsppol)
 real(dp),intent(out) :: dos(nomega)
 complex(dpc),intent(in) :: omega(nomega)
 complex(dpc),intent(in) :: opt_cvk(lomo_min:nbnds,lomo_min:nbnds,Kmesh%nbz,nsppol,nq)
 complex(dpc),intent(out) :: eps_rpa(nomega,nq)

!Local variables ------------------------------
!scalars
 integer :: iw,ib_v,ib_c,ik_bz,ik_ibz,spin,iq
 real(dp) :: fact,arg,ediff
 real(dp) :: linewidth
 complex(dpc) :: ctemp
 logical :: do_linewidth

!************************************************************************

 ! TODO: four_pi comes from the bare Coulomb term hence the
 ! present implementation is not compatible with the cutoff technique.
 fact=four_pi/(ucvol*Kmesh%nbz)
 if (nsppol==1) fact=two*fact ! two accounts for the occupation factors.

 eps_rpa=czero; dos=zero

 do_linewidth = .FALSE.
 do_linewidth = allocated(BSt%linewidth)

 !write(std_out,*)nsppol,Kmesh%nbz,lomo_min,homo,nbnds
 !
 ! Sum over all QP transitions.
 do spin=1,nsppol
   do ik_bz=1,Kmesh%nbz
     ik_ibz = Kmesh%tab(ik_bz)
     do ib_v=lomo_spin(spin),homo_spin(spin)
       do ib_c=homo_spin(spin)+1,nbnds
         !
         ! TODO here energies are always assumed to be real.
         ediff = BSt%eig(ib_c,ik_ibz,spin) - BSt%eig(ib_v,ik_ibz,spin)

         !
         if(do_linewidth) then
           linewidth = BSt%linewidth(1,ib_c,ik_ibz,spin) + BSt%linewidth(1,ib_v,ik_ibz,spin)
           do iq=1,nq
             ctemp = opt_cvk(ib_c,ib_v,ik_bz,spin,iq)
             do iw=1,nomega
               eps_rpa(iw,iq) = eps_rpa(iw,iq)  + ctemp * CONJG(ctemp) *&
&             (one/(ediff-j_dpc*linewidth-omega(iw)) + one/(ediff+j_dpc*linewidth+omega(iw)))
             end do
           end do
           !
           ! The JDOS at q=0
           !if (ediff*Ha_eV < 0.3) then
           !  write(std_out,*)"Small transition ",ik_ibz,ib_v,ib_c
           !end if

           do iw=1,nomega
             arg = DBLE(omega(iw)) - ediff
             dos(iw) = dos(iw) + gaussian(arg, linewidth)
           end do
         else
           do iq=1,nq
             ctemp = opt_cvk(ib_c,ib_v,ik_bz,spin,iq)
             do iw=1,nomega
               eps_rpa(iw,iq) = eps_rpa(iw,iq)  + ctemp * CONJG(ctemp) *&
               (one/(ediff-omega(iw)) + one/(ediff+omega(iw)))
             end do
           end do
           !
           ! The JDOS at q=0
           !if (ediff*Ha_eV < 0.3) then
           !  write(std_out,*)"Small transition ",ik_ibz,ib_v,ib_c
           !end if

           do iw=1,nomega
             arg = DBLE(omega(iw)) - ediff
             dos(iw) = dos(iw) + gaussian(arg, broad)
           end do
         end if
         !
       end do !ib_c
     end do !ib_v
   end do !ik_bz
 end do !spin

 dos = dos/Kmesh%nbz
 eps_rpa = cone + fact*eps_rpa

end subroutine exc_eps_rpa
!!***

!----------------------------------------------------------------------

!!****f* m_exc_spectra/exc_eps_resonant
!! NAME
!!  exc_eps_resonant
!!
!! FUNCTION
!!  This routine builds the macroscopic dielectric function with excitonic effects.
!!
!! INPUTS
!! Bsp
!! lomo_min,max_band
!! nkbz=Number of points in the BZ
!! nsppol=Number of independent polarizations.
!! nomega=Number of frequencies
!! omega(nomega)=frequency mesh (complex shift is already included)
!! ucvol=Volume of the unit cell.
!! opt_cvk(lomo_min:max_band,mib:max_band,nkbz,nsppol,Bsp%nq)=Matrix elements <b k|e^{-iqr}|b" k> for a given q in the full BZ.
!!
!! OUTPUT
!!  eps_exc(nomega,Bsp%nq)=Macroscopic dielectric function with excitonic effects.
!!  dos_exc(nomega)=The DOS of the excitonic Hamiltonian
!!
!! PARENTS
!!      m_exc_spectra
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine exc_eps_resonant(Bsp,filbseig,ost_fname,lomo_min,max_band,nkbz,nsppol,opt_cvk,&
&    ucvol,nomega,omega,eps_exc,dos_exc,elph_lifetime)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lomo_min,max_band,nkbz,nomega,nsppol
 real(dp),intent(in) :: ucvol
 type(excparam),intent(in) :: BSp
 character(len=fnlen),intent(in) :: filbseig,ost_fname
 logical,optional,intent(in) :: elph_lifetime
!arrays
 real(dp),intent(out) :: dos_exc(nomega)
 complex(dpc),intent(in) :: opt_cvk(lomo_min:max_band,lomo_min:max_band,nkbz,nsppol,BSp%nq),omega(nomega)
 complex(dpc),intent(out) :: eps_exc(nomega,BSp%nq)

!Local variables ------------------------------
!scalars
 integer :: ll,it,iw,ib_v,ib_c,ik_bz,neig_read,eig_unt,exc_size,iq !fform,
 integer :: spin,spad,hsize_read,nstates,ost_unt
 logical :: do_ep_lifetime, file_do_lifetime
 real(dp) :: fact,arg
 complex(dpc) :: dotprod
 character(len=500) :: msg,frm,errmsg
 !type(Hdr_type) :: tmp_Hdr
!arrays
 real(dp),allocatable :: exc_ene(:)
 complex(dpc) :: ctemp(BSp%nq),dtemp(BSp%nq)
 complex(dpc),allocatable :: ostrength(:,:),exc_ene_cplx(:),exc_state(:),exc_state2(:)

!************************************************************************

 call wrtout(std_out," Calculating excitonic epsilon with antiresonant","COLL")

 if (nsppol==2) then
   MSG_WARNING("nsppol==2 still under development")
 end if

 exc_size = SUM(BSp%nreh)
 nstates  = BSp%nstates

 do_ep_lifetime = .FALSE.
 if(PRESENT(elph_lifetime)) then
   do_ep_lifetime = elph_lifetime
 end if

 if (ANY(Bsp%nreh/=Bsp%nreh(1))) then
   write(msg,'(a,2(i0,1x))')"BSE does not support different number of transitions for the two spin channels. nreh: ",Bsp%nreh
   MSG_WARNING(msg)
 end if
 !
 ! TODO:
 ! four_pi comes from the bare Coulomb term hence the
 ! present implementation is not compatible with the cutoff technique.
 fact=four_pi/(ucvol*nkbz); if (nsppol==1) fact=two*fact ! two to account for the occupation numbers.

 call wrtout(std_out," Reading excitonic eigenstates from file: "//TRIM(filbseig),"COLL")
 if (open_file(filbseig,msg,newunit=eig_unt,form="unformatted",status="old",action="read") /= 0) then
   MSG_ERROR(msg)
 end if

 read(eig_unt, err=10, iomsg=errmsg) file_do_lifetime

 if(do_ep_lifetime .and. .not. file_do_lifetime) then
  MSG_ERROR("Cannot do lifetime as the data is not present in the file !")
 end if

 read(eig_unt, err=10, iomsg=errmsg) hsize_read,neig_read

 if (hsize_read /= exc_size) then
   write(msg,'(2(a,i0))')" Wrong size of the Hamiltonian: read: ",hsize_read," expected= ",exc_size
   MSG_ERROR(msg)
 end if

 if (neig_read /= nstates) then
   write(msg,'(2(a,i0))')" Wrong number of eigenstates: read: ",neig_read," expected= ",nstates
   MSG_ERROR(msg)
 end if
 !
 ! Read eigenvalues, ignore possibly small imaginary part.
 ABI_MALLOC(exc_ene_cplx,(neig_read))
 read(eig_unt, err=10, iomsg=errmsg) exc_ene_cplx

 ABI_MALLOC(exc_ene,(neig_read))
 exc_ene = DBLE(exc_ene_cplx)
 !ABI_FREE(exc_ene_cplx)
 !
 ! Calculate oscillator strength.
 ABI_MALLOC(exc_state,(exc_size))
 ABI_MALLOC(ostrength,(neig_read,BSp%nq))

 if (do_ep_lifetime) then
   ABI_MALLOC(exc_state2,(exc_size))
   do ll=1,neig_read ! Loop over excitonic eigenstates reported on file.
     read(eig_unt, err=10, iomsg=errmsg) exc_state(:) ! Righteigenvector
     read(eig_unt, err=10, iomsg=errmsg) exc_state2(:) ! Lefteigenvector

     ! Here assuming that eigenvectors are such as Xl_i' Xr_j = delta_ij
     ! Otherwise, I need to invert the overlap matrix !

     ! Rescale the vectors so that they are "normalized" with respect to the other one !
     dotprod = xdotc(exc_size,exc_state2(:),1,exc_state(:),1)
     exc_state2(:) = exc_state2(:)/CONJG(dotprod)

     ctemp(:) = czero
     dtemp(:) = czero
     do spin=1,nsppol
       spad=(spin-1)*BSp%nreh(1) ! Loop over spin channels.
       do it=1,BSp%nreh(spin)    ! Loop over resonant transition t = (k,v,c,s)
         ik_bz = Bsp%Trans(it,spin)%k
         ib_v  = Bsp%Trans(it,spin)%v
         ib_c  = Bsp%Trans(it,spin)%c
         do iq=1,BSp%nq
           ctemp(iq) = ctemp(iq) + CONJG(opt_cvk(ib_c,ib_v,ik_bz,spin,iq)) * exc_state(it+spad)
           dtemp(iq) = dtemp(iq) + CONJG(opt_cvk(ib_c,ib_v,ik_bz,spin,iq)) * exc_state2(it+spad)
         end do
       end do ! it
     end do
     ostrength(ll,:) = ctemp(:)*CONJG(dtemp(:))
   end do ! ll
   ABI_FREE(exc_state2)
 else
   do ll=1,neig_read ! Loop over excitonic eigenstates reported on file.
     read(eig_unt, err=10, iomsg=errmsg) exc_state(:)
     if(file_do_lifetime) read(eig_unt, err=10, iomsg=errmsg)

     ctemp(:) = czero
     do spin=1,nsppol
       spad=(spin-1)*BSp%nreh(1) ! Loop over spin channels.
       do it=1,BSp%nreh(spin)    ! Loop over resonant transition t = (k,v,c,s)
         ik_bz = Bsp%Trans(it,spin)%k
         ib_v  = Bsp%Trans(it,spin)%v
         ib_c  = Bsp%Trans(it,spin)%c
         do iq=1,BSp%nq
           ctemp(iq) = ctemp(iq) + CONJG(opt_cvk(ib_c,ib_v,ik_bz,spin,iq)) * exc_state(it+spad)
         end do
       end do ! it
     end do
     ostrength(ll,:) = ctemp(:)*CONJG(ctemp(:))
   end do ! ll
 end if


 close(eig_unt, err=10, iomsg=errmsg)
 ABI_FREE(exc_state)

 if(do_ep_lifetime) then
   eps_exc = one
   do ll=1,neig_read ! Sum over all excitonic eigenstates read from file.
     do iq=1,BSp%nq
        do iw=1,nomega
          eps_exc(iw,iq) = eps_exc(iw,iq) +  &
  &         fact * ostrength(ll,iq) * (one/(exc_ene_cplx(ll) - omega(iw)) - one/(-DCONJG(exc_ene_cplx(ll)) - omega(iw)))
        end do
     end do !ll
   end do !iw
 else
   eps_exc = one
   do ll=1,neig_read ! Sum over all excitonic eigenstates read from file.
     do iq=1,BSp%nq
        do iw=1,nomega
          eps_exc(iw,iq) = eps_exc(iw,iq) +  &
  &         fact * ostrength(ll,iq) * (one/(exc_ene(ll) - omega(iw)) - one/(-exc_ene(ll) - omega(iw)))
        end do
     end do !ll
   end do !iw
 end if

 !
 ! The excitonic DOS.
 dos_exc=zero
 do ll=1,neig_read ! Sum over the calculate excitonic eigenstates.
   do iw=1,nomega
     arg = ( DBLE(omega(iw)) - exc_ene(ll))
     if(do_ep_lifetime) then
       dos_exc(iw) = dos_exc(iw) + gaussian(arg, AIMAG(exc_ene_cplx(ll)))
     else
       dos_exc(iw) = dos_exc(iw) + gaussian(arg, Bsp%broad)
     end if
   end do
 end do
 !
 ! Write the oscillator strengths to file.
 if (open_file(ost_fname,msg,newunit=ost_unt,form="formatted",action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 write(ost_unt,'("# Oscillator strengths of the excitonic states for the different q-polarizations.")')
 !
 ! Write the list of q-points.
 write(ost_unt,*)"# List of q-points for the optical limit"
 do iq=1,BSp%nq
   write(ost_unt,'(a,3(f9.6,","),a)')'# q = ',BSp%q(:,iq),' [Reduced coords] '
 end do

 write(ost_unt,*)"# E_lambda [eV]     ostrength(q=1) ostrength(q=2) .... "
 write(frm,*)'(f8.4,',BSp%nq,'es12.4)'
 do ll=1,neig_read
   write(ost_unt,frm)exc_ene(ll)*Ha_eV,(ostrength(ll,iq), iq=1,BSp%nq)
 end do

 close(ost_unt)

 ABI_FREE(ostrength)
 ABI_FREE(exc_ene)
 ABI_FREE(exc_ene_cplx)

 !call exc_amplitude(Bsp,filbseig,1,(/(ll,ll=1,10)/),"TEST_AMPLITUDE")
 !call exc_amplitude(Bsp,filbseig,1,(/30/),"TEST_AMPLITUDE")

 return

 ! Handler IO-error
10 continue
 MSG_ERROR(errmsg)

end subroutine exc_eps_resonant
!!***

!----------------------------------------------------------------------

!!****f* m_exc_spectra/exc_eps_coupling
!! NAME
!!  exc_eps_coupling
!!
!! FUNCTION
!!  Make epsilon EXCITONIC with full COUPLING.
!!
!! INPUTS
!! Bsp
!! nkbz=Number of points in the BZ
!! lomo_min,max_band
!! nomega=Number of frequencies
!! omega(nomega)=frequency mesh.
!! nsppol=Number of independent spin polarizations.
!! ucvol=Unit cell volume.
!! BS_files<excfiles>File names used in the Bethe-Salpeter code.
!! opt_cvk(lomo_min:max_band,lomo_min:max_band,nkbz,nsppol)=Matrix elements <b k|e^{-iqr}|b" k> for a given q in the full BZ.
!!
!! OUTPUT
!!  eps_exc(nomega)=Macroscopic dielectric function with excitonic effects calculated including the COUPLING.
!!  dos_exc(nomega)=The DOS of the excitonic Hamiltonian
!!
!! PARENTS
!!      m_exc_spectra
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine exc_eps_coupling(Bsp,BS_files,lomo_min,max_band,nkbz,nsppol,opt_cvk,ucvol,nomega,omega,eps_exc,dos_exc)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lomo_min,max_band,nkbz,nomega,nsppol
 real(dp),intent(in) :: ucvol
 type(excfiles),intent(in) :: BS_files
 type(excparam),intent(in) :: BSp
!arrays
 real(dp),intent(out) :: dos_exc(nomega)
 complex(dpc),intent(in) :: opt_cvk(lomo_min:max_band,lomo_min:max_band,nkbz,nsppol,BSp%nq),omega(nomega)
 complex(dpc),intent(out) :: eps_exc(nomega,BSp%nq)

!Local variables ------------------------------
!scalars
 integer :: mi,it,ii,ib_v,ib_c,ik_bz,exc_size_read,nstates_read,eig_unt !,fform
 integer :: exc_size,iq,spin,tr_idx,tar_idx,nstates,iw,ll,ierr
 real(dp) :: fact,arg
 complex(dpc) :: eps,fam,famp
 character(len=500) :: msg,errmsg
 character(len=fnlen) :: filbseig
 logical :: do_lifetime
 !type(Hdr_type) :: tmp_Hdr
!arrays
 complex(dpc),allocatable :: Ami(:),exc_ene(:),Sm1mi(:)
 complex(dpc),allocatable :: msfap(:,:),fa(:,:),fap(:,:)

!************************************************************************

 call wrtout(std_out," Calculating absorption strength with full coupling","COLL")

 if (nsppol==2) then
   MSG_WARNING("nsppol==2 is still under development")
 end if

 ! Rank of the entire excitonic Hamiltonian including the coupling block.
 exc_size = 2*SUM(BSp%nreh); if (nsppol==2) exc_size = 2*(SUM(BSp%nreh) + BSp%nreh(2))
 nstates  = BSp%nstates

 ! TODO: four_pi comes from the bare Coulomb term hence the
 ! present implementation is not compatible with the cutoff technique.
 ! factor two is due to the occupation factors.
 fact=four_pi/(ucvol*nkbz); if (nsppol==1) fact=two*fact

 if (BS_files%in_eig /= BSE_NOFILE) then
   filbseig = BS_files%in_eig
 else
   filbseig = BS_files%out_eig
 end if

 call wrtout(std_out," Reading excitonic eigenstates from file: "//trim(filbseig),"COLL")
 if (open_file(filbseig,msg,newunit=eig_unt,form="unformatted", status="old", action="read") /= 0) then
   MSG_ERROR(msg)
 end if

 read(eig_unt, err=10, iomsg=errmsg) do_lifetime

 if (do_lifetime) then
   ABI_CHECK(.not. do_lifetime, "Finite lifetime with coupling is not supported yet !")
 end if

 read(eig_unt, err=10, iomsg=errmsg) exc_size_read, nstates_read
 ABI_CHECK(exc_size_read==exc_size,"wrong file")
 ABI_CHECK(nstates_read==nstates,"Partial diago not supported yet")
 !
 ! Read eigenvalues
 ABI_MALLOC(exc_ene,(nstates))
 read(eig_unt, err=10, iomsg=errmsg) exc_ene(:)

 ABI_MALLOC(fa,(nstates,BSp%nq))
 ABI_MALLOC(fap,(nstates,BSp%nq))
 ABI_MALLOC_OR_DIE(Ami,(exc_size), ierr)

 do mi=1,nstates ! Loop on excitonic eigenvalues mi
   read(eig_unt, err=10, iomsg=errmsg) Ami(:)

   do iq=1,BSp%nq
     fam  = czero
     famp = czero
     do spin=1,nsppol
       do it=1,BSp%nreh(spin) ! Loop over transition t = (k,v,c)
         ik_bz = Bsp%Trans(it,spin)%k
         ib_v  = Bsp%Trans(it,spin)%v
         ib_c  = Bsp%Trans(it,spin)%c
         tr_idx  = it + (spin-1)*Bsp%nreh(1)
         if (nsppol==1) then
           tar_idx = it + Bsp%nreh(1)
         else
           if (spin==1) tar_idx = it + SUM(Bsp%nreh)
           if (spin==2) tar_idx = it + 2*Bsp%nreh(1)+Bsp%nreh(2)
         end if

         fam = fam + CONJG(opt_cvk(ib_c,ib_v,ik_bz,spin,iq)) * Ami(tr_idx) &
&                  + CONJG(opt_cvk(ib_v,ib_c,ik_bz,spin,iq)) * Ami(tar_idx)

         famp = famp - opt_cvk(ib_c,ib_v,ik_bz,spin,iq) * CONJG(Ami(tr_idx)) &
&                    + opt_cvk(ib_v,ib_c,ik_bz,spin,iq) * CONJG(Ami(tar_idx))
       end do
     end do
     ! Save results.
     fa (mi,iq) = fam
     fap(mi,iq) = famp
   end do
 end do ! mi

 ABI_FREE(Ami)

 ! Read O{-1} and sum over the eigenstates.
 ABI_MALLOC(msfap,(nstates,BSp%nq))
 ABI_MALLOC(Sm1mi,(nstates))

 do mi=1,nstates
   read(eig_unt, err=10, iomsg=errmsg) Sm1mi
   Sm1mi = DCONJG(Sm1mi) ! This gives the row since O^{-1} is Hermitian.
   do iq=1,BSp%nq
     msfap(mi,iq) = xdotu(exc_size,Sm1mi,1,fap(:,iq),1)
   end do
 end do

 ABI_FREE(Sm1mi)

 close(eig_unt, err=10, iomsg=errmsg)
 !
 ! === Calculate excitonic epsilon with coupling ===
 do iq=1,BSp%nq
   !
   do ii=1,nomega
     eps = czero
     do mi=1,nstates ! sum over all exciton eigenstates
       eps = eps - fa(mi,iq) * msfap(mi,iq) / (exc_ene(mi) - omega(ii))
     end do
     eps_exc(ii,iq) = one + fact * eps
   end do
   !
 end do

 ABI_FREE(fa)
 ABI_FREE(msfap)
 ABI_FREE(fap)
 !
 ! The excitonic DOS.
 dos_exc=zero
 do ll=1,nstates ! Sum over the calculate excitonic eigenstates.
   do iw=1,nomega
     arg = DBLE(omega(iw) - exc_ene(ll))
     dos_exc(iw) = dos_exc(iw) + gaussian(arg, Bsp%broad)
   end do
 end do

 ABI_FREE(exc_ene)

 return

10 continue
 MSG_ERROR(errmsg)

end subroutine exc_eps_coupling
!!***

!----------------------------------------------------------------------

!!****f* m_exc_spectra/exc_write_tensor
!! NAME
!!  exc_write_tensor
!!
!! FUNCTION
!!  This routine drives the writing of complex dielectric tensor
!!
!! INPUTS
!! BSp<excparam>=Bethe-Salpeter Parameters.
!! what= "EXC_TSR_CART" or "EXC_TSR_RED"
!!       "RPA_NLF_TSR_CART" or "RPA_NLF_TSR_RED"
!!       "GW_NLF_TSR_CART" or "GW_NLF_TSR_RED"
!!
!! OUTPUT
!!  Only writing.
!!
!! SIDE EFFECTS
!!  tensor(BSp%nomega,6) = Complex dielectric tensor to be written
!!
!! PARENTS
!!      m_haydock
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine exc_write_tensor(BSp,BS_files,what,tensor)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: what
 type(excparam),intent(in) :: BSp
 type(excfiles),intent(in) :: BS_files
!arrays
 complex(dpc),intent(in) :: tensor(BSp%nomega,6)

!Local variables ------------------------------
!scalars
 integer :: io,iq,funt
 real(dp) :: omegaev
!arrays
 character(len=500) :: lf_type,block_type,wgg_type,frm,str_type, msg
 character(len=fnlen) :: fname

!************************************************************************

 fname = strcat(BS_files%out_basename,'_',toupper(what))
 if (open_file(fname,msg,newunit=funt,form="formatted", action="write") /= 0) then
   MSG_ERROR(msg)
 end if

 select case (toupper(what))

 case ("EXC_TSR_CART")
   write(funt,'("# Complex dielectric tensor (cart. coord.) obtained with the BS equation.")')

   lf_type = 'WITHOUT LOCAL FIELD EFFECTS'
   if (BSp%exchange_term>0) lf_type='LOCAL FIELD EFFECTS INCLUDED'
   call bsp_calctype2str(Bsp,str_type)
   write(funt,'("# ",a,"     " ,a)') TRIM(str_type), TRIM(lf_type)

   block_type = 'RESONANT-ONLY calculation'
   if (BSp%use_coupling>0) block_type = 'RESONANT+COUPLING calculation'
   write(funt,'("# ",a)') TRIM(block_type)

   if (BSp%use_coulomb_term) then
     wgg_type = "Coulomb term constructed with full W(G1,G2)"
     if ( BSp%use_diagonal_Wgg ) wgg_type = "Coulomb term constructed with diagonal approximation W(G1,G1)"
     write(funt,'("# ",a)') TRIM(wgg_type)
   end if

   write(funt,'(a,f7.4,a)')'# Scissor operator energy = ',BSp%mbpt_sciss*Ha_eV,' [eV]'

 case ("EXC_TSR_RED")
   write(funt,'("# Complex dielectric tensor (red. coord.) obtained with the BS equation.")')

   lf_type = 'WITHOUT LOCAL FIELD EFFECTS'
   if (BSp%exchange_term>0) lf_type='LOCAL FIELD EFFECTS INCLUDED'
   call bsp_calctype2str(Bsp,str_type)
   write(funt,'("# ",a,"     " ,a)') TRIM(str_type), TRIM(lf_type)

   block_type = 'RESONANT-ONLY calculation'
   if (BSp%use_coupling>0) block_type = 'RESONANT+COUPLING calculation'
   write(funt,'("# ",a)') TRIM(block_type)

   if (BSp%use_coulomb_term) then
     wgg_type = "Coulomb term constructed with full W(G1,G2)"
     if ( BSp%use_diagonal_Wgg ) wgg_type = "Coulomb term constructed with diagonal approximation W(G1,G1)"
     write(funt,'("# ",a)') TRIM(wgg_type)
   end if

   write(funt,'(a,f7.4,a)')'# Scissor operator energy = ',BSp%mbpt_sciss*Ha_eV,' [eV]'

 case ("RPA_NLF_TSR_CART")
   write(funt,'("# RPA complex dielectric tensor (cart. coord.) without local fields")')

 case ("RPA_NLF_TSR_RED")
   write(funt,'("# RPA complex dielectric tensor (red. coord.) without local fields")')

 case ("GW_NLF_TSR_CART")
   write(funt,'("# GW complex dielectric tensor (cart. coord.) without local field effects ")')
   write(funt,'(a,f7.4,a)')'# Scissor operator energy = ',BSp%mbpt_sciss*Ha_eV,' [eV]'

 case ("GW_NLF_TSR_RED")
   write(funt,'("# GW complex dielectric tensor (red. coord.) without local field effects ")')
   write(funt,'(a,f7.4,a)')'# Scissor operator energy = ',BSp%mbpt_sciss*Ha_eV,' [eV]'

 case default
   MSG_ERROR("Unknown value for what: "//TRIM(what))
 end select
 !
 ! Paramaters common to the different calculations.
 if (BSp%algorithm /= BSE_ALGO_HAYDOCK) then
   write(funt,'(a,i0)')"# nstates included in the diagonalization = ",BSp%nstates
 end if

 if (BSp%algorithm == BSE_ALGO_HAYDOCK) then
   write(funt,'(a,2f7.4)')'# Tolerance = ',BSp%haydock_tol
 end if

 write(funt,'(a,i0)')"# npweps  = ",BSp%npweps
 write(funt,'(a,i0)')"# npwwfn  = ",BSp%npwwfn
 write(funt,'(a,i0)')"# nbands  = ",BSp%nbnds
 write(funt,'(a,i0)')"# loband  = ",BSp%lomo_spin(1)
 if (Bsp%nsppol==2) write(funt,'(a,i0)')"# loband(spin=2) = ",BSp%lomo_spin(2)
 write(funt,'(a,i0)')"# nkibz   = ",BSp%nkibz
 write(funt,'(a,i0)')"# nkbz    = ",BSp%nkbz
 write(funt,'(a,f7.4,a)')'# Lorentzian broadening = ',BSp%broad*Ha_eV,' [eV]'

 !
 ! Write tensor.
 write(funt,'(3a)') "# omega [eV] RE(eps_11) IM(eps_11) RE(eps_22)", &
&  "IM(eps_22) RE(eps_33) IM(eps_33) RE(eps_12) IM(eps_12)", &
&  "RE(eps_13) IM(eps_13) RE(eps_23) IM(eps_23))"
 write(frm,*) '(f7.3,12es14.6)'
 do io=1,BSp%nomega
   omegaev = DBLE(BSp%omega(io))*Ha_eV
   write(funt,frm) omegaev,(tensor(io,iq), iq=1,6)
 end do

 close(funt)

end subroutine exc_write_tensor
!!***

!----------------------------------------------------------------------

!!****f* m_exc_spectra/mdfs_ncwrite
!! NAME
!! mdfs_ncwrite
!!
!! FUNCTION
!!  Writes the MDF.nc file with the final results.
!!
!! INPUTS
!!  ncid =NC file handle
!!  Bsp<excparam>=Data type gathering the paramenters used for the Bethe-Salpeter calculation.
!!  eps_exc = Excitonic MDF
!!  eps_rpanlf = KS-RPA MDF without local-field effects.
!!  eps_gwnlf = GW-RPA MDF without local-field effects.
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!      m_exc_spectra,m_haydock
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine mdfs_ncwrite(ncid,Bsp,eps_exc,eps_rpanlf,eps_gwnlf)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 type(excparam),intent(in) :: BSp
!arrays
 complex(dpc),target,intent(in) :: eps_exc(BSp%nomega,BSp%nq)
 complex(dpc),target,intent(in) :: eps_rpanlf(BSp%nomega,BSp%nq)
 complex(dpc),target,intent(in) :: eps_gwnlf(BSp%nomega,BSp%nq)

!Local variables-------------------------------
!scalars
#ifdef HAVE_NETCDF
 integer :: ncerr
 real(dp), ABI_CONTIGUOUS pointer :: rvals(:,:,:)

! *************************************************************************
 ! =========================
 ! === Write the dimensions
 ! =========================

 ncerr = nctk_defnwrite_ivars(ncid, [character(len=nctk_slen) :: &
&  "mdf_version", "nsppol", "npwwfn", "npweps", "nkibz", "nkbz",&
&  "nkibz_iterp", "nkbz_interp", "wtype", "interp_mode"],&
&  [1, Bsp%nsppol, Bsp%npwwfn, Bsp%npweps, Bsp%nkibz, Bsp%nkbz, &
&   Bsp%nkibz_interp, Bsp%nkbz_interp,Bsp%wtype, Bsp%interp_mode])
 NCF_CHECK(ncerr)

 ncerr = nctk_defnwrite_dpvars(ncid, [character(len=nctk_slen) :: &
&  "ecutwfn", "ecuteps", "mbpt_sciss", "broad", "eps_inf"],&
&  [Bsp%ecutwfn, Bsp%ecuteps, Bsp%mbpt_sciss, Bsp%broad, Bsp%eps_inf])
 NCF_CHECK(ncerr)

 ncerr = nctk_def_dims(ncid, [nctkdim_t("two", 2), nctkdim_t("three", 3), nctkdim_t("number_of_qpoints", Bsp%nq),&
   nctkdim_t("number_of_frequencies", Bsp%nomega), nctkdim_t("number_of_spins", bsp%nsppol)], defmode=.True.)
 NCF_CHECK(ncerr)

! Define variables.

!arrays
 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t('qpoints', "dp", 'three, number_of_qpoints'),&
   nctkarr_t('wmesh', "dp", 'number_of_frequencies'),&
   nctkarr_t('nreh', "i", "number_of_spins"),&
   nctkarr_t('lomo_spin', "i", "number_of_spins"),&
   nctkarr_t('humo_spin', "i", "number_of_spins"),&
   nctkarr_t('exc_mdf', "dp", 'two, number_of_frequencies, number_of_qpoints'),&
   nctkarr_t('rpanlf_mdf', "dp", 'two, number_of_frequencies, number_of_qpoints'),&
   nctkarr_t('gwnlf_mdf', "dp", 'two, number_of_frequencies, number_of_qpoints')])
 NCF_CHECK(ncerr)

! Write data.
 NCF_CHECK(nctk_set_datamode(ncid))
 NCF_CHECK(nf90_put_var(ncid, vid('qpoints'), Bsp%q))
 NCF_CHECK(nf90_put_var(ncid, vid('nreh'), bsp%nreh))
 NCF_CHECK(nf90_put_var(ncid, vid('lomo_spin'), bsp%lomo_spin))
 NCF_CHECK(nf90_put_var(ncid, vid('humo_spin'), bsp%humo_spin))

 ! Write frequency in mesh in eV.
 NCF_CHECK(nf90_put_var(ncid, vid('wmesh'), REAL(Bsp%omega)*Ha_eV))

 call c_f_pointer(c_loc(eps_exc(1,1)), rvals, shape=[2, bsp%nomega, bsp%nq])
 NCF_CHECK(nf90_put_var(ncid, vid('exc_mdf'), rvals))

 call c_f_pointer(c_loc(eps_rpanlf(1,1)), rvals, shape=[2, bsp%nomega, bsp%nq])
 NCF_CHECK(nf90_put_var(ncid, vid('rpanlf_mdf'), rvals))

 call c_f_pointer(c_loc(eps_gwnlf(1,1)), rvals, shape=[2, bsp%nomega, bsp%nq])
 NCF_CHECK(nf90_put_var(ncid, vid("gwnlf_mdf"), rvals))

#else
 MSG_ERROR("ETSF-IO support is not activated.")
#endif

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end subroutine mdfs_ncwrite
!!***

!!****f* m_exc_spectra/check_kramerskronig
!! NAME
!!  check_kramerskronig
!!
!! FUNCTION
!!   check Kramers Kronig
!!   \int_0^\infty d\omega' frac{\omega'}{\omega'^2 - \omega^2}
!!   Im \epsilon(\omega') = Re \epsilon(\omega)
!!
!! INPUTS
!!  n=Number of frequency points.
!!  eps(n)=Dielectric function.
!!  o(n)=Frequency mesh.
!!
!! OUTPUT
!!  Only checking.
!!
!! PARENTS
!!      m_exc_spectra
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine check_kramerskronig(n,o,eps)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
!arrays
 complex(dpc),intent(in) :: eps(n)
 real(dp),intent(in) :: o(n)

!Local variables ------------------------------
!scalars
 integer :: ii,ip
 real(dp) :: omega,omegap,domega,kk,kkrms,eav
 complex(dpc) c
 character(len=500) :: msg
!arrays
 real(dp) :: e1kk(n)
 complex(dpc) :: intg(n)

!************************************************************************
! init jmb
 e1kk=zero
 intg=(zero,zero)

! calculate domega step and verify all
 domega = (o(n) - o(1)) / (n-1)

 do ii=2,n
  if (domega-(o(ii)-o(ii-1)) > tol3) then
    MSG_WARNING("Frequency mesh not linear. Returning")
    return
  end if
 end do

 if(o(1) > 0.1/Ha_eV) then
   MSG_WARNING("First frequency is not zero. Returning")
   return
 end if

 if (aimag(eps(n)) > 0.1) then
   write(msg,'(a,f12.6,3a,f12.6,2a)')&
&   ' Im epsilon for omega= ',o(n)*Ha_eV,'eV',ch10,&
&   ' is not yet zero, epsilon_2= ',aimag(eps(n)),ch10,&
&   ' Kramers Kronig test could give wrong results. '
   MSG_WARNING(msg)
 end if

! Fill array for kramers kronig.
 do ii=1,n
   omega=o(ii)
   c = (0.0,0.0)
   do ip=1,n
     if(ip == ii) cycle
     omegap = o(ip)
     c = c + omegap / (omegap**2-omega**2) * aimag(eps(ip))
   end do
   e1kk(ii) = one + two/pi * domega*real(c)
 end do

!perform kramers kronig with simpson integration
 do ii=1,n
   omega=o(ii)
   do ip=1,n
     if (ip==ii) cycle
     omegap = o(ip)
     intg(ip) = omegap / (omegap**2 - omega**2) * aimag(eps(ip))
   end do
   c = simpson_cplx(n,domega,intg)
   e1kk(ii) = one + two/pi * real(c)
 end do

!verify kramers kronig
 eav=zero; kk=zero; kkrms=zero
 do ii=1,n
   kk = kk + abs(real(eps(ii)) - e1kk(ii))
   kkrms = kkrms +(real(eps(ii)) - e1kk(ii))*(real(eps(ii)) - e1kk(ii))
   eav = eav + abs(real(eps(ii)))
 end do

 eav = eav/n
 kk = (kk/n)/eav
 kkrms = (kkrms/n) / (eav*eav)

 kk = abs(real(eps(1)) - e1kk(1)) / real(eps(1))

! write data
 write(msg,'(a,f7.2,a)')" The Kramers-Kronig is verified within ",100*kk,"%"
 call wrtout(std_out,msg,"COLL")

! write(std_out,'("# Kramers Kronig calculation of epsilon1")')
! write(std_out,'("# omega   epsilon1  epsilon1kk")')
! do ii=1,n
!   write(std_out,'(f7.3,2e15.7)') o(ii)*Ha_eV, real(eps(ii)), e1kk(ii)
! end do

end subroutine check_kramerskronig
!!***

!----------------------------------------------------------------------

!!****f* m_exc_spectra/check_fsumrule
!! NAME
!!  check_fsumrule
!!
!! FUNCTION
!!   check f-sum rule
!!   \int_0^\infty d\omega \omega Im \epsilon_GG'(q,\omega) =
!!   = \frac{1}{2} \pi \omega_p^2  \frac{\rho(G-G')}{\rho(0)}
!!   versor(q+G) \dot versor(q+G')
!!   for q = G = G' = 0, it reads:
!!   \int_0^\infty d\omega \omega Im \epsilon_00(q=0,\omega) =
!!   = \pi \omega_p^2 / 2
!!   calculate only the second one
!!   calculate the integral to evaluate an omega_plasma^eff to compare with omega_plasma
!!
!! INPUTS
!!  n=Number of frequencies.
!!  o(n)=Frequency mesh.
!!  e2(n)=imaginary part of epsilon_00
!!  omegaplasma=Drude plasma frequency.
!!
!! OUTPUT
!!  Only checking.
!!
!! PARENTS
!!      m_exc_spectra
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine check_fsumrule(n,o,e2,omegaplasma)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n
 real(dp),intent(in) :: omegaplasma
!arrays
 real(dp),intent(in) :: o(n),e2(n)

!Local variables ------------------------------
!scalars
 integer :: ii,ip
 real(dp) :: omegap,domega,integral,omegaplasmaeff,fsumrule
 character(len=500) :: msg
!arrays
 complex(dpc) :: intg(n)

!************************************************************************

! calculate domega step and verify
 domega = (o(n) - o(1)) / (n-1)

 do ii=2,n
   if (domega-(o(ii)-o(ii-1)) > tol3) then
     MSG_WARNING("Frequency mesh not linear. Returning")
     return
   end if
 end do

 if (o(1) > 0.1/Ha_eV) then
   MSG_WARNING("First frequency is not zero. Returning")
   return
 end if

 if (e2(n) > 0.1) then
   write(msg,'(a,f12.6,3a,f12.6,2a)')&
&   ' Im epsilon for omega= ',o(n)*Ha_eV,' eV ',ch10,&
&   ' is not yet zero, epsilon_2= ',e2(n),ch10,&
&   ' f-sum rule test could give wrong results.'
   MSG_WARNING(msg)
 end if

! integrate to obtain f-sum rule
 integral=zero
 do ip=1,n
   omegap=o(ip)
   integral = integral + omegap * e2(ip)
 end do
 integral = domega * integral

!integrate with simpson to obtain f-sum rule
 do ip = 1, n
   omegap = o(ip)
   intg(ip) = omegap * e2(ip)
 end do

 integral = real(simpson_cplx(n,domega,intg))
 if(integral < 0) then
   MSG_ERROR("The integral of the imaginary of dielectric function is negative !!!")
 else
   omegaplasmaeff = sqrt(integral*two/pi)
 end if

 fsumrule = abs((omegaplasmaeff - omegaplasma)) / omegaplasma

! write data
 write(msg,'(3(a,f6.2,2a))')&
&  " omega_plasma     = ",omegaplasma*Ha_eV,   " [eV]",ch10,&
&  " omega_plasma^eff = ",omegaplasmaeff*Ha_eV," [eV]",ch10,&
&  " the f-sum rule is verified within ",fsumrule*100,"%",ch10
 call wrtout(std_out,msg,"COLL")

end subroutine check_fsumrule
!!***

END MODULE m_exc_spectra
