!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_elphon
!! NAME
!! m_elphon
!!
!! FUNCTION
!! This routine extracts the electron phonon coupling matrix
!! elements and calculates related properties - Tc, phonon linewidths...
!!
!! COPYRIGHT
!! Copyright (C) 2004-2018 ABINIT group (MVer,BXu,MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

module m_elphon

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_elphon
 use m_profiling_abi
 use m_kptrank
 use m_errors
 use m_xmpi
 use m_hdr
 use m_ebands

 use m_io_tools,        only : open_file, is_open, get_unit
 use m_time,            only : timein
 use m_numeric_tools,   only : wrap2_pmhalf, simpson, simpson_int
 use m_pptools,         only : printvtk
 use m_dynmat,          only : ftgam_init, ftgam
 use m_geometry,        only : phdispl_cart2red, ifc_fourq
 use m_crystal,         only : crystal_t
 use m_ifc,             only : ifc_type
 use m_nesting,         only : mknesting, bfactor
 use m_anaddb_dataset,  only : anaddb_dataset_type
 use m_eliashberg_1d,   only : eliashberg_1d

 implicit none

 private
!!***

 public :: elphon

contains

!!****f* m_elphon/elphon
!!
!! NAME
!! elphon
!!
!! FUNCTION
!! This routine extracts the electron phonon coupling matrix
!! elements and calculates related properties - Tc, phonon linewidths...
!!
!! INPUTS
!!   anaddb_dtset=dataset with input variables
!!     anaddb_dtset%a2fsmear = smearing for alpha2F function
!!     anaddb_dtset%brav = type of Bravais lattice
!!     anaddb_dtset%elphsmear = smearing width for gaussian integration
!!           or buffer in energy for calculations with tetrahedra (telphint=0)
!!     anaddb_dtset%elph_fermie = input value of Fermi energy
!!           0 means use value from wfk file
!!     anaddb_dtset%enunit = governs the units to be used for the output of
!!           the phonon frequencies and e-ph quantities
!!     anaddb_dtset%gkk2write= flag to write out gkk2 matrix elements to disk
!!     anaddb_dtset%gkk_rptwrite= flag to write out real space gkk_rpt matrix elements to disk
!!     anaddb_dtset%gkqwrite= flag to write out gkq matrix elements to disk
!!     anaddb_dtset%ep_b_min= first band taken into account in FS integration (if telphint==2)
!!     anaddb_dtset%ep_b_max= last band taken into account in FS integration (if telphint==2)
!!     anaddb_dtset%prtfsurf = integer flag for the output of the Fermi surface (XCrysden file format)
!!     anaddb_dtset%prtnest = integer flag for the calculation of the nesting function
!!     anaddb_dtset%ifcflag = flag for IFC matrices in anaddb calling routine
!!           the IFCs are presumed to be known!
!!     anaddb_dtset%ifltransport= flag for transport properties (no=0: yes_LOVA=1; yes_nonLOVA=2 )
!!     anaddb_dtset%kptrlatt=kpoint grid generating vectors, as in abinit
!!     anaddb_dtset%kptrlatt_fine=kpoint grid generating vectors, for fine grid used in FS integration
!!     anaddb_dtset%mustar = parameter for Coulombic pseudo-potential in McMillan T_c calculation
!!     anaddb_dtset%ngqpt(3)=integers defining the number of points in the qpt sampling
!!     anaddb_dtset%nqpath=number of vertices in the path in reciprocal space, for band structure
!!           and phonon linewidth output
!!     anaddb_dtset%nqshft= number of shift vectors for defining the sampling of q points
!!     anaddb_dtset%ntemper = number of temperature points to calculate, from tempermin to
!!           tempermin+ntemper*temperinc
!!     anaddb_dtset%qpath=vertices in the path in reciprocal space, for band structure
!!           and phonon linewidth output
!!     anaddb_dtset%q1shft(3,4) =qpoint shifts considered
!!     anaddb_dtset%telphint = flag for integration over the FS with 0=tetrahedra 1=gaussians
!!     anaddb_dtset%tempermin = minimum temperature at which resistivity etc are calculated (in K)
!!     anaddb_dtset%temperinc = interval temperature grid on which resistivity etc are calculated (in K)
!!     anaddb_dtset%ep_keepbands = flag to keep gamma matrix dependence on electronic bands
!! Cryst<crystal_t>=data type gathering info on the crystalline structure.
!! Ifc<ifc_type>=Object containing the interatomic force constants.
!!     atmfrc  = inter-atomic force constants from anaddb
!!     rpt(3,nprt) =canonical positions of R points in the unit cell
!!     nrpt =number of real space points used to integrate IFC (for interpolation of dynamical matrices)
!!     wghatm(natom,natom,nrpt) =Weight for the pair of atoms and the R vector
!! filnam(7)=character strings giving file names
!! comm=MPI communicator.
!!
!! OUTPUT
!!
!! NOTES
!!  inspired to a large extent by epcouple.f from the DecAFT package by J. Kay Dewhurst
!!  most inputs taken from mkifc.f
!!  in anaddb anaddb_dtset%ifcflag must be 1 such that the IFC are calculated in atmfrc prior to calling elphon
!!
!!  brav not taken into account propely in all of the code. (MG?)
!!
!!  could choose to make a full 3 dimensional kpt array (:,:,:). Easier for many operations
!!
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!      complete_gamma,complete_gamma_tr,copy_kptrank,d2c_weights,ebands_free
!!      ebands_update_occ,eliashberg_1d,elph_ds_clean,elph_k_procs
!!      elph_tr_ds_clean,ep_fs_weights,ep_setupqpt,ftgam,ftgam_init
!!      get_all_gkk2,get_all_gkq,get_all_gkr,get_fs_bands,get_nv_fs_en
!!      get_nv_fs_temp,get_rank_1kpt,get_tau_k,get_veloc_tr,hdr_bcast
!!      hdr_fort_read,hdr_free,integrate_gamma,integrate_gamma_tr
!!      integrate_gamma_tr_lova,mka2f,mka2f_tr,mka2f_tr_lova,mka2fqgrid
!!      mkfskgrid,mknesting,mkph_linwid,mkqptequiv,order_fs_kpts,outelph
!!      printvtk,rchkgsheader,read_el_veloc,symkpt,timein,wrap2_pmhalf,wrtout
!!      xmpi_bcast
!!
!! NOTES
!!
!! SOURCE

subroutine elphon(anaddb_dtset,Cryst,Ifc,filnam,comm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'elphon'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_77_ddb
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(anaddb_dataset_type),intent(inout) :: anaddb_dtset
 type(crystal_t),intent(in) :: Cryst
 type(ifc_type),intent(inout) :: Ifc
 integer,intent(in) :: comm
!arrays
 character(len=fnlen),intent(in) :: filnam(7)

!Local variables-------------------------------
!scalars
 integer,parameter :: timrev2=2,space_group0=0,master=0
 integer :: ikpt_fine,ierr,unitgkk, unit_epts,iband,ibandp,ii
 integer :: ikpt,jkpt,kkpt, ik1,ik2,ik3,nk1, nk2, nk3
 integer :: iqpt,isppol,n1wf,nband,natom,onegkksize
 integer :: timrev,unitfskgrid,qtor,idir,iFSkpq,symrankkpt,ikpt_irr
 integer :: ep_prt_wtk ! eventually to be made into an input variable
 integer :: fform,ie,ie1,ie2,i_start,i_end
 integer :: ssp,s1,s2,tmp_nenergy, top_vb,nproc,me
 integer :: nkpt_tmp
 real(dp) :: max_occ,realdp_ex,res !,ss
 real(dp) :: tcpu, twall, tcpui, twalli
 real(dp) :: e1, e2, btocm3,diff, omega_max
 real(dp) :: e_vb_max, e_cb_min, etemp_vb
 logical :: make_gkk2,use_afm,use_tr
 character(len=500) :: message
 character(len=fnlen) :: fname,elph_base_name,ddkfilename,gkk_fname
 character(len=fnlen) :: nestname
 type(elph_tr_type) :: elph_tr_ds
 type(elph_type) :: elph_ds
 type(hdr_type) :: hdr,hdr1
 type(ebands_t) :: Bst
!arrays
 integer :: s1ofssp(4), s2ofssp(4)
 integer :: qptrlatt(3,3),kptrlatt_fine(3,3)
 integer,allocatable :: indkpt1(:)
 integer,allocatable :: FSfullpqtofull(:,:)
 integer,allocatable :: qpttoqpt(:,:,:)
 integer,allocatable :: pair2red(:,:), red2pair(:,:)
 !real(dp) :: acell_in(3),rprim_in(3,3),rprim(3,3),acell(3),
 real(dp) :: kpt(3),shiftk(3)
 real(dp),allocatable :: wtk_fullbz(:),wtk_folded(:)
 real(dp),allocatable :: a2f_1d(:),dos_phon(:)
 real(dp),allocatable :: eigenGS(:,:,:),eigenGS_fine(:,:,:)
 real(dp),allocatable :: v_surf(:,:,:,:,:,:)
 real(dp),allocatable :: tmp_veloc_sq1(:,:), tmp_veloc_sq2(:,:)
 real(dp),allocatable :: coskr(:,:), sinkr(:,:)

! *************************************************************************

 write(message, '(a,a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
& ' Properties based on electron-phonon coupling ',ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 call timein(tcpui,twalli)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-begin elphon at tcpu',tcpui,'  and twall',twalli,' sec'
 call wrtout(std_out,message,'COLL')

 nproc = xmpi_comm_size(comm); me = xmpi_comm_rank(comm)

 write(message, '(a,i0,a,i0)' )'- running on ', nproc,'  cpus me = ', me
 call wrtout(std_out,message,'PERS')
 write(std_out,*) message

!==================================
!Initialization of some variables
!==================================

 if (master == me) then
   gkk_fname = filnam(5)
   if (open_file(gkk_fname,message,newunit=unitgkk,form="unformatted",status="old",action="read") /=0) then
     MSG_ERROR(message)
   end if
 end if

 elph_base_name=trim(filnam(2))//"_ep"
 ddkfilename=trim(filnam(7))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 natom = Cryst%natom
 elph_ds%mustar       = anaddb_dtset%mustar        ! input mustar
 elph_ds%nbranch      = 3*natom                    ! number of phonon modes = 3 * natom
 elph_ds%natom        = natom                      !
 elph_ds%ep_keepbands = anaddb_dtset%ep_keepbands  ! flag to sum over bands
 elph_ds%a2fsmear     = anaddb_dtset%a2fsmear      ! smearing for Eliashberg functions
 elph_ds%elphsmear    = anaddb_dtset%elphsmear     ! smearing for Eliashberg functions
 elph_ds%ep_b_min     = anaddb_dtset%ep_b_min
 elph_ds%ep_b_max     = anaddb_dtset%ep_b_max
 elph_ds%telphint     = anaddb_dtset%telphint
 elph_ds%kptrlatt     = anaddb_dtset%kptrlatt
 elph_ds%kptrlatt_fine= anaddb_dtset%kptrlatt_fine
 elph_ds%tempermin    = anaddb_dtset%tempermin
 elph_ds%temperinc    = anaddb_dtset%temperinc
 elph_ds%ntemper      = anaddb_dtset%ntemper
 elph_ds%use_k_fine   = anaddb_dtset%use_k_fine
 elph_ds%ep_int_gkk   = anaddb_dtset%ep_int_gkk
 elph_ds%ep_nspline   = anaddb_dtset%ep_nspline
 elph_ds%ep_scalprod  = anaddb_dtset%ep_scalprod
 elph_ds%prtbltztrp   = anaddb_dtset%prtbltztrp

 elph_ds%tuniformgrid = 1
 elph_ds%na2f         = 400                        ! maximum number of Matsubara frequencies.
 elph_ds%ep_lova      = 0                          ! 1 for lova and 0 for general
 elph_ds%nenergy      = 8
 btocm3 = 1.4818474347690475d-25

!The nenergy needs to be 1) large enough to converge the integral, 2) greater
!than the max phonon energy.
!elph_ds%nenergy      = INT(8*(anaddb_dtset%tempermin+anaddb_dtset%ntemper*anaddb_dtset%temperinc)/ &
!&                              (anaddb_dtset%tempermin+anaddb_dtset%temperinc))  ! number of energy levels

 write(message,'(a,i6)')' The initial number of energy levels above/below Ef is set to be :',elph_ds%nenergy
 call wrtout(std_out,message,'COLL')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!The precise number used depends on the value of Tc:
!they span $w_n = (2n+1) \pi T_c$  where $abs(w_n) < w_{cutoff}$
!ie $|n| < n_{cutoff} = ( \frac{w_{cutoff}}{\pi T_c} ) / 2$

!save gkk data for full kpoints to file on disk

 elph_ds%gkqwrite     = anaddb_dtset%gkqwrite
 elph_ds%gkk_rptwrite = anaddb_dtset%gkk_rptwrite
 elph_ds%gkk2write    = anaddb_dtset%gkk2write

!This should never be turned off: symmetrization of elphon matrix elements in complete_gkk. See get_all_gkq
 elph_ds%symgkq=anaddb_dtset%symgkq

 elph_ds%elph_base_name = trim(elph_base_name)

 !MG: @Matthieu: Why this? Now we should always use the value of rprim and acell reported in IFC
 !rprim_in  = Ifc%rprim
 !acell_in = Ifc%acell

!normalize input rprim and acell.
 !do ii=1,3
 !  ss = sqrt(rprim_in(1,ii)**2+rprim_in(2,ii)**2+rprim_in(3,ii)**2)
 !  rprim(:,ii) = rprim_in(:,ii)/ss
 !  acell(ii) = acell_in(ii) * ss
 !end do

!make dimension-ful rprimd and gprimd for transformation of derivatives to cartesian coordinates.
 !call mkrdim(acell,rprim,rprimd)
 !call matr3inv(rprimd,gprimd)

 !rprimd = cryst%rprimd
 !gprimd = cryst%gprimd

!===================
!Check some inputs
!===================
 if (Cryst%nsym==1) then
   write (message,'(7a)')ch10,&
&   ' elphon: COMMENT- ',ch10,&
&   ' Symmetries are not used! ',ch10,&
&   ' Full matrix elements must be supplied for all perturbations and qpoints!',ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   if ( ANY( ABS(Cryst%tnons(:,1)) > tol10) ) then
     MSG_ERROR('nsym==1 but the symmetry is not the identity')
   end if
 end if

 if (anaddb_dtset%ifcflag/=1) then
   write(message,'(a,i0)')&
&   ' ifcflag should be set to 1 since the IFC matrices are supposed to exist but ifcflag= ',anaddb_dtset%ifcflag
   MSG_ERROR(message)
 end if

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon begin setup after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

!=================================
!Set up the full grid of qpoints
!=================================
!use time reversal symmetry always when possible for kpoint reduction,
!and suppose it has been used in WF generation
!not used for the moment: values are always taken from input files.
 timrev = 1
 call ep_setupqpt(elph_ds,cryst,anaddb_dtset,qptrlatt,timrev)

!====================================
!Read the GS header of the GKK file
!this will give the phon grid of k
!and the Fermi surface integration weights
!====================================
 call wrtout (std_out,' elphon: reading and checking the GS header of the GKK file','COLL')

 if (master == me) then
   call rchkGSheader(hdr,natom,nband,unitgkk)
 end if

!the following is for the non master nodes
 call hdr_bcast(hdr,master,me,comm)
 call xmpi_bcast(nband, master,comm,ierr)
 elph_ds%nband = nband

 elph_ds%nsppol =hdr%nsppol
 elph_ds%nspinor=hdr%nspinor

!in spinor or spin polarized case, orbitals have occupation <= 1 instead of 2
 max_occ = one
 if (hdr%nspinor == 2) max_occ = half ! this accounts for the doubling of the num of bands, even though spin channels are not well defined
 if (elph_ds%nsppol > 1) max_occ = one
 write (std_out,*) ' max_occ factor  ', max_occ

 elph_ds%occ_factor = one
 if (hdr%nspinor == 1 .and. hdr%nsppol == 1) then
   elph_ds%occ_factor = one
 else if (hdr%nspinor == 2) then
   elph_ds%occ_factor = two
 else if (hdr%nsppol == 2) then
   elph_ds%occ_factor = one
 end if

!==================================================
!Read GS eigenvalues for each irreducible kpt and
!number of 1WF files contributing to the GKK file
!==================================================

 ABI_ALLOCATE(eigenGS,(nband,hdr%nkpt,elph_ds%nsppol))

 if (master == me) then
   do isppol=1,elph_ds%nsppol
     do ikpt=1,hdr%nkpt
       read(unitgkk) eigenGS(:,ikpt,isppol)
     end do
   end do

!  read number of 1WF files contributing to the GKK file
   read(unitgkk) n1wf
   write(message,'(a,i0)')' elphon : number of perturbations in the gkk file = ',n1wf
   call wrtout(std_out,message,'COLL')
 end if
 call xmpi_bcast(n1wf, master, comm, ierr)
 call xmpi_bcast(eigenGS, master, comm, ierr)

!==================================================
!Set elph_ds%fermie: either comes from anaddb input file or from wfk file
!==================================================
 elph_ds%fermie = hdr%fermie
 !elph_ds%nelect = hdr_get_nelect_byocc(Hdr)
 elph_ds%nelect = Hdr%nelect
 if (abs(anaddb_dtset%elph_fermie) > tol10) then
   elph_ds%fermie = anaddb_dtset%elph_fermie
   write(message,'(a,E20.12)')' Fermi level set by the user at :',elph_ds%fermie
   call wrtout(std_out,message,'COLL')
   Bst = ebands_from_hdr(Hdr,nband,eigenGS)
 else if (abs(anaddb_dtset%ep_extrael) > tol10) then
   if (abs(anaddb_dtset%ep_extrael) > 1.0d2) then
     write(message,'(a,E20.12)')' Doping set by the user is (negative for el doping) :',&
&     anaddb_dtset%ep_extrael
     call wrtout(std_out,message,'COLL')
     anaddb_dtset%ep_extrael = anaddb_dtset%ep_extrael*cryst%ucvol*btocm3*(-1.0d0)
   end if
   write(message,'(a,E20.12)')' Additional electrons per unit cell set by the user at :',&
&   anaddb_dtset%ep_extrael
   call wrtout(std_out,message,'COLL')
   elph_ds%nelect = elph_ds%nelect + anaddb_dtset%ep_extrael
   bst = ebands_from_hdr(Hdr,nband,eigenGS,nelect=elph_ds%nelect)

!  set Bst to use FD occupations:
   Bst%occopt = 3
!   Bst%tsmear = 0.00001_dp ! is this small etol9 Bst%tsmeatol90001_dp ! last used
   Bst%tsmear = tol9 ! is this small etol9 Bst%tsmeatol90001_dp ! last used
!  Calculate occupation numbers.
   call ebands_update_occ(Bst,-99.99_dp)
   write(message,'(a,E20.12)')' Fermi level is now calculated to be :',Bst%fermie
   call wrtout(std_out,message,'COLL')
   elph_ds%fermie = BSt%fermie
 else
   bst = ebands_from_hdr(Hdr,nband,eigenGS)
 end if
 call wrtout(std_out,message,'COLL')

!====================================================================
!Setup of the phon k-grid :
!1) get bands near Ef
!====================================================================
 call get_fs_bands(eigenGS,hdr,elph_ds%fermie,anaddb_dtset%ep_b_min, anaddb_dtset%ep_b_max,&
& elph_ds%minFSband,elph_ds%maxFSband,elph_ds%k_phon%nkptirr)

 elph_ds%nFSband = elph_ds%maxFSband - elph_ds%minFSband + 1

!Modify the band gap by sissor shift of the CB
 if (abs(anaddb_dtset%band_gap) < 10.0d0) then
   anaddb_dtset%band_gap = anaddb_dtset%band_gap*0.036749309 ! eV2Ha
   do isppol=1,elph_ds%nsppol

!First find where the gap is
     etemp_vb = 999.0d0
     top_vb = elph_ds%minFSband
     do iband = elph_ds%minFSband, elph_ds%maxFSband
       e_vb_max = maxval(eigenGS(iband,:,isppol))
       if (dabs(e_vb_max-elph_ds%fermie) < etemp_vb) then
         etemp_vb = dabs(e_vb_max-elph_ds%fermie)
         top_vb = iband
       end if
     end do
     do iband = top_vb, elph_ds%maxFSband
       e_vb_max = maxval(eigenGS(iband,:,isppol))
       if (dabs(e_vb_max-maxval(eigenGS(top_vb,:,isppol))) < tol6) then
         etemp_vb = dabs(e_vb_max-elph_ds%fermie)
         top_vb = iband
       end if
     end do
     e_vb_max = maxval(eigenGS(top_vb,:,isppol))
     e_cb_min = minval(eigenGS(top_vb+1,:,isppol))
     write(message,'(a,E20.12,2x,E20.12)')' elphon : original fermi energy = ', elph_ds%fermie
     call wrtout(std_out,message,'COLL')
     write(message,'(a,E20.12,2x,E20.12)')' elphon : top of VB, bottom of CB = ',e_vb_max, e_cb_min
     call wrtout(std_out,message,'COLL')

     do iband = top_vb+1, elph_ds%maxFSband
       eigenGS(iband,:,isppol) = eigenGS(iband,:,isppol) + (anaddb_dtset%band_gap-(e_cb_min-e_vb_max))
     end do
   end do !nsppol

!! recalculate Fermi level
   !elph_ds%nelect = hdr_get_nelect_byocc(Hdr)
   elph_ds%nelect = Hdr%nelect
   if (abs(anaddb_dtset%elph_fermie) > tol10) then
     elph_ds%fermie = anaddb_dtset%elph_fermie
     write(message,'(a,E20.12)')' Fermi level set by the user at :',elph_ds%fermie
     call wrtout(std_out,message,'COLL')
     bst = ebands_from_hdr(Hdr,nband,eigenGS)
   else if (abs(anaddb_dtset%ep_extrael) > tol10) then
     write(message,'(a,E20.12)')' Additional electrons per unit cell set by the user at :',anaddb_dtset%ep_extrael
     call wrtout(std_out,message,'COLL')
     elph_ds%nelect = elph_ds%nelect + anaddb_dtset%ep_extrael
     bst = ebands_from_hdr(Hdr,nband,eigenGS,nelect=elph_ds%nelect)

!    set Bst to use FD occupations:
     Bst%occopt = 3
!     Bst%tsmear = 0.00001_dp ! is this small etol9 Bst%tsmeatol90001_dp ! last used
     Bst%tsmear = tol9 ! is this small etol9 Bst%tsmeatol90001_dp ! last used
!    Calculate occupation numbers.
     call ebands_update_occ(Bst,-99.99_dp)
     write(message,'(a,E20.12)')' Fermi level is now calculated to be :',Bst%fermie
     call wrtout(std_out,message,'COLL')
     elph_ds%fermie = BSt%fermie
   else
     bst = ebands_from_hdr(Hdr,nband,eigenGS)
   end if
   call wrtout(std_out,message,'COLL')
 end if !modify band_gap

 if (elph_ds%ep_keepbands == 0) then !we are summing over bands
   elph_ds%ngkkband = 1
 else if (elph_ds%ep_keepbands == 1) then
!  keep the band dependency btw elph_ds%minFSband and elph_ds%maxFSband
   elph_ds%ngkkband = elph_ds%nFSband
 else
   write(message,'(a,i0)')' ep_keepbands must be 0 or 1 while it is: ',elph_ds%ep_keepbands
   MSG_BUG(message)
 end if

 write(message,'(a,i0,2x,i0)')' elphon : minFSband, maxFSband = ',elph_ds%minFSband,elph_ds%maxFSband
 call wrtout(std_out,message,'COLL')


 ABI_ALLOCATE(elph_ds%k_phon%kptirr,(3,elph_ds%k_phon%nkptirr))
 ABI_ALLOCATE(elph_ds%k_phon%irredtoGS,(elph_ds%k_phon%nkptirr))

!====================================================================
!2) order irred k-points
!====================================================================
 if (master == me) then
   call order_fs_kpts(hdr%kptns, hdr%nkpt, elph_ds%k_phon%kptirr,elph_ds%k_phon%nkptirr,elph_ds%k_phon%irredtoGS)
 end if
 call xmpi_bcast(elph_ds%k_phon%nkptirr, master, comm, ierr)
 call xmpi_bcast(elph_ds%k_phon%kptirr, master, comm, ierr)
 call xmpi_bcast(elph_ds%k_phon%irredtoGS, master, comm, ierr)

!==========================================
!3) reconstruct full kgrid from irred kpoints,
!==========================================
 call mkFSkgrid (elph_ds%k_phon, Cryst%nsym, Cryst%symrec, timrev)

! check that kptrlatt is coherent with kpt found here
 nkpt_tmp = elph_ds%kptrlatt(1,1)*elph_ds%kptrlatt(2,2)*elph_ds%kptrlatt(3,3)
 if (sum(abs(elph_ds%kptrlatt(:,:))) /= nkpt_tmp) then
   MSG_WARNING(' the input kptrlatt is not diagonal... ')
 end if
 if (anaddb_dtset%ifltransport > 1 .and. nkpt_tmp /= elph_ds%k_phon%nkpt) then
   write(message,'(a,i0,a,i0)')&
&   ' the input kptrlatt is inconsistent  ', nkpt_tmp, " /= ", elph_ds%k_phon%nkpt
   MSG_ERROR(message)
 end if

 if (anaddb_dtset%ifltransport==3 ) then
!====================================================================
! The real irred kpt, now only used by get_tau_k
!====================================================================

   ABI_ALLOCATE(indkpt1,(elph_ds%k_phon%nkpt))
   ABI_ALLOCATE(wtk_fullbz,(elph_ds%k_phon%nkpt))
   ABI_ALLOCATE(wtk_folded,(elph_ds%k_phon%nkpt))

   wtk_fullbz(:) = one/dble(elph_ds%k_phon%nkpt) !weights normalized to unity
   call symkpt(0,cryst%gmet,indkpt1,0,elph_ds%k_phon%kpt,elph_ds%k_phon%nkpt,elph_ds%k_phon%new_nkptirr,&
&   Cryst%nsym,Cryst%symrec,timrev,wtk_fullbz,wtk_folded)

   write (message,'(2a,i0)')ch10,' Number of irreducible k-points = ',elph_ds%k_phon%new_nkptirr
   call wrtout(std_out,message,'COLL')

   ABI_ALLOCATE(elph_ds%k_phon%new_kptirr,(3,elph_ds%k_phon%new_nkptirr))
   ABI_ALLOCATE(elph_ds%k_phon%new_wtkirr,(elph_ds%k_phon%new_nkptirr))
   ABI_ALLOCATE(elph_ds%k_phon%new_irredtoGS,(elph_ds%k_phon%new_nkptirr))

   ikpt_irr = 0
   do ikpt=1,elph_ds%k_phon%nkpt
     if (wtk_folded(ikpt) /= zero) then
       ikpt_irr = ikpt_irr + 1
       elph_ds%k_phon%new_kptirr(:,ikpt_irr) = elph_ds%k_phon%kpt(:,ikpt)
       elph_ds%k_phon%new_wtkirr(ikpt_irr) = wtk_folded(ikpt)
       elph_ds%k_phon%new_irredtoGS(ikpt_irr) = ikpt
     end if
   end do
   if (ikpt_irr .ne. elph_ds%k_phon%new_nkptirr) then
     write (message,'(a)')' The number of irred nkpt does not match! '
     MSG_ERROR(message)
   end if

   ABI_DEALLOCATE(indkpt1)
   ABI_DEALLOCATE(wtk_fullbz)
   ABI_DEALLOCATE(wtk_folded)
 end if

!====================================================================
!4) setup weights for integration (gaussian or tetrahedron method)
!====================================================================
 elph_ds%k_phon%nband = elph_ds%nFSband
 elph_ds%k_phon%nsppol = elph_ds%nsppol
 elph_ds%k_phon%nsym = Cryst%nsym
 ABI_ALLOCATE(elph_ds%k_phon%wtk,(elph_ds%nFSband,elph_ds%k_phon%nkpt,elph_ds%k_phon%nsppol))

 call ep_fs_weights(anaddb_dtset%ep_b_min, anaddb_dtset%ep_b_max, eigenGS, anaddb_dtset%elphsmear, &
& elph_ds%fermie, cryst%gprimd, elph_ds%k_phon%irredtoGS, elph_ds%kptrlatt, max_occ, elph_ds%minFSband, nband, elph_ds%nFSband, &
& elph_ds%nsppol, anaddb_dtset%telphint, elph_ds%k_phon)

!distribute k-points among processors, if any
 call elph_k_procs(nproc, elph_ds%k_phon)

!=====================================================
!get kpt info from the fine grid part
!=====================================================
 if (anaddb_dtset%use_k_fine == 1) then

   if (abs(anaddb_dtset%band_gap) < 10.0d0) then
     write (message,'(a)')' Not coded yet when use_k_fine and band_gap are both used'
     MSG_ERROR(message)
   end if

   if (master == me) then
     if (open_file("densergrid_GKK",message,newunit=unitfskgrid,form="unformatted",status="old") /=0) then
       MSG_ERROR(message)
     end if
     !read the header of file
     call hdr_fort_read(hdr1, unitfskgrid, fform)
     ABI_CHECK(fform/=0,'denser grid GKK header was mis-read. fform == 0')
   end if
   call hdr_bcast(hdr1,master,me,comm)

   ABI_ALLOCATE(eigenGS_fine,(nband,hdr1%nkpt,elph_ds%nsppol))

   if (master == me) then
     do isppol=1,elph_ds%nsppol
       do ikpt=1,hdr1%nkpt
         read(unitfskgrid) eigenGS_fine(:,ikpt,isppol)
       end do
     end do
     close(unitfskgrid)
   end if
   call xmpi_bcast(eigenGS_fine, master, comm, ierr)

!  Reinit the structure storing the eigevalues.
!  Be careful. This part has not been tested.
   call ebands_free(Bst)
   bst = ebands_from_hdr(hdr1,nband,eigenGS_fine)

   elph_ds%k_fine%nkptirr = hdr1%nkpt
   ABI_ALLOCATE(elph_ds%k_fine%kptirr,(3,elph_ds%k_fine%nkptirr))
   ABI_ALLOCATE(elph_ds%k_fine%irredtoGS,(elph_ds%k_fine%nkptirr))

   call order_fs_kpts(hdr1%kptns, hdr1%nkpt, elph_ds%k_fine%kptirr,&
&   elph_ds%k_fine%nkptirr,elph_ds%k_fine%irredtoGS)

   call hdr_free(hdr1)

   call mkFSkgrid (elph_ds%k_fine, Cryst%nsym, Cryst%symrec, timrev)

   elph_ds%k_fine%nband = elph_ds%nFSband
   elph_ds%k_fine%nsppol = elph_ds%nsppol
   elph_ds%k_fine%nsym = Cryst%nsym

   ABI_ALLOCATE(elph_ds%k_fine%wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt,elph_ds%nsppol))

   kptrlatt_fine = elph_ds%kptrlatt_fine

   call ep_fs_weights(anaddb_dtset%ep_b_min, anaddb_dtset%ep_b_max, &
&   eigenGS_fine, anaddb_dtset%elphsmear, &
&   elph_ds%fermie, cryst%gprimd, elph_ds%k_fine%irredtoGS, kptrlatt_fine, &
&   max_occ, elph_ds%minFSband, nband, elph_ds%nFSband, &
&   elph_ds%nsppol, anaddb_dtset%telphint, elph_ds%k_fine)

 else ! not using k_fine
   elph_ds%k_fine%nband = elph_ds%k_phon%nband
   elph_ds%k_fine%nsppol = elph_ds%k_phon%nsppol
   elph_ds%k_fine%nsym = elph_ds%k_phon%nsym

   elph_ds%k_fine%nkpt = elph_ds%k_phon%nkpt
   elph_ds%k_fine%nkptirr = elph_ds%k_phon%nkptirr

   elph_ds%k_fine%my_nkpt = elph_ds%k_phon%my_nkpt

   ABI_ALLOCATE(elph_ds%k_fine%my_kpt,(elph_ds%k_fine%nkpt))
   elph_ds%k_fine%my_kpt = elph_ds%k_phon%my_kpt

   ABI_ALLOCATE(elph_ds%k_fine%my_ikpt,(elph_ds%k_fine%my_nkpt))
   elph_ds%k_fine%my_ikpt = elph_ds%k_phon%my_ikpt

   ABI_ALLOCATE(elph_ds%k_fine%kptirr,(3,elph_ds%k_fine%nkptirr))
   elph_ds%k_fine%kptirr = elph_ds%k_phon%kptirr
   ABI_ALLOCATE(elph_ds%k_fine%wtkirr,(elph_ds%k_fine%nkptirr))
   elph_ds%k_fine%wtkirr = elph_ds%k_phon%wtkirr

   ABI_ALLOCATE(elph_ds%k_fine%wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt,elph_ds%k_fine%nsppol))
   elph_ds%k_fine%wtk = elph_ds%k_phon%wtk
   ABI_ALLOCATE(elph_ds%k_fine%kpt,(3,elph_ds%k_fine%nkpt))
   elph_ds%k_fine%kpt = elph_ds%k_phon%kpt

   call copy_kptrank(elph_ds%k_phon%kptrank_t, elph_ds%k_fine%kptrank_t)

   ABI_ALLOCATE(elph_ds%k_fine%irr2full,(elph_ds%k_fine%nkptirr))
   elph_ds%k_fine%irr2full = elph_ds%k_phon%irr2full
   ABI_ALLOCATE(elph_ds%k_fine%full2irr,(3,elph_ds%k_fine%nkpt))
   elph_ds%k_fine%full2irr = elph_ds%k_phon%full2irr
   ABI_ALLOCATE(elph_ds%k_fine%full2full,(2,elph_ds%k_fine%nsym,elph_ds%k_fine%nkpt))
   elph_ds%k_fine%full2full = elph_ds%k_phon%full2full

   ABI_ALLOCATE(elph_ds%k_fine%irredtoGS,(elph_ds%k_fine%nkptirr))
   elph_ds%k_fine%irredtoGS = elph_ds%k_phon%irredtoGS

!  call elph_k_copy(elph_ds%k_phon, elph_ds%k_fine)

   kptrlatt_fine = elph_ds%kptrlatt

   ABI_ALLOCATE(eigenGS_fine,(nband,elph_ds%k_fine%nkptirr,elph_ds%nsppol))

   eigenGS_fine = eigenGS
 end if ! k_fine or not

 if (elph_ds%kptrlatt_fine(1,1) == 0) then ! when there is not input for kptrlatt_fine
   elph_ds%kptrlatt_fine = kptrlatt_fine
 end if

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon k and q grids have been setup after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

!====================================================================
!5) calculate DOS at Ef
!====================================================================
 ABI_ALLOCATE(elph_ds%n0,(elph_ds%nsppol))

!SPPOL sum over spin channels to get total DOS
!channels decoupled => use separate values for DOS_up(Ef) resp down
 do isppol=1,elph_ds%nsppol
   elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt
 end do

 if (elph_ds%nsppol == 1) then
   write (std_out,*) ' elphon : the estimated DOS(E_Fermi) = ', elph_ds%n0(1), ' states/Ha/spin '
   write (std_out,*) ' elphon : the total FS weight and # of kpoints = ',sum(elph_ds%k_fine%wtk),elph_ds%k_fine%nkpt
 else if (elph_ds%nsppol == 2) then
   write (std_out,*) ' elphon : the spin up   DOS(E_Fermi) = ', elph_ds%n0(1), ' states/Ha/spin '
   write (std_out,*) ' elphon : the spin down DOS(E_Fermi) = ', elph_ds%n0(2), ' states/Ha/spin '
   write (std_out,*) ' elphon : total DOS(E_Fermi) = ', elph_ds%n0(1)+elph_ds%n0(2), ' states/Ha '
   write (std_out,*) ' elphon : the spin up   FS weight and # of kpoints = ',&
&   sum(elph_ds%k_fine%wtk(:,:,1)),elph_ds%k_fine%nkpt
   write (std_out,*) ' elphon : the spin down FS weight and # of kpoints = ',&
&   sum(elph_ds%k_fine%wtk(:,:,2)),elph_ds%k_fine%nkpt
 else
   write (message,'(a,i0)') 'bad value for nsppol ', elph_ds%nsppol
   MSG_ERROR(message)
 end if

 ABI_ALLOCATE(elph_ds%gkk_intweight,(elph_ds%ngkkband,elph_ds%k_phon%nkpt,elph_ds%nsppol))

 if (elph_ds%ep_keepbands == 0) then
!  use trivial integration weights  for single band,
!  since average over bands is done in normsq_gkk
   elph_ds%gkk_intweight(1,:,:) = one

 else if (elph_ds%ep_keepbands == 1) then
!  use elph_ds%k_fine%wtk since average over bands is not done in normsq_gkk
   if (elph_ds%use_k_fine == 1) then
     call d2c_weights(elph_ds)
   end if
   elph_ds%gkk_intweight(:,:,:) = elph_ds%k_phon%wtk(:,:,:)
 else
   write(message,'(a,i0)')' ep_keepbands must be 0 or 1 while it is : ',elph_ds%ep_keepbands
   MSG_ERROR(message)
 end if

 ep_prt_wtk = 0
 if (ep_prt_wtk == 1) then
   do iband=1, elph_ds%ngkkband
     do ikpt_fine=1, elph_ds%k_fine%nkpt
       write (300,*) ikpt_fine, elph_ds%gkk_intweight(iband,ikpt_fine,1)
     end do
   end do
 end if


 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon weights and DOS setup after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

!Output of the Fermi Surface
 if (anaddb_dtset%prtfsurf == 1 .and. master == me) then
   fname=trim(elph_ds%elph_base_name) // '_BXSF'
   if (ebands_write_bxsf(Bst, Cryst, fname) /= 0) then
     MSG_WARNING("Cannot produce file for Fermi surface, check log file for more info")
   end if
 end if

!=========================================================
!Get equivalence between a kpt_phon pair and a qpt in qpt_full
!only works if the qpt grid is complete (identical to
!the kpt one, with a basic shift of (0,0,0)
!=========================================================

!mapping of k + q onto k' for k and k' in full BZ
 ABI_ALLOCATE(FSfullpqtofull,(elph_ds%k_phon%nkpt,elph_ds%nqpt_full))

!qpttoqpt(itim,isym,iqpt) = qpoint index which transforms to iqpt under isym and with time reversal itim.
 ABI_ALLOCATE(qpttoqpt,(2,Cryst%nsym,elph_ds%nqpt_full))

 call wrtout(std_out,'elphon: calling mkqptequiv to set up the FS qpoint set',"COLL")

 call mkqptequiv (FSfullpqtofull,Cryst,elph_ds%k_phon%kpt,elph_ds%k_phon%nkpt,&
& elph_ds%nqpt_full,qpttoqpt,elph_ds%qpt_full)

!==========================================
!Set up dataset for phonon interpolations
!==========================================

!transfer ifltransport flag to structure
 elph_tr_ds%ifltransport=anaddb_dtset%ifltransport
!transfer name of files file for ddk
 elph_tr_ds%ddkfilename=ddkfilename

!reduce qpt_full to correct zone
 do iqpt=1,elph_ds%nqpt_full
   call wrap2_pmhalf(elph_ds%qpt_full(1,iqpt),kpt(1),res)
   call wrap2_pmhalf(elph_ds%qpt_full(2,iqpt),kpt(2),res)
   call wrap2_pmhalf(elph_ds%qpt_full(3,iqpt),kpt(3),res)
   elph_ds%qpt_full(:,iqpt)=kpt
 end do

!test density of k+q grid: the following should be close to n0 squared
!FIXME: generalize for sppol
 res = zero
 do ikpt_fine = 1, elph_ds%k_phon%nkpt
   do iqpt = 1, elph_ds%nqpt_full
     kpt = elph_ds%k_phon%kpt(:,ikpt_fine) + elph_ds%qpt_full(:,iqpt)
     call get_rank_1kpt (kpt,symrankkpt,elph_ds%k_phon%kptrank_t)
     iFSkpq = elph_ds%k_phon%kptrank_t%invrank(symrankkpt)
     do iband = 1, elph_ds%ngkkband
       do ibandp = 1, elph_ds%ngkkband
         res = res + elph_ds%gkk_intweight(iband,ikpt_fine,1)*elph_ds%gkk_intweight(ibandp,iFSkpq,1)
       end do
     end do
   end do
 end do
 res = res / elph_ds%k_phon%nkpt/elph_ds%k_phon%nkpt
 write (std_out,*) 'elphon: integrated value of intweight for given k and q grid : ', res, res / elph_ds%n0(1)**2

 res = zero
 do ikpt_fine = 1, elph_ds%k_phon%nkpt
   do iqpt = 1, elph_ds%k_phon%nkpt
     kpt = elph_ds%k_phon%kpt(:,ikpt_fine) + elph_ds%k_phon%kpt(:,iqpt)
     call get_rank_1kpt (kpt,symrankkpt,elph_ds%k_phon%kptrank_t)
     iFSkpq = elph_ds%k_phon%kptrank_t%invrank(symrankkpt)
     do iband = 1, elph_ds%ngkkband
       do ibandp = 1, elph_ds%ngkkband
         res = res + elph_ds%gkk_intweight(iband,ikpt_fine,1)*elph_ds%gkk_intweight(ibandp,iFSkpq,1)
       end do
     end do
   end do
 end do
 res = res / elph_ds%k_phon%nkpt/elph_ds%k_phon%nkpt
 write (std_out,*) 'elphon: integrated value of intweight for double k grid : ', res, res / elph_ds%n0(1)**2

!===================================================
!Allocate all important arrays for FS integrations
!===================================================

!Record sizes for matrices on disk: complex and real versions (for real and recip space resp!)
 onegkksize = 2*elph_ds%nbranch*elph_ds%nbranch*&
& elph_ds%ngkkband*elph_ds%ngkkband*&
& elph_ds%nsppol*kind(realdp_ex)

 elph_tr_ds%onegkksize=onegkksize

 write (message,'(4a)')&
& ' elphon : preliminary setup completed ',ch10,&
& '          calling get_all_gkq to read in all the e-ph matrix elements',ch10
 call wrtout(std_out,message,'COLL')

!flag to do scalar product in gkq before interpolation:
!should also used in interpolate_gkk and mkph_linwid
 if (elph_ds%ep_scalprod==0) then
   write (std_out,*) ' elphon: will NOT perform scalar product with phonon'
   write (std_out,*) '  displacement vectors in read_gkk. ep_scalprod==0'
 else if (elph_ds%ep_scalprod==1) then
   write (std_out,*) ' elphon: will perform scalar product with phonon'
   write (std_out,*) '  displacement vectors in read_gkk. ep_scalprod==1'
 else
   MSG_ERROR('illegal value for ep_scalprod')
 end if

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon begin gkq construction after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

 call get_all_gkq (elph_ds,Cryst,ifc,Bst,FSfullpqtofull,nband,n1wf,onegkksize,&
& qpttoqpt,anaddb_dtset%ep_prt_yambo,unitgkk,elph_tr_ds%ifltransport)

 if (master == me) then
   close (unitgkk)
 end if

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon end gkq construction after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

 if (elph_tr_ds%ifltransport==1 .or. elph_tr_ds%ifltransport==2 .or. elph_tr_ds%ifltransport==3)then

!  check inputs
!  TODO: should be done at earlier stage of initialization and checking
   if (elph_ds%ngkkband /= elph_ds%nFSband) then
     write (message,'(a)') 'need to keep electron band dependency in memory for transport calculations'
     MSG_ERROR(message)
   end if

!  bxu, moved the allocation from get_veloc_tr to elphon
   if (anaddb_dtset%use_k_fine == 1) then
     ABI_ALLOCATE(elph_tr_ds%el_veloc,(elph_ds%k_fine%nkpt,nband,3,elph_ds%nsppol))
   else
     ABI_ALLOCATE(elph_tr_ds%el_veloc,(elph_ds%k_phon%nkpt,nband,3,elph_ds%nsppol))
   end if
   ABI_ALLOCATE(elph_tr_ds%FSelecveloc_sq,(3,elph_ds%nsppol))

!  this only needs to be read in once - the fermi level average is later done many times with get_veloc_tr
   if (me == master) then
     if (anaddb_dtset%use_k_fine == 1) then
       call read_el_veloc(nband,elph_ds%k_fine%nkpt,elph_ds%k_fine%kpt,elph_ds%nsppol,elph_tr_ds)
     else
       call read_el_veloc(nband,elph_ds%k_phon%nkpt,elph_ds%k_phon%kpt,elph_ds%nsppol,elph_tr_ds)
     end if
   end if
   call xmpi_bcast (elph_tr_ds%el_veloc, master, comm, ierr)

   call get_veloc_tr(elph_ds,elph_tr_ds)
 end if

!Output of the Fermi velocities
!to be used for Mayavi visualization
 if (anaddb_dtset%prtfsurf == 1 .and. master == me) then
   fname = trim(elph_ds%elph_base_name) // '_VTK'

!  FIXME
!  shiftk is defined neither in the anaddb nor in the hdr data type
!  an incorrect FS will be produced in case of a shifted k-grid used during the GS calculation
!  check if we are using a unshifthed kgrid, obviously doesnt work in case
!  of multiple shifts containg a zero translation but in this case prtbxsf should work
   shiftk=one
   do ii=1,hdr%nkpt
     if (all(hdr%kptns(:,ii) == zero)) shiftk=zero
   end do

   use_afm=(hdr%nsppol==1.and.hdr%nspden==2)
!  MG FIXME warning time reversal is always assumed to be present.
!  the header should report this information.

   use_tr=(timrev==1)

   nk1 = elph_ds%kptrlatt_fine(1,1)
   nk2 = elph_ds%kptrlatt_fine(2,2)
   nk3 = elph_ds%kptrlatt_fine(3,3)

   ABI_ALLOCATE(v_surf,(nband,nk1+1,nk2+1,nk3+1,3,elph_ds%nsppol))
   v_surf = zero
   do isppol=1,elph_ds%nsppol
     do iband=1,nband
       do ikpt = 1, nk1+1
         do jkpt = 1, nk2+1
           do kkpt = 1, nk3+1
             ik1 = ikpt
             ik2 = jkpt
             ik3 = kkpt
             if (ikpt > nk1) ik1 = ikpt - nk1
             if (jkpt > nk2) ik2 = jkpt - nk2
             if (kkpt > nk3) ik3 = kkpt - nk3
             ikpt_fine = (ik1-1)*nk2*nk3 + (ik2-1)*nk3 + ik3
!            v_surf(iband,ikpt,jkpt,kkpt,:,isppol)=elph_tr_ds%el_veloc(ikpt_fine,iband,:,isppol)*elph_ds%k_fine%wtk(iband,ikpt_fine,isppol)
             v_surf(iband,ikpt,jkpt,kkpt,:,isppol)=elph_tr_ds%el_veloc(ikpt_fine,iband,:,isppol)
           end do
         end do
       end do
     end do
   end do

   call printvtk(eigenGS,v_surf,zero,elph_ds%fermie,Cryst%gprimd,&
&   elph_ds%kptrlatt_fine,nband,hdr%nkpt,hdr%kptns,&
&   Cryst%nsym,use_afm,Cryst%symrec,Cryst%symafm,use_tr,elph_ds%nsppol,shiftk,1,fname,ierr)

   ABI_DEALLOCATE(v_surf)

 end if !anaddb_dtset%prtfsurf

!============================================================================
!Evaluate lambda and omega_log using the weighted sum over the irred q-points
!found in the GKK file. All the data we need are stored in elph_ds%qgrid_data
!============================================================================

 if (master == me) then
   fname=trim(elph_ds%elph_base_name) // '_QPTS'
   call outelph(elph_ds,anaddb_dtset%enunit,fname)
 end if

!========================================================
!Get FS averaged gamma matrices and Fourier transform to real space
!========================================================

 ABI_ALLOCATE(coskr, (elph_ds%nqpt_full,Ifc%nrpt))
 ABI_ALLOCATE(sinkr, (elph_ds%nqpt_full,Ifc%nrpt))
 call ftgam_init(ifc%gprim, elph_ds%nqpt_full,Ifc%nrpt, elph_ds%qpt_full, Ifc%rpt, coskr, sinkr)

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon begin integration of gkq after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall

 call integrate_gamma(elph_ds,FSfullpqtofull)

 if (elph_ds%symgkq ==1) then
!  complete the gamma_qpt here instead of the gkk previously
   call complete_gamma(Cryst,elph_ds%nbranch,elph_ds%nsppol,elph_ds%nqptirred,elph_ds%nqpt_full,&
&   elph_ds%ep_scalprod,elph_ds%qirredtofull,qpttoqpt,elph_ds%gamma_qpt)
 end if

!Now FT to real space too
!NOTE: gprim (not gprimd) is used for all FT interpolations,
!to be consistent with the dimensions of the rpt, which come from anaddb.
 ABI_ALLOCATE(elph_ds%gamma_rpt, (2,elph_ds%nbranch**2,elph_ds%nsppol,Ifc%nrpt))
 elph_ds%gamma_rpt = zero

 qtor = 1 ! q --> r
 do isppol=1,elph_ds%nsppol
   call ftgam(Ifc%wghatm,elph_ds%gamma_qpt(:,:,isppol,:),elph_ds%gamma_rpt(:,:,isppol,:),natom,&
&   elph_ds%nqpt_full,Ifc%nrpt,qtor, coskr, sinkr)
 end do

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-elphon end integration and completion of gkq after tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(std_out,message,'COLL')
 tcpui = tcpu
 twalli = twall


!==========================================================
!calculate transport matrix elements, integrated over FS
!==========================================================

 if (elph_tr_ds%ifltransport == 1)then ! LOVA

   call integrate_gamma_tr_lova(elph_ds,FSfullpqtofull,elph_tr_ds)

   call complete_gamma_tr(cryst,elph_ds%ep_scalprod,elph_ds%nbranch,elph_ds%nqptirred,&
&   elph_ds%nqpt_full,elph_ds%nsppol,elph_tr_ds%gamma_qpt_trout,elph_ds%qirredtofull,qpttoqpt)

   call complete_gamma_tr(cryst,elph_ds%ep_scalprod,elph_ds%nbranch,elph_ds%nqptirred,&
&   elph_ds%nqpt_full,elph_ds%nsppol,elph_tr_ds%gamma_qpt_trin,elph_ds%qirredtofull,qpttoqpt)

   ABI_ALLOCATE(elph_tr_ds%gamma_rpt_trout,(2,9,elph_ds%nbranch**2,elph_ds%nsppol,Ifc%nrpt))
   elph_tr_ds%gamma_rpt_trout = zero

   ABI_ALLOCATE(elph_tr_ds%gamma_rpt_trin,(2,9,elph_ds%nbranch**2,elph_ds%nsppol,Ifc%nrpt))
   elph_tr_ds%gamma_rpt_trin = zero

!  Now FT to real space too
   qtor = 1 ! q --> r
   do isppol=1,elph_ds%nsppol
     do idir=1,9
       call ftgam(Ifc%wghatm,elph_tr_ds%gamma_qpt_trout(:,idir,:,isppol,:),&
&       elph_tr_ds%gamma_rpt_trout(:,idir,:,isppol,:),natom,&
&       elph_ds%nqpt_full,Ifc%nrpt,qtor, coskr, sinkr)

       call ftgam(Ifc%wghatm,elph_tr_ds%gamma_qpt_trin(:,idir,:,isppol,:),&
&       elph_tr_ds%gamma_rpt_trin(:,idir,:,isppol,:),natom,&
&       elph_ds%nqpt_full,Ifc%nrpt,qtor, coskr, sinkr)
     end do
   end do

 else if (elph_tr_ds%ifltransport==2) then ! non-LOVA case

!  Get Ef, DOS(Ef), veloc(Ef) for looped temperatures
   call get_nv_fs_temp(elph_ds,BSt,eigenGS_fine,cryst%gprimd,max_occ,elph_tr_ds)

!  Get DOS(E), veloc(E) for looped energy levels
   call get_nv_fs_en(cryst,ifc,elph_ds,eigenGS_fine,max_occ,elph_tr_ds,omega_max)

!  Save the E, N(E), v^2(E), dE
   if (master == me) then
     fname = trim(elph_ds%elph_base_name) // '_EPTS'
     if (open_file(fname,message,newunit=unit_epts,status="unknown") /=0) then
       MSG_ERROR(message)
     end if
     do isppol = 1, elph_ds%nsppol
       write(unit_epts,"(a,i6)") '# E, N(E), v^2(E), dE for spin channel ', isppol
       do ie1 = 1, elph_ds%nenergy
         write(unit_epts,"(4E20.12)") elph_tr_ds%en_all(isppol,ie1), elph_tr_ds%dos_n(ie1,isppol),&
&         elph_tr_ds%veloc_sq(1,isppol,ie1), elph_tr_ds%de_all(isppol,ie1)
       end do
     end do
     close(unit=unit_epts)
   end if

   ABI_ALLOCATE(tmp_veloc_sq1,(3,elph_ds%nsppol))
   ABI_ALLOCATE(tmp_veloc_sq2,(3,elph_ds%nsppol))
   ABI_ALLOCATE(elph_tr_ds%tmp_gkk_intweight1,(elph_ds%ngkkband,elph_ds%k_phon%nkpt,elph_ds%nsppol))
   ABI_ALLOCATE(elph_tr_ds%tmp_gkk_intweight2,(elph_ds%ngkkband,elph_ds%k_phon%nkpt,elph_ds%nsppol))
   ABI_ALLOCATE(elph_tr_ds%tmp_velocwtk1,(elph_ds%ngkkband,elph_ds%k_phon%nkpt,3,elph_ds%nsppol))
   ABI_ALLOCATE(elph_tr_ds%tmp_velocwtk2,(elph_ds%ngkkband,elph_ds%k_phon%nkpt,3,elph_ds%nsppol))
   ABI_ALLOCATE(elph_tr_ds%tmp_vvelocwtk1,(elph_ds%ngkkband,elph_ds%k_phon%nkpt,3,3,elph_ds%nsppol))
   ABI_ALLOCATE(elph_tr_ds%tmp_vvelocwtk2,(elph_ds%ngkkband,elph_ds%k_phon%nkpt,3,3,elph_ds%nsppol))

   tmp_veloc_sq1 = zero
   tmp_veloc_sq2 = zero
   elph_tr_ds%tmp_gkk_intweight1 = zero
   elph_tr_ds%tmp_gkk_intweight2 = zero
   elph_tr_ds%tmp_velocwtk1 = zero
   elph_tr_ds%tmp_velocwtk2 = zero
   elph_tr_ds%tmp_vvelocwtk1 = zero
   elph_tr_ds%tmp_vvelocwtk2 = zero

   if (elph_ds%ep_lova .eq. 1) then
     tmp_nenergy = 1
   else if (elph_ds%ep_lova .eq. 0) then
     tmp_nenergy = elph_ds%nenergy
   else
     write(message,'(a,i0)')' ep_lova must be 0 or 1 while it is : ', elph_ds%ep_lova
     MSG_ERROR(message)
   end if

!  This only works for ONE temperature!! for test only
   elph_ds%n0(:) = elph_tr_ds%dos_n0(1,:)

!  bxu, no need for complete sets of ie1 and ie2
!  Only save those within the range of omega_max from Ef
   ABI_ALLOCATE(pair2red,(tmp_nenergy,tmp_nenergy))
   pair2red = zero

   elph_ds%n_pair = 0
   do ie1=1,tmp_nenergy
     e1 = elph_tr_ds%en_all(1,ie1)
     e2 = e1 - omega_max
     if (e2 .lt. elph_tr_ds%en_all(1,1)) then
       i_start = 1
     else
       i_start = 1
       diff = dabs(e2-elph_tr_ds%en_all(1,1))
       do ie2 = 2, tmp_nenergy
         if (dabs(e2-elph_tr_ds%en_all(1,ie2)) .lt. diff) then
           diff = dabs(e2-elph_tr_ds%en_all(1,ie2))
           i_start = ie2
         end if
       end do
     end if
     e2 = e1 + omega_max
     if (e2 .gt. elph_tr_ds%en_all(1,tmp_nenergy)) then
       i_end = tmp_nenergy
     else
       i_end = 1
       diff = dabs(e2-elph_tr_ds%en_all(1,1))
       do ie2 = 2, tmp_nenergy
         if (dabs(e2-elph_tr_ds%en_all(1,ie2)) .lt. diff) then
           diff = dabs(e2-elph_tr_ds%en_all(1,ie2))
           i_end = ie2
         end if
       end do
     end if
     do ie2 = i_start, i_end
       elph_ds%n_pair = elph_ds%n_pair + 1
       pair2red(ie1,ie2) = elph_ds%n_pair
     end do
   end do

!  symmetrize paire2red
   elph_ds%n_pair = 0
   do ie1 = 1, tmp_nenergy
     do ie2 = 1, tmp_nenergy
       if (pair2red(ie1,ie2) .ne. 0 .or. pair2red(ie2,ie1) .ne. 0) then
         elph_ds%n_pair = elph_ds%n_pair + 1
         pair2red(ie1,ie2) = elph_ds%n_pair
       end if
     end do
   end do

   write(message,'(a,i3,a)')' There are  ', elph_ds%n_pair, '  energy pairs. '
   call wrtout(std_out,message,'COLL')

   ABI_ALLOCATE(red2pair,(2,elph_ds%n_pair))
   red2pair = zero
   elph_ds%n_pair = 0
   do ie1 = 1, tmp_nenergy
     do ie2 = 1, tmp_nenergy
       if (pair2red(ie1,ie2) .ne. 0 .or. pair2red(ie2,ie1) .ne. 0) then
         elph_ds%n_pair = elph_ds%n_pair + 1
         red2pair(1,elph_ds%n_pair) = ie1
         red2pair(2,elph_ds%n_pair) = ie2
       end if
     end do
   end do

!  moved from integrate_gamma_tr to here
   ABI_ALLOCATE(elph_tr_ds%gamma_qpt_tr,(2,9,elph_ds%nbranch**2,elph_ds%nsppol,elph_ds%nqpt_full))
   ABI_ALLOCATE(elph_tr_ds%gamma_rpt_tr,(2,9,elph_ds%nbranch**2,elph_ds%nsppol,Ifc%nrpt,4,elph_ds%n_pair))
   elph_tr_ds%gamma_rpt_tr = zero

   s1ofssp = (/1,1,-1,-1/)
   s2ofssp = (/1,-1,1,-1/)

!  Get gamma
   do ie=1,elph_ds%n_pair
     ie1 = red2pair(1,ie)
     ie2 = red2pair(2,ie)

     tmp_veloc_sq1(:,:)=elph_tr_ds%veloc_sq(:,:,ie1)
     elph_tr_ds%tmp_gkk_intweight1(:,:,:) = elph_tr_ds%tmp_gkk_intweight(:,:,:,ie1)
     elph_tr_ds%tmp_velocwtk1(:,:,:,:) = elph_tr_ds%tmp_velocwtk(:,:,:,:,ie1)
     elph_tr_ds%tmp_vvelocwtk1(:,:,:,:,:) = elph_tr_ds%tmp_vvelocwtk(:,:,:,:,:,ie1)

     tmp_veloc_sq2(:,:)=elph_tr_ds%veloc_sq(:,:,ie2)
     elph_tr_ds%tmp_gkk_intweight2(:,:,:) = elph_tr_ds%tmp_gkk_intweight(:,:,:,ie2)
     elph_tr_ds%tmp_velocwtk2(:,:,:,:) = elph_tr_ds%tmp_velocwtk(:,:,:,:,ie2)
     elph_tr_ds%tmp_vvelocwtk2(:,:,:,:,:) = elph_tr_ds%tmp_vvelocwtk(:,:,:,:,:,ie2)

     do ssp=1,4  ! (s,s'=+/-1, condense the indices)
       s1=s1ofssp(ssp)
       s2=s2ofssp(ssp)
       elph_tr_ds%gamma_qpt_tr = zero

       call integrate_gamma_tr(elph_ds,FSfullpqtofull,s1,s2, &
&       tmp_veloc_sq1,tmp_veloc_sq2,elph_tr_ds)

       call complete_gamma_tr(cryst,elph_ds%ep_scalprod,elph_ds%nbranch,elph_ds%nqptirred,&
&       elph_ds%nqpt_full,elph_ds%nsppol,elph_tr_ds%gamma_qpt_tr,elph_ds%qirredtofull,qpttoqpt)

!      Now FT to real space too
       qtor = 1 ! q --> r
       do isppol=1,elph_ds%nsppol
         do idir=1,9
           call ftgam(Ifc%wghatm,elph_tr_ds%gamma_qpt_tr(:,idir,:,isppol,:),&
&           elph_tr_ds%gamma_rpt_tr(:,idir,:,isppol,:,ssp,ie),natom,&
&           elph_ds%nqpt_full,Ifc%nrpt,qtor,coskr, sinkr)
         end do
       end do

     end do !ss
   end do !ie

   ABI_DEALLOCATE(tmp_veloc_sq1)
   ABI_DEALLOCATE(tmp_veloc_sq2)
 end if ! ifltransport

 ABI_DEALLOCATE(qpttoqpt)
 ABI_DEALLOCATE(FSfullpqtofull)


!==============================================================
!Calculate phonon linewidths, interpolating on chosen qpoints
!==============================================================

 call mkph_linwid(Cryst,ifc,elph_ds,anaddb_dtset%nqpath,anaddb_dtset%qpath)

!==============================================================
!the nesting factor calculation
!FIXME: this could go higher up, before the call to get_all_gkq
!you only need the kpt and weight info
!==============================================================
 if (any(anaddb_dtset%prtnest==[1,2])) then

   nestname = trim(elph_ds%elph_base_name) // "_NEST"
   call mknesting(elph_ds%k_phon%nkpt,elph_ds%k_phon%kpt,elph_ds%kptrlatt,elph_ds%nFSband,&
&   elph_ds%k_phon%wtk,anaddb_dtset%nqpath,anaddb_dtset%qpath,elph_ds%nqpt_full, &
&   elph_ds%qpt_full,nestname,cryst%gprimd,cryst%gmet,anaddb_dtset%prtnest,qptrlatt)
 end if

!======================================================
!Calculate alpha^2 F integrating over fine kpt_phon grid
!======================================================

 ABI_ALLOCATE(a2f_1d,(elph_ds%na2f))
 ABI_ALLOCATE(dos_phon,(elph_ds%na2f))

 call mka2f(Cryst,Ifc,a2f_1d,dos_phon,elph_ds,elph_ds%kptrlatt_fine,elph_ds%mustar)

!calculate transport spectral function and coefficients
 if (elph_tr_ds%ifltransport==1 )then ! LOVA

   call mka2f_tr_lova(cryst,ifc,elph_ds,elph_ds%ntemper,elph_ds%tempermin,elph_ds%temperinc,elph_tr_ds)

 else if (elph_tr_ds%ifltransport==2 )then ! non LOVA

   call mka2f_tr(cryst,ifc,elph_ds,elph_ds%ntemper,elph_ds%tempermin,elph_ds%temperinc,pair2red,elph_tr_ds)

   ABI_DEALLOCATE(pair2red)
   ABI_DEALLOCATE(red2pair)

 else if (elph_tr_ds%ifltransport==3 )then ! get k-dependent tau

   call get_tau_k(Cryst,ifc,Bst,elph_ds,elph_tr_ds,eigenGS,max_occ)
   !call trans_rta(elph_ds,elph_tr_ds,cryst%gprimd,eigenGS,max_occ,cryst%ucvol)
 end if ! ifltransport

 ABI_DEALLOCATE(eigenGS)
 ABI_DEALLOCATE(eigenGS_fine)


!evaluate a2F only using the input Q-grid (without using interpolated matrices)
!SCOPE: test the validity of the Fourier interpolation
 call wrtout(std_out,' elphon : calling mka2fQgrid',"COLL")

 fname=trim(elph_ds%elph_base_name) // '_A2F_QGRID'
 call mka2fQgrid(elph_ds,fname)

!=============================================
!Eliashberg equation in 1-D (isotropic case)
!=============================================

 call eliashberg_1d(a2f_1d,elph_ds,anaddb_dtset%mustar)

 ABI_DEALLOCATE(a2f_1d)
 ABI_DEALLOCATE(dos_phon)

!MJV: 20070805 should exit here. None of the rest is tested or used yet to my knowledge

!========================================================================
!Now gkk contains the matrix elements of dH(1)/dxi i=1,2,3
!for kpoints on the FS but qpoints only in the given grid {Q}.
!
!1.) Need to complete the gkk elements for q and k\prime=k+q not
!in the set of {k+Q} by Fourier interpolation on the Q.
!
!2.) Need to complete the dynamical matrices and phonon freqs for
!all q between points on the FS.
!
!3.) With the eigenvectors e_ph of the dyn mats, do the scalar product
!e_ph . gkk, which implies the gkk are turned to the eigenbasis of
!the phonons. Before the (non eigen-) modes are ordered
!atom1 xred1 atom1 xred2 atom1 xred3
!atom2 xred1 atom2 xred2 atom2 xred3 ...
!=======================================================================

 make_gkk2=.false.

 if (.not. make_gkk2) then
   call wrtout(std_out,' elphon : skipping full g(k,k") interpolation ',"COLL")
 else

!  ==========================================================
!  FT of recip space gkk matrices to real space (gkk_rpt)
!  NOTE: could be made into FFT, couldnt it? If shifts are
!  used with a homogeneous grid
!  ==========================================================
   write (message,'(2a,i0)')ch10,&
&   ' elphon : Fourier transform (q --> r) of the gkk matrices using nrpt = ',Ifc%nrpt
   call wrtout(std_out,message,'COLL')

   call get_all_gkr(elph_ds,ifc%gprim,natom,Ifc%nrpt,onegkksize,Ifc%rpt,elph_ds%qpt_full,Ifc%wghatm)

!  =========================================================
!  complete gkk2 for all qpts between points
!  on full kpt grid (interpolation from real space values)
!  =========================================================

   write(message,'(2a)')ch10,&
&   ' elphon : Calling get_all_gkk2 to calculate gkk2 for q points over the full k grid'
   call wrtout(std_out,message,'COLL')

   call get_all_gkk2(cryst,ifc,elph_ds,elph_ds%k_phon%kptirr,elph_ds%k_phon%kpt)
 end if

!=====================================================
!Here should be the anisotropic Eliashberg equations.
!=====================================================

!clean and deallocate junk
 call ebands_free(Bst)
 call elph_ds_clean(elph_ds)
 call elph_tr_ds_clean(elph_tr_ds)
 call hdr_free(hdr)

 ABI_DEALLOCATE(coskr)
 ABI_DEALLOCATE(sinkr)

 if (is_open(elph_ds%unitgkq)) close(elph_ds%unitgkq)

end subroutine elphon
!!***

!!****f* m_elphon/outelph
!! NAME
!! outelph
!!
!! FUNCTION
!!  Output to stdout and file the data for electron phonon coupling,
!!  on the q-points which were really calculated by abinit (no interpolation yet)
!!
!! INPUTS
!!  elph_ds  the elph_type structured variable
!!  enunit   from the anaddb dataset 0 ==> Hartree and cm-1;
!!                                   1 ==> meV and Thz;
!!
!! OUTPUT
!!  only write
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      bfactor,destroy_kptrank,mkkptrank,wrtout
!!
!! SOURCE

subroutine outelph(elph_ds,enunit,fname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'outelph'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: enunit
 character(len=fnlen),intent(in) :: fname
 type(elph_type),intent(in) :: elph_ds

!Local variables-------------------------------
!scalars
 integer :: ibranch,ii,iqfull,iqirr,isppol,jj,nfile,qmax,qnest_max,qnest_min
 integer :: nbranch,nsppol,nqptirred
 real(dp) :: lambda_q_max,lambda_qbranch_max,lambda_tot,nest_max,nest_min
 real(dp) :: omegalog_q,omegalog_qgrid,tc_macmill
 character(len=500) :: msg
 type(kptrank_type) :: kptrank_t
!arrays
 integer :: qbranch_max(2)
 real(dp),allocatable :: lambda_q(:,:),nestfactor(:),qirred(:,:)

! *************************************************************************

 if ( ALL (enunit /= (/0,1,2/)) )  then
   write(msg,'(a,i0)')' enunit should be 0 or 1 or 2 while it is ',enunit
   MSG_BUG(msg)
 end if

 nbranch   = elph_ds%nbranch
 nsppol    = elph_ds%nsppol
 nqptirred = elph_ds%nqptirred

!==========================================================
!write header
!==========================================================
 if (open_file(fname,msg,newunit=nfile,form="formatted",status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if

 write(msg,'(2a,80a,4a,80a)')ch10,' ',('=',ii=1,80),ch10,&
& ' Values of the parameters that define the electron-phonon calculation',ch10,&
& ' ',('=',ii=1,80)
 call wrtout(nfile,msg,'COLL')

 write(msg,'(a,i10,a,i10,a,i10)')&
& ' nkpt_phon    = ',elph_ds%k_phon%nkpt,   ' nkpt_phonirred = ',elph_ds%k_phon%nkptirr,&
& ' nqpt      = ',elph_ds%nqpt_full
 call wrtout(nfile,msg,'COLL')

 if (nsppol==1) then
   write(msg,'(2a,f10.7,a,f10.6,a,f10.7)')ch10,&
&   ' Fermi DOS = ',elph_ds%n0(1),       ' Fermi level = ',elph_ds%fermie,&
&   ' mustar    = ',elph_ds%mustar
   call wrtout(nfile,msg,'COLL')
 else if (nsppol==2) then
   write(msg,'(2a,f10.7,f10.7,a,f10.6,a,f10.7)')ch10,&
&   ' Fermi DOS (up/dn) = ',elph_ds%n0(1),elph_ds%n0(2),       ' Fermi level = ',elph_ds%fermie,&
&   ' mustar    = ',elph_ds%mustar
   call wrtout(nfile,msg,'COLL')
 else
   MSG_BUG("bad value for nsppol")
 end if

 write(msg,'(2a,i10,a,i10,a,i10)')ch10,&
& ' minFSband = ',elph_ds%minFSband,' maxFSband   = ',elph_ds%maxFSband,&
& ' ngkkband  = ',elph_ds%ngkkband
 call wrtout(nfile,msg,'COLL')

 write(msg,'(80a,a)')('=',ii=1,80),ch10
 call wrtout(nfile,msg,'COLL')

!==========================================================
!evaluate lambda and omega_log as a weighted sum over the q grid
!NOTE: in this part of the code atomic units are used
!==========================================================

 ABI_ALLOCATE(lambda_q,(nqptirred,nsppol))
 lambda_q=zero
 lambda_tot=zero ; lambda_q_max=zero
 qmax=0          ; lambda_qbranch_max=zero
 qbranch_max(:)=1; omegalog_qgrid=zero

 do iqirr=1,nqptirred
   omegalog_q=zero

   do isppol=1,nsppol
     do ibranch=1,nbranch
!      find Max lambda(q,n)
       if (elph_ds%qgrid_data(iqirr,ibranch,isppol,3) > lambda_qbranch_max) then
         lambda_qbranch_max=elph_ds%qgrid_data(iqirr,ibranch,isppol,3)
         qbranch_max(1)=iqirr
         qbranch_max(2)=ibranch
       end if
       lambda_q(iqirr,isppol)=lambda_q(iqirr,isppol)+elph_ds%qgrid_data(iqirr,ibranch,isppol,3)
       if (abs(elph_ds%qgrid_data(iqirr,ibranch,isppol,1)) <= tol10) cycle
       omegalog_q=omegalog_q + elph_ds%qgrid_data(iqirr,ibranch,isppol,3)*log(abs(elph_ds%qgrid_data(iqirr,ibranch,isppol,1)))
     end do

     lambda_tot=lambda_tot+elph_ds%wtq(elph_ds%qirredtofull(iqirr))*lambda_q(iqirr,isppol)
     omegalog_qgrid=omegalog_qgrid+elph_ds%wtq(elph_ds%qirredtofull(iqirr))*omegalog_q


!    find Max lambda(q)
     if (lambda_q(iqirr,isppol) > lambda_q_max) then
       lambda_q_max=lambda_q(iqirr,isppol)
       qmax=iqirr
     end if
   end do

 end do !iqirr

 omegalog_qgrid=exp(omegalog_qgrid/lambda_tot)

 write (msg,'(3a,2(a,es16.8))')                                                                              &
& ' Values of Lambda, Omega_log and Tc obtained using the weighted sum over the input Q-grid',ch10,ch10,&
& ' Isotropic Lambda = ',lambda_tot,'  Input mustar     = ',elph_ds%mustar
 call wrtout(nfile,msg,'COLL')

 if (enunit==0) then !use hartree and cm-1
   write (msg,'(2a,es16.8,a,es16.8,a)')ch10,&
&   ' Omega_log        = ',omegalog_qgrid,' (Ha) ',omegalog_qgrid*Ha_cmm1,' (cm-1)'
   call wrtout(nfile,msg,'COLL')
 else if (enunit==1) then !mev Thz
   write (msg,'(2a,es16.8,a,es16.8,a)')ch10,&
&   ' Omega_log        = ',omegalog_qgrid*Ha_eV/1000._dp,' (meV) ',omegalog_qgrid*Ha_THz,' (THz)'
   call wrtout(nfile,msg,'COLL')
 else !hartree,cm-1,mev,Thz,kelvin
   write (msg,'(2a,es16.8,a,es16.8,3a,es16.8,a,es16.8,3a,es16.8,a)')ch10,                              &
&   ' Omega_log        = ',omegalog_qgrid,' (Ha)  ',omegalog_qgrid*Ha_cmm1,' (cm-1)',ch10,             &
&   '                  = ',omegalog_qgrid*Ha_eV/1000._dp,' (meV) ',omegalog_qgrid*Ha_THz,' (THz)',ch10,&
&   '                  = ',omegalog_qgrid*Ha_K,' (K) '
   call wrtout(nfile,msg,'COLL')
 end if

 tc_macmill = omegalog_qgrid/1.2_dp&
& *exp((-1.04_dp*(one+lambda_tot)) / (lambda_tot-elph_ds%mustar*(one+0.62_dp*lambda_tot)))

 if (enunit==0) then !use hartree and cm-1
   write (msg,'(2a,es16.8,a,es16.8,2a)')ch10,&
&   ' MacMillan Tc     = ',tc_macmill,' (Ha) ',tc_macmill*Ha_cmm1,' (cm-1) ',ch10
   call wrtout(nfile,msg,'COLL')
 else if (enunit==1) then !use mev and Thz
   write (msg,'(2a,es16.8,a,es16.8,2a)')ch10,&
&   ' MacMillan Tc     = ',tc_macmill*Ha_eV/1000._dp,' (meV) ',tc_macmill*Ha_THz,' (THz) ',ch10
   call wrtout(nfile,msg,'COLL')
 else !use hartree,cm-1,mev,Thz,kelvin
   write (msg,'(2a,es16.8,a,es16.8,3a,es16.8,a,es16.8,3a,es16.8,2a)')ch10,                 &
&   ' MacMillan Tc     = ',tc_macmill,' (Ha)  ',tc_macmill*Ha_cmm1,' (cm-1) ',ch10,            &
&   '                  = ',tc_macmill*Ha_eV/1000._dp,' (meV) ',tc_macmill*Ha_THz,' (THz) ',ch10,&
&   '                  = ',tc_macmill*Ha_K,' (K) ',ch10
   call wrtout(nfile,msg,'COLL')
 end if

!==========================================================
!output lambda(q) values for each q point in the irred grid
!==========================================================

 write(msg,'(2a)')' Irreducible q-points and corresponding Lambda(q)',ch10
 call wrtout(nfile,msg,'COLL')

 do isppol=1,nsppol
   write(msg,'(a,i3,2a)')'  === isppol ', isppol,' === ',ch10
   call wrtout(nfile,msg,'COLL')
!
   do iqirr=1,nqptirred
     iqfull=elph_ds%qirredtofull(iqirr)
     write(msg,'(i5,a,3(es16.8,1x),a,es16.8,a)')&
&     iqfull,') ',elph_ds%qpt_full(:,iqfull),'(',lambda_q(iqirr,isppol),'  )'
     call wrtout(nfile,msg,'COLL')
   end do
!
 end do

!use same indexing as that used for the full q-grid
 qmax=elph_ds%qirredtofull(qmax)
 qbranch_max(1)=elph_ds%qirredtofull(qbranch_max(1))

 write (msg,'(2a,es16.8,a,i6,3a,es16.8,a,i6,a,i4)')ch10,            &
& ' Max lambda(q)      = ',lambda_q_max,      ' at qpt ',qmax,')',ch10, &
& ' Max lambda(q,n)    = ',lambda_qbranch_max,' at qpt ',qbranch_max(1),&
& ') and Mode number ',qbranch_max(2)
 call wrtout(nfile,msg,'COLL')

!==========================================================
!evaluation of the nesting-factor over the irreducible q grid.
!==========================================================

!fill irreducile q-grid
 ABI_ALLOCATE(qirred,(3,nqptirred))
 qirred(:,:)=zero

 do iqirr=1,nqptirred
   qirred(:,iqirr)=elph_ds%qpt_full(:,elph_ds%qirredtofull(iqirr))
 end do

 call mkkptrank (elph_ds%k_phon%kpt,elph_ds%k_phon%nkpt,kptrank_t)

 ABI_ALLOCATE(nestfactor,(nqptirred))

!NOTE: weights are not normalised, the normalisation factor in reintroduced in bfactor
 call bfactor(elph_ds%k_phon%nkpt,elph_ds%k_phon%kpt,nqptirred,qirred,kptrank_t,&
& elph_ds%k_phon%nkpt,elph_ds%k_phon%wtk,elph_ds%nFSband,nestfactor)

 ABI_DEALLOCATE(qirred)
 call destroy_kptrank (kptrank_t)


!find Max and min of the nesting factor
!NOTE maxloc and minloc are arrays so they cannot be used in the formatted output
!anyway the size of nestfactor is not so huge!!!
 nest_max=maxval(nestfactor); nest_min=minval(nestfactor)

 qnest_max=0
 do iqirr=1,nqptirred
   if (nestfactor(iqirr)==nest_max) then
     qnest_max=iqirr
     exit
   end if
 end do

 qnest_min=0
 do iqirr=1,nqptirred
   if (nestfactor(iqirr)==nest_min) then
     qnest_min=iqirr
     exit
   end if
 end do

 write (std_out,*) maxloc(nestfactor),minloc(nestfactor)
 write(msg,'(a,(a,es16.8,a,i6,a),a,(a,es16.8,a,i6,a))')ch10,  &
& ' Max nesting factor = ',nest_max,' at qpt ',qnest_max,') ',ch10,&
& ' min nesting factor = ',nest_min,' at qpt ',qnest_min,') '
 call wrtout(nfile,msg,'COLL')

!==========================================================
!Write ph-linewidths and lambda(q,n) obtained before the
!Fourier interpolation
!==========================================================

 write (msg,'(2a)')ch10,&
& ' Phonon frequencies, linewidths and e-ph coefficients for each irreducible q point '
 call wrtout(nfile,msg,'COLL')

 do isppol=1,nsppol
   write (msg,'(a,i3,a)') '========= quantities for isppol = ', isppol, ' ================='
   call wrtout(nfile,msg,'COLL')
   do iqirr=1,nqptirred
!    same numbering as that used for irred q points
     iqfull=elph_ds%qirredtofull(iqirr)
!    write(std_out,*) 'iqfull = ', iqfull
     write(msg,'(64a,i6,a,3(es16.8),3a,es16.8,a,es16.8,2a,es16.8,a,f8.3,65a)')ch10,&
&     ' ',('=',jj=1,60),ch10,&
&     ' qpt ',iqfull,') ',elph_ds%qpt_full(:,iqfull),ch10,ch10,&
&     ' Weight    = ',elph_ds%wtq(iqfull),'    Lambda(q,isppol) = ',lambda_q(iqirr,isppol),ch10,&
&     ' Nest fact = ',nestfactor(iqirr),'    (',100*nestfactor(iqirr)/nest_max,' % of max_value )',ch10,&
&     ' ',('=',jj=1,60),ch10,' Mode number    Frequency       Linewidth        Lambda(q,n)'
     call wrtout(nfile,msg,'COLL')

!    use units according to enunit
     if (enunit==0 .or. enunit==2) then !hartree and cm-1
       write(msg,'(63a)')' ',('-',jj=1,60),ch10,&
       '                  (Ha)             (Ha)'
       call wrtout(nfile,msg,'COLL')
       do ibranch=1,nbranch
!        branch index, frequency, linewidth, lamda(q,n) (hartree units)
         write(msg,'(i6,5x,3(es16.8,1x))' )ibranch,(elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,3)
         call wrtout(nfile,msg,'COLL')
       end do
       write(msg,'(63a)')' ',('-',jj=1,60),ch10,&
&       '                 (cm-1)           (cm-1)'
       call wrtout(nfile,msg,'COLL')
       do ibranch=1,nbranch
!        branch index, frequency, linewidth (in cm-1)
         write(msg,'(i6,5x,2(es16.8,1x))' )ibranch,(Ha_cmm1*elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,2)
         call wrtout(nfile,msg,'COLL')
       end do
     end if !hartree and cm-1

     if (enunit==2 .or. enunit==1) then !write also meV Thz and Kelvin
       write(msg,'(63a)')' ',('-',jj=1,60),ch10,&
&       '                 (meV)             (meV)'
       call wrtout(nfile,msg,'COLL')
       if (enunit == 1 ) then !write also lambda values
         do ibranch=1,nbranch
!          branch index, frequency, linewidth, lamda(q,n) (mev units)
           write(msg,'(i6,5x,3(es16.8,1x))' )ibranch,((Ha_eV/1000._dp)*elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,2),&
&           elph_ds%qgrid_data(iqirr,ibranch,isppol,3)
           call wrtout(nfile,msg,'COLL')
         end do
       else !do not write lambda values
         do ibranch=1,nbranch
!          branch index, frequency, linewidth (in meV)
           write(msg,'(i6,5x,2(es16.8,1x))' )ibranch,((Ha_eV/1000._dp)*elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,2)
           call wrtout(nfile,msg,'COLL')
         end do
       end if

       write(msg,'(63a)')' ',('-',jj=1,60),ch10,&
&       '                 (Thz)             (Thz)'
       call wrtout(nfile,msg,'COLL')
       do ibranch=1,nbranch
!        branch index, frequency, linewidth (in Thz)
         write(msg,'(i6,5x,2(es16.8,1x))' )ibranch,(Ha_THz*elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,2)
         call wrtout(nfile,msg,'COLL')
       end do

       if (enunit == 2 ) then !kelvin
         write(msg,'(63a)')' ',('-',jj=1,60),ch10,&
&         '                  (K)               (K)'
         call wrtout(nfile,msg,'COLL')
         do ibranch=1,nbranch
!          branch index, frequency, linewidth (in Kelvin)
           write(msg,'(i6,5x,2(es16.8,1x))' )ibranch,(Ha_K*elph_ds%qgrid_data(iqirr,ibranch,isppol,jj),jj=1,2)
           call wrtout(nfile,msg,'COLL')
         end do
       end if !kelvin

     end if  !end write also meV Thz and Kelvin

     write(msg,'(62a)')' ',('=',jj=1,60),ch10
     call wrtout(nfile,msg,'COLL')

   end do !nqptirred
 end do !nsppol

 ABI_DEALLOCATE(nestfactor)
 ABI_DEALLOCATE(lambda_q)

 close (nfile)

end subroutine outelph
!!***

!!****f* m_elphon/rchkGSheader
!!
!! NAME
!! rchkGSheader
!!
!! FUNCTION
!! This routine reads the GS header information in the GKK file and checks it
!!
!! INPUTS
!!  natom = number of atoms from DDB, for check
!!  kptirr_phon = coordinates of the irreducible kpoints close to the FS
!!
!! OUTPUT
!!  hdr = header information
!!  nband = number of bands for rest of calculation
!!          should be the same for all kpts
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      hdr_echo,hdr_fort_read
!!
!! SOURCE

subroutine rchkGSheader (hdr,natom,nband,unitgkk)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'rchkGSheader'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,unitgkk
 integer,intent(out) :: nband
 type(hdr_type),intent(inout) :: hdr

!Local variables-------------------------------
!scalars
 integer :: fform
 character(len=500) :: message

! *************************************************************************
!
!read in general header of _GKK file
!this is where we get nkpt, ngkpt(:,:)... which are also read in
!rdddb9 and inprep8. Probably should do some checking to avoid
!using ddb files from other configurations
!
 rewind(unitgkk)
 call hdr_fort_read(hdr, unitgkk, fform)
 ABI_CHECK(fform/=0," GKK header mis-read. fform == 0")

 if (hdr%natom /= natom) then
   MSG_ERROR('natom in gkk file is different from anaddb input')
 end if

 if (any(hdr%nband(:) /= hdr%nband(1))) then
   write(message,'(3a)')&
&   'Use the same number of bands for all kpts: ',ch10,&
&   'could have spurious effects if efermi is too close to the last band '
   MSG_ERROR(message)
 end if

 call hdr_echo(hdr, fform, 4, unit=std_out)

 nband=hdr%nband(1)

end subroutine rchkGSheader
!!***

!!****f* m_elphon/mkfskgrid
!!
!! NAME
!! mkfskgrid
!!
!! FUNCTION
!! This routine sets up the full FS kpt grid by symmetry
!!
!! INPUTS
!!  nsym    = number of symmetries for the full system
!!  symrec  = reciprocal space symmetries (those for the kpts)
!!  timrev  = 1 if time reversal symmetry is to be used
!!
!! OUTPUT
!!  elph_k datastructure:
!!  elph_k%nkpt           = full number of kpoints close to the FS
!!  elph_k%kpt            = full set of kpoints close to the FS
!!  elph_k%wtkirr         = weights of the irreducible kpoints
!!  elph_k%kphon_irr2full = indices of irred kpoints in full array
!!
!! NOTES
!!  WARNING: supposes kpt grid has full symmetry!! Not always true!!!
!!    but should be for Monkhorst-Pack, efficient grids.
!!    otherwise you get an error message in interpolate_gkk because
!!    an FS kpt can not be found in the gkk file.
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      destroy_kptrank,get_rank_1kpt,mkkptrank,sort_int,wrap2_pmhalf,wrtout
!!
!! SOURCE

subroutine mkFSkgrid (elph_k, nsym, symrec, timrev)

 use m_sort

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkFSkgrid'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym,timrev
 type(elph_kgrid_type),intent(inout) :: elph_k
!arrays
 integer,intent(in) :: symrec(3,3,nsym)

!Local variables-------------------------------
!scalars
 integer :: ikpt1,ikpt2,isym,itim,new,symrankkpt
 real(dp) :: timsign, res
 character(len=500) :: message

!arrays
 real(dp) :: kpt(3),redkpt(3)
 integer, allocatable :: sortindexing(:), rankallk(:)

 integer, allocatable :: tmpkphon_full2irr(:,:)
 real(dp), allocatable :: tmpkpt(:,:)

! *************************************************************************

 if(timrev /= 1 .and. timrev /= 0)then
   write (message,'(a,i0)')' timrev must be 1 or 0 but found timrev= ',timrev
   MSG_BUG(message)
 end if

 ABI_ALLOCATE(tmpkphon_full2irr,(3,2*elph_k%nkptirr*nsym))
 tmpkphon_full2irr = -1

 ABI_ALLOCATE(tmpkpt,(3,2*elph_k%nkptirr*nsym))

 ABI_ALLOCATE(elph_k%wtkirr,(elph_k%nkptirr))
 elph_k%wtkirr(:) = zero

!first allocation for irred kpoints - will be destroyed below
 call mkkptrank (elph_k%kptirr,elph_k%nkptirr,elph_k%kptrank_t)
 ABI_ALLOCATE(rankallk,(elph_k%kptrank_t%max_rank))

!elph_k%kptrank_t%invrank is used as a placeholder in the following loop
 rankallk = -1
 elph_k%kptrank_t%invrank = -1

!replicate all irred kpts by symmetry to get the full k grid.
 elph_k%nkpt=0 !zero k-points found so far
 do isym=1,nsym
   do itim=0,1
     timsign = one-two*itim
     do ikpt1=1,elph_k%nkptirr
!      generate symmetrics of kpt ikpt1
       kpt(:) = timsign*(symrec(:,1,isym)*elph_k%kptirr(1,ikpt1) + &
&       symrec(:,2,isym)*elph_k%kptirr(2,ikpt1) + &
&       symrec(:,3,isym)*elph_k%kptirr(3,ikpt1))

       call get_rank_1kpt (kpt,symrankkpt,elph_k%kptrank_t)

!      is the kpt on the full grid (may have lower symmetry than full spgroup)
!      is kpt among the full FS kpts found already?
       if (elph_k%kptrank_t%invrank(symrankkpt) == -1) then
         elph_k%wtkirr(ikpt1)=elph_k%wtkirr(ikpt1)+1
         elph_k%nkpt=elph_k%nkpt+1

         call wrap2_pmhalf(kpt(1),redkpt(1),res)
         call wrap2_pmhalf(kpt(2),redkpt(2),res)
         call wrap2_pmhalf(kpt(3),redkpt(3),res)
         tmpkpt(:,elph_k%nkpt) = redkpt
         tmpkphon_full2irr(1,elph_k%nkpt) = ikpt1
!        save sym that sends irred kpt ikpt1 onto full kpt
         tmpkphon_full2irr(2,elph_k%nkpt) = isym
         tmpkphon_full2irr(3,elph_k%nkpt) = itim

         elph_k%kptrank_t%invrank(symrankkpt) = elph_k%nkpt
         rankallk(elph_k%nkpt) = symrankkpt
       end if

     end do !end loop over irred k points
   end do !end loop over timrev
 end do !end loop over symmetry

 write(message,'(a,i0)')'mkfskgrid: after first evaluation, elph_k%nkpt= ', elph_k%nkpt
 call wrtout(std_out,message,"COLL")

 elph_k%wtkirr(:) = elph_k%wtkirr(:) / elph_k%nkpt

!copy the kpoints and full --> irred kpt map
!reorder the kpts to get rank increasing monotonically with a sort
!also reorder tmpkphon_full2irr
 ABI_ALLOCATE(elph_k%kpt,(3,elph_k%nkpt))
 ABI_ALLOCATE(elph_k%full2irr,(3,elph_k%nkpt))
 ABI_ALLOCATE(sortindexing,(elph_k%nkpt))

 do ikpt1=1,elph_k%nkpt
   sortindexing(ikpt1)=ikpt1
 end do
 call sort_int(elph_k%nkpt, rankallk, sortindexing)
 do ikpt1=1,elph_k%nkpt
   if (sortindexing(ikpt1) < 1 .or. sortindexing(ikpt1) > elph_k%nkpt) then
     MSG_BUG('sorted k ranks are out of bounds: 1 to nkpt')
   end if
   elph_k%kpt(:,ikpt1) = tmpkpt(:,sortindexing(ikpt1))
   elph_k%full2irr(:,ikpt1) = tmpkphon_full2irr(:,sortindexing(ikpt1))
 end do

 ABI_DEALLOCATE(sortindexing)
 ABI_DEALLOCATE(rankallk)
 ABI_DEALLOCATE(tmpkphon_full2irr)
 ABI_DEALLOCATE(tmpkpt)
 call destroy_kptrank (elph_k%kptrank_t)


!make proper full rank arrays
 call mkkptrank (elph_k%kpt,elph_k%nkpt,elph_k%kptrank_t)


!find correspondence table between irred FS kpoints and a full one
 ABI_ALLOCATE(elph_k%irr2full,(elph_k%nkptirr))
 elph_k%irr2full(:) = 0

 do ikpt1=1,elph_k%nkptirr
   call get_rank_1kpt (elph_k%kptirr(:,ikpt1),symrankkpt,elph_k%kptrank_t)
   elph_k%irr2full(ikpt1) = elph_k%kptrank_t%invrank(symrankkpt)
 end do

!find correspondence table between FS kpoints under symmetry
 ABI_ALLOCATE(elph_k%full2full,(2,nsym,elph_k%nkpt))
 elph_k%full2full(:,:,:) = -999

 do ikpt1=1,elph_k%nkpt
!  generate symmetrics of kpt ikpt1
   do isym=1,nsym
     do itim=0,timrev
       timsign = one-two*itim
       kpt(:) = timsign*(symrec(:,1,isym)*elph_k%kpt(1,ikpt1) + &
&       symrec(:,2,isym)*elph_k%kpt(2,ikpt1) + &
&       symrec(:,3,isym)*elph_k%kpt(3,ikpt1))

!      which kpt is it among the full FS kpts
       call get_rank_1kpt (kpt,symrankkpt,elph_k%kptrank_t)
       ikpt2 = elph_k%kptrank_t%invrank(symrankkpt)
       new=1
       if (ikpt2 /= -1) then
         elph_k%full2full(itim+1,isym,ikpt2) = ikpt1
         new = 0
       end if

       if (new == 1) then
         write(std_out,*) ' mkfskgrid Error: FS kpt ',ikpt1,' has no symmetric under sym', isym,' with itim ',itim
         write(std_out,*) ' redkpt = ', redkpt
         write(std_out,*) ' symrankkpt,ikpt2 = ', symrankkpt,ikpt2
         MSG_ERROR("Fatal error, cannot continue")
       end if
     end do
   end do
 end do

!got nkpt, tmpkpt, kphon_full2irr, kphon_full2full, and wtkirr

end subroutine mkFSkgrid
!!***

!!****f* m_elphon/mka2f
!!
!! NAME
!! mka2f
!!
!! FUNCTION
!!  calculate the FS averaged alpha^2F function
!!
!! INPUTS
!! Cryst<crystal_t>=data type gathering info on the crystalline structure.
!! Ifc<ifc_type>=Object containing the interatomic force constants.
!!  elph_ds
!!    elph_ds%gkk2 = gkk2 matrix elements on full FS grid for each phonon mode
!!    elph_ds%nbranch = number of phonon branches = 3*natom
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%k_phon%nkpt = number of kpts included in the FS integration
!!    elph_ds%k_phon%kpt = coordinates of all FS kpoints
!!    elph_ds%k_phon%wtk = integration weights on the FS
!!    elph_ds%n0 = DOS at the Fermi level calculated from the k_phon integration weights (event. 2 spin pol)
!!  mustar = coulomb pseudopotential parameter
!!  natom = number of atoms
!!
!! OUTPUT
!!  a2f_1d = 1D alpha
!!  dos_phon = density of states for phonons
!!  elph_ds
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      d2c_wtq,ep_ph_weights,ftgam,ftgam_init,gam_mult_displ,ifc_fourq
!!      phdispl_cart2red,simpson_int,wrtout,zgemm
!!
!! NOTES
!!   copied from ftiaf9.f
!!
!! SOURCE

subroutine mka2f(Cryst,ifc,a2f_1d,dos_phon,elph_ds,kptrlatt,mustar)

 use m_special_funcs,  only : fermi_dirac, bose_einstein

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mka2f'
 use interfaces_14_hidewrite
 use interfaces_77_ddb, except_this_one => mka2f
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: mustar
 type(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: Cryst
 type(elph_type),target,intent(inout) :: elph_ds
!arrays
 integer, intent(in) :: kptrlatt(3,3)
 real(dp),intent(out) :: a2f_1d(elph_ds%na2f),dos_phon(elph_ds%na2f)

!Local variables -------------------------
!scalars
 integer :: natom,iFSqpt,ibranch,iomega,nbranch,na2f,nsppol,nkpt,nrpt
 integer :: isppol,jbranch,unit_a2f,unit_phdos,ep_scalprod
 integer :: itemp, ntemp = 100
 real(dp) :: temp
 real(dp) :: a2fprefactor,avgelphg,avglambda,avgomlog,diagerr
 real(dp) :: lambda_2,lambda_3,lambda_4,lambda_5
 real(dp) :: spinfact
 real(dp) :: lambda_iso(elph_ds%nsppol)
 real(dp) :: lqn,omega
 real(dp) :: omegalog(elph_ds%nsppol)
 real(dp) :: omlog_qn
 real(dp) :: tc_macmill,a2fsmear,domega,omega_min,omega_max
 real(dp) :: gaussval, gaussprefactor, gaussfactor, gaussmaxval, xx
 character(len=500) :: msg
 character(len=fnlen) :: fname,base_name
!arrays
 real(dp) :: displ_cart(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: eigval(elph_ds%nbranch)
 real(dp) :: gam_now(2,elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: imeigval(elph_ds%nbranch)
! real(dp) :: pheigvec(2*elph_ds%nbranch*elph_ds%nbranch),phfrq(elph_ds%nbranch)
 real(dp) :: tmp_a2f(elph_ds%na2f)
 real(dp) :: tmp_gam1(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_gam2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_phondos(elph_ds%na2f),n0(elph_ds%nsppol)
 real(dp),pointer :: kpt(:,:)
 real(dp),allocatable :: phfrq(:,:)
 real(dp),allocatable :: pheigvec(:,:)
 real(dp),allocatable :: tmp_wtq(:,:,:)
 real(dp),allocatable :: a2f1mom(:),a2f2mom(:),a2f3mom(:),a2f4mom(:)
 real(dp),allocatable :: a2f_1mom(:),a2flogmom(:)
 real(dp),allocatable :: a2flogmom_int(:)
 real(dp),allocatable :: coskr(:,:)
 real(dp),allocatable :: sinkr(:,:)
 real(dp),allocatable :: linewidth_of_t(:)
 real(dp),allocatable :: linewidth_integrand(:,:)

! *********************************************************************
!calculate a2f for frequencies between 0 and elph_ds%omega_max

 DBG_ENTER("COLL")

!might need kptrlatt for finer interpolation later
 ABI_UNUSED(kptrlatt(1,1))

 ! nrpt = number of real-space points for FT interpolation
 nrpt = Ifc%nrpt
 natom = Cryst%natom

 nbranch   =  elph_ds%nbranch
 na2f      =  elph_ds%na2f
 nsppol    =  elph_ds%nsppol
 base_name =  elph_ds%elph_base_name
 a2fsmear  =  elph_ds%a2fsmear
 nkpt      =  elph_ds%k_phon%nkpt
 kpt       => elph_ds%k_phon%kpt

 ep_scalprod = elph_ds%ep_scalprod
 n0        = elph_ds%n0

!spinfact should be 1 for a normal non sppol calculation without spinorbit
!for spinors it should also be 1 as bands are twice as numerous but n0 has been divided by 2
!for sppol 2 it should be 0.5 as we have 2 spin channels to sum
 spinfact = one/elph_ds%nsppol !/elph_ds%nspinor

!maximum value of frequency (a grid has to be chosen for the representation of alpha^2 F)
!WARNING! supposes this value has been set in mkelph_linwid.
 domega = (elph_ds%omega_max-elph_ds%omega_min)/(na2f-one)
 elph_ds%domega  = domega  ! MG Why do we need to store domega in elph_ds?
 omega_min       = elph_ds%omega_min
 omega_max       = elph_ds%omega_max

 gaussprefactor = sqrt(piinv) / a2fsmear
 gaussfactor = one / a2fsmear
 gaussmaxval = sqrt(-log(1.d-100))

 ! only open the file for the first sppol
 fname = trim(base_name) // '_A2F'
 if (open_file(fname,msg,newunit=unit_a2f,status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if

 !write (std_out,*) ' a2f function integrated over the FS'

!output the a2f_1d header
 write (unit_a2f,'(a)')                 '#'
 write (unit_a2f,'(a)')                 '# ABINIT package : a2f file'
 write (unit_a2f,'(a)')                 '#'
 write (unit_a2f,'(a)')                 '# a2f function integrated over the FS. omega in a.u.'
 write (unit_a2f,'(a,I10)')             '#  number of kpoints integrated over : ',nkpt
 write (unit_a2f,'(a,I10)')             '#  number of energy points : ',na2f
 write (unit_a2f,'(a,E16.6,a,E16.6,a)') '#  between omega_min = ',omega_min,' Ha and omega_max = ',omega_max,' Ha'
 write (unit_a2f,'(a,E16.6)')           '#  and the smearing width for gaussians is ',a2fsmear

 ! Open file for PH DOS
 fname = trim(base_name) // '_PDS'
 if (open_file(fname,msg,newunit=unit_phdos,status="replace") /= 0) then
   MSG_ERROR(msg)
 end if

 ! output the phonon DOS header
 write (unit_phdos,'(a)')                '#'
 write (unit_phdos,'(a)')                '# ABINIT package : phonon DOS file'
 write (unit_phdos,'(a)')                '#'
 write (unit_phdos,'(a)')                '# Phonon DOS integrated over the FS. omega in a.u. EXPERIMENTAL!!!'
 write (unit_phdos,'(a,I10)')            '# number of kpoints integrated over : ',nkpt
 write (unit_phdos,'(a,I10)')            '# number of energy points : ',na2f
 write (unit_phdos,'(a,E16.6,a,E16.6,a)')'# between omega_min = ',omega_min,' Ha and omega_max = ',omega_max,' Ha'
 write (unit_phdos,'(a,i4,a,E16.6)')     '# The DOS at Fermi level for spin ', 1, ' is ', n0(1)
 if (nsppol==2) then
   write (unit_phdos,'(a,i4,a,E16.6)')   '# The DOS at Fermi level for spin ', 2, ' is ', n0(2)
 end if
 write (unit_phdos,'(a,E16.6)')          '# and the smearing width for gaussians is ',a2fsmear
 write (unit_phdos,'(a)') '#'

!Get the integration weights, using tetrahedron method or gaussian
 ABI_ALLOCATE(tmp_wtq,(nbranch,elph_ds%k_fine%nkpt,na2f+1))
 ABI_ALLOCATE(elph_ds%k_fine%wtq,(nbranch,elph_ds%k_fine%nkpt,na2f))
 ABI_ALLOCATE(elph_ds%k_phon%wtq,(nbranch,nkpt,na2f))

 ABI_ALLOCATE(phfrq,(nbranch,elph_ds%k_fine%nkpt))
 ABI_ALLOCATE(pheigvec,(2*nbranch*nbranch,elph_ds%k_fine%nkpt))

 do iFSqpt=1,elph_ds%k_fine%nkpt
   call ifc_fourq(ifc,cryst,elph_ds%k_fine%kpt(:,iFSqpt),phfrq(:,iFSqpt),displ_cart,out_eigvec=pheigvec(:,iFSqpt))
 end do

 omega_min = omega_min - domega

 call ep_ph_weights(phfrq,elph_ds%a2fsmear,omega_min,omega_max,na2f+1,Cryst%gprimd,elph_ds%kptrlatt_fine, &
& elph_ds%nbranch,elph_ds%telphint,elph_ds%k_fine,tmp_wtq)
!call ep_ph_weights(phfrq,elph_ds%a2fsmear,omega_min,omega_max,na2f+1,Cryst%gprimd,elph_ds%kptrlatt_fine, &
!& elph_ds%nbranch,1,elph_ds%k_fine,tmp_wtq)
 omega_min = omega_min + domega

 do iomega = 1, na2f
   elph_ds%k_fine%wtq(:,:,iomega) = tmp_wtq(:,:,iomega+1)
 end do
 ABI_DEALLOCATE(tmp_wtq)

 if (elph_ds%use_k_fine == 1) then
   call d2c_wtq(elph_ds)
 end if

 ABI_ALLOCATE(coskr, (nkpt,nrpt))
 ABI_ALLOCATE(sinkr, (nkpt,nrpt))
 call ftgam_init(Ifc%gprim, nkpt, nrpt, kpt, Ifc%rpt, coskr, sinkr)

 ABI_DEALLOCATE(phfrq)
 ABI_DEALLOCATE(pheigvec)

 do isppol=1,nsppol
   write (std_out,*) '##############################################'
   write (std_out,*) 'mka2f : Treating spin polarization ', isppol
   write (std_out,*) '##############################################'

!  Average of electron phonon coupling over the whole BZ
   avgelphg = zero
!  MG20060607 Do the same for lambda and omega_log
   avglambda = zero
   avgomlog = zero

   a2f_1d(:) = zero
   dos_phon(:) = zero

!  reduce the dimenstion from fine to phon for phfrq and pheigvec
   ABI_ALLOCATE(phfrq,(nbranch,elph_ds%k_phon%nkpt))
   ABI_ALLOCATE(pheigvec,(2*nbranch*nbranch,elph_ds%k_phon%nkpt))

!  loop over qpoint in full kpt grid (presumably dense)
!  MG TODO : This loop can be performed using the IBZ and appropriated weights.
   do iFSqpt=1,nkpt
!
!    This reduced version of ftgkk supposes the kpoints have been integrated
!    in integrate_gamma. Do FT from real-space gamma grid to 1 qpt.

     if (elph_ds%ep_int_gkk == 1) then
       gam_now(:,:) = elph_ds%gamma_qpt(:,:,isppol,iFSqpt)
     else
       call ftgam(Ifc%wghatm,gam_now,elph_ds%gamma_rpt(:,:,isppol,:),natom,1,nrpt,0, &
&       coskr(iFSqpt,:), sinkr(iFSqpt,:))
     end if

     call ifc_fourq(ifc,cryst,kpt(:,iFSqpt),phfrq(:,iFSqpt),displ_cart,out_eigvec=pheigvec)

!    Diagonalize gamma matrix at qpoint (complex matrix).

!    if ep_scalprod==0 we have to dot in the displacement vectors here
     if (ep_scalprod==0) then

       call phdispl_cart2red(natom,Cryst%gprimd,displ_cart,displ_red)

       tmp_gam2 = reshape (gam_now, (/2,nbranch,nbranch/))
       call gam_mult_displ(nbranch, displ_red, tmp_gam2, tmp_gam1)

       do jbranch=1,nbranch
         eigval(jbranch) = tmp_gam1(1, jbranch, jbranch)
         imeigval(jbranch) = tmp_gam1(2, jbranch, jbranch)

         if (abs(imeigval(jbranch)) > tol8) then
           write (msg,'(a,i0,a,es16.8)')" imaginary values  branch = ",jbranch,' imeigval = ',imeigval(jbranch)
           MSG_WARNING(msg)
         end if

       end do

!      if ep_scalprod==1 we have to diagonalize the matrix we interpolated.
     else if (ep_scalprod == 1) then

!      MJV NOTE : gam_now is being recast as a (3*natom)**2 matrix here
       call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, cone, gam_now, 3*natom,&
&       pheigvec, 3*natom, czero, tmp_gam1, 3*natom)

       call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, cone, pheigvec, 3*natom,&
&       tmp_gam1, 3*natom, czero, tmp_gam2, 3*natom)

       diagerr = zero
       do ibranch=1,nbranch
         eigval(ibranch) = tmp_gam2(1,ibranch,ibranch)
         do jbranch=1,ibranch-1
           diagerr = diagerr + abs(tmp_gam2(1,jbranch,ibranch))
         end do
         do jbranch=ibranch+1,nbranch
           diagerr = diagerr + abs(tmp_gam2(1,jbranch,ibranch))
         end do
       end do

       if (diagerr > tol12) then
         write(msg,'(a,es15.8)') 'mka2f: residual in diagonalization of gamma with phon eigenvectors: ', diagerr
         MSG_WARNING(msg)
       end if

     else
       write (msg,'(a,i0)')' Wrong value for ep_scalprod = ',ep_scalprod
       MSG_BUG(msg)
     end if

!    MG20060603MG
!    there was a bug in the calculation of the phonon DOS
!    since frequencies with small e-ph interaction were skipped inside the loop
!    In this new version all the frequencies (both positive and negative) are taken into account.
!    IDEA: it could be useful to calculate the PH-dos and the a2f
!    using several smearing values to perform a convergence study
!    Now the case ep_scalprod=1 is treated in the right way although it is not default anymore
!    FIXME to be checked
!    ENDMG

!    Add all contributions from the phonon modes at this qpoint to a2f and the phonon dos.
     do ibranch=1,nbranch

!      if (abs(phfrq(ibranch,iFSqpt)) < tol10) then
       if (abs(phfrq(ibranch,iFSqpt)) < tol7) then
         a2fprefactor= zero
         lqn         = zero
         omlog_qn    = zero
       else
         a2fprefactor = eigval(ibranch)/(two_pi*abs(phfrq(ibranch,iFSqpt))*n0(isppol))
         lqn          = eigval(ibranch)/(pi*phfrq(ibranch,iFSqpt)**2*n0(isppol))
         omlog_qn     = lqn*log(abs(phfrq(ibranch,iFSqpt)))
       end if

!      Add contribution to average elphon coupling
!      MANY ISSUES WITH FINITE T SUMS. THIS IS DEFINITELY
!      NOT A CORRECT FORMULATION YET.

!      Added avglambda and avgomglog to calculate lamda and omega_log using the sum over the kpt-grid.
!      If the k-grid is dense enough, these values should be better than the corresponding quantities
!      evaluated through the integration over omega that depends on the a2fsmear

       avgelphg = avgelphg + eigval(ibranch)
       avglambda = avglambda + lqn
       avgomlog= avgomlog + omlog_qn
!      ENDMG

       omega = omega_min
       tmp_a2f(:) = zero
       tmp_phondos(:) = zero
       do iomega=1,na2f
         xx = (omega-phfrq(ibranch,iFSqpt))*gaussfactor
         omega = omega + domega
         if (abs(xx) > gaussmaxval) cycle

         gaussval = gaussprefactor*exp(-xx*xx)
         tmp_a2f(iomega) = tmp_a2f(iomega) + gaussval*a2fprefactor
         tmp_phondos(iomega) = tmp_phondos(iomega) + gaussval
       end do

!      tmp_a2f(:) = zero
!      tmp_phondos(:) = zero
!      do iomega=1,na2f
!      tmp_a2f(iomega) = tmp_a2f(iomega) + a2fprefactor*elph_ds%k_phon%wtq(ibranch,iFSqpt,iomega)
!      tmp_phondos(iomega) = tmp_phondos(iomega) + elph_ds%k_phon%wtq(ibranch,iFSqpt,iomega)
!      end do

       a2f_1d(:) = a2f_1d(:) + tmp_a2f(:)
       dos_phon(:) = dos_phon(:) + tmp_phondos(:)

     end do ! ibranch
   end do  ! iFSqpt do


!  second 1 / nkpt factor for the integration weights
   a2f_1d(:) = a2f_1d(:) / nkpt
   dos_phon(:) = dos_phon(:) / nkpt

!  MG
   avglambda = avglambda/nkpt
   avgomlog= avgomlog/nkpt
   avgomlog = exp (avgomlog/avglambda)
   write(std_out,*) ' from mka2f: for spin ', isppol
   write(std_out,*) ' w/o interpolation lambda = ',avglambda,' omega_log= ',avgomlog
!  ENDMG

   write (std_out,'(a,I4,a,E16.6)') '# The DOS at Fermi level for spin ',isppol,' is ',n0(isppol)

   write (unit_a2f,'(a,I4,a,E16.6)') '# The DOS at Fermi level for spin ',isppol,' is ',n0(isppol)
   write (unit_a2f,'(a)') '#'

   omega = omega_min
   do iomega=1,na2f
     write (unit_a2f,*) omega, a2f_1d(iomega)
     omega=omega + domega
   end do
   write (unit_a2f,*)
!
!  output the phonon DOS, but only for the first sppol case
   if (isppol == 1) then
     omega = omega_min
     do iomega=1,na2f
       write (unit_phdos,*) omega, dos_phon(iomega)
       omega=omega + domega
     end do
   end if
!
!  Do isotropic calculation of lambda and output lambda, Tc(MacMillan)
!
   ABI_ALLOCATE(a2f_1mom,(na2f))
   ABI_ALLOCATE(a2f1mom,(na2f))
   ABI_ALLOCATE(a2f2mom,(na2f))
   ABI_ALLOCATE(a2f3mom,(na2f))
   ABI_ALLOCATE(a2f4mom,(na2f))
   ABI_ALLOCATE(linewidth_integrand,(na2f,ntemp))
   ABI_ALLOCATE(linewidth_of_t,(ntemp))

   a2f_1mom=zero
   a2f1mom=zero;  a2f2mom=zero
   a2f3mom=zero;  a2f4mom=zero
   linewidth_integrand = zero

   omega = omega_min
   do iomega=1,na2f
     if (abs(omega) > tol10) then
       a2f_1mom(iomega) =    two*spinfact*a2f_1d(iomega)/abs(omega)   ! first inverse moment of alpha2F
       a2f1mom(iomega)  =    two*spinfact*a2f_1d(iomega)*abs(omega)   ! first positive moment of alpha2F
       a2f2mom(iomega)  =     a2f1mom(iomega)*abs(omega)  ! second positive moment of alpha2F
       a2f3mom(iomega)  =     a2f2mom(iomega)*abs(omega)  ! third positive moment of alpha2F
       a2f4mom(iomega)  =     a2f3mom(iomega)*abs(omega)  ! fourth positive moment of alpha2F
!
!  electron lifetimes eq 4.48 in Grimvall electron phonon coupling in Metals (with T dependency). Also 5.69-5.72, 5.125, section 3.4
!  phonon lifetimes eq 19 in Savrasov PhysRevB.54.16487 (T=0)
!  a first T dependent expression in Allen PRB 6 2577 eq 10. Not sure about the units though
!
       do itemp = 1, ntemp
         temp = (itemp-1)*10._dp*kb_HaK
         linewidth_integrand(iomega, itemp) = a2f_1d(iomega) * (fermi_dirac(omega,zero,temp) + bose_einstein(omega,temp))
       end do
     end if
     omega=omega + domega
   end do
!
!  From Allen PRL 59 1460
!  \lambda <\omega^n> = 2 \int_0^{\infty} d\omega [\alpha^2F / \omega] \omega^n
!
   lambda_iso(isppol) = simpson(domega,a2f_1mom)
   lambda_2 = simpson(domega,a2f1mom)
   lambda_3 = simpson(domega,a2f2mom)
   lambda_4 = simpson(domega,a2f3mom)
   lambda_5 = simpson(domega,a2f4mom)
   do itemp = 1, ntemp
     linewidth_of_t(itemp) = simpson(domega,linewidth_integrand(:,itemp))
! print out gamma(T) here
     temp = (itemp-1)*10._dp*kb_HaK
     write (std_out,*) 'mka2f: T, average linewidth', temp, linewidth_of_t(itemp)
   end do


   ABI_DEALLOCATE(phfrq)
   ABI_DEALLOCATE(pheigvec)
   ABI_DEALLOCATE(a2f_1mom)
   ABI_DEALLOCATE(a2f1mom)
   ABI_DEALLOCATE(a2f2mom)
   ABI_DEALLOCATE(a2f3mom)
   ABI_DEALLOCATE(a2f4mom)
   ABI_DEALLOCATE(linewidth_integrand)
   ABI_DEALLOCATE(linewidth_of_t)

   write (std_out,*) 'mka2f: elphon coupling lambdas for spin = ', isppol
   write (std_out,*) 'mka2f: isotropic lambda', lambda_iso(isppol)
   write (std_out,*) 'mka2f: positive moments of alpha2F:'
   write (std_out,*) 'lambda <omega^2> = ', lambda_2
   write (std_out,*) 'lambda <omega^3> = ', lambda_3
   write (std_out,*) 'lambda <omega^4> = ', lambda_4
   write (std_out,*) 'lambda <omega^5> = ', lambda_5
!
!  Get log moment of alpha^2F
   ABI_ALLOCATE(a2flogmom,(na2f))
   ABI_ALLOCATE(a2flogmom_int,(na2f))
   omega = omega_min
   a2flogmom(:) = zero
   do iomega=1,na2f
     if (abs(omega) > tol10) then
       a2flogmom(iomega) = a2f_1d(iomega)*log(abs(omega))/abs(omega)
     end if
     omega=omega + domega
   end do
   call simpson_int(na2f,domega,a2flogmom,a2flogmom_int)

!  NOTE: omegalog actually stores the log moment of a2F, which is the quantity to sum over spins, instead of
!  exp(moment/lambda) which is an actual frequency
   omegalog(isppol) = two*spinfact*a2flogmom_int(na2f)

   ABI_DEALLOCATE(a2flogmom)
   ABI_DEALLOCATE(a2flogmom_int)

   if (nsppol > 1) then
     write (msg, '(3a)' ) ch10,&
&     ' Warning : some of the following quantities should be integrated over spin', ch10
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   end if

   write (msg, '(3a)' ) ch10,&
&   ' Superconductivity : isotropic evaluation of parameters from electron-phonon coupling.',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   if (elph_ds%nsppol > 1) then
     write (msg, '(a,i6,a,es16.6)' )' mka2f: isotropic lambda for spin ', isppol, ' = ', lambda_iso(isppol)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   end if

   write (msg, '(a,es16.6)' )' mka2f: lambda <omega^2> = ', lambda_2
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write (msg, '(a,es16.6)' )' mka2f: lambda <omega^3> = ', lambda_3
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write (msg, '(a,es16.6)' )' mka2f: lambda <omega^4> = ', lambda_4
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write (msg, '(a,es16.6)' )' mka2f: lambda <omega^5> = ', lambda_5
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   if (elph_ds%nsppol > 1) then
     write (msg, '(a,i6,a,es16.6,a,es16.6,a)' )' mka2f: omegalog for spin ', isppol, ' = ',&
&     exp(omegalog(isppol)/lambda_iso(isppol)), ' (Ha) ', exp(omegalog(isppol)/lambda_iso(isppol))/kb_HaK, ' (Kelvin) '
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   end if

 end do ! isppol



!also print out spin-summed quantities
 lambda_2 = sum(lambda_iso(1:elph_ds%nsppol))
 write (msg, '(a,es16.6)' )' mka2f: isotropic lambda = ', lambda_2
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 omega = exp( sum(omegalog(1:elph_ds%nsppol))/lambda_2 )
 write (msg, '(a,es16.6,a,es16.6,a)' )' mka2f: omegalog  = ', omega, ' (Ha) ', omega/kb_HaK, ' (Kelvin) '
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 write (msg, '(a,es16.6)' )' mka2f: input mustar = ', mustar
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 tc_macmill = omega/1.2_dp * exp((-1.04_dp*(one+lambda_2)) / (lambda_2-mustar*(one+0.62_dp*lambda_2)))
 write ( msg, '(a,es16.6,a,es16.6,a)')'-mka2f: MacMillan Tc = ', tc_macmill, ' (Ha) ', tc_macmill/kb_HaK, ' (Kelvin) '
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 close(unit=unit_a2f)
 close(unit=unit_phdos)

 ABI_DEALLOCATE(elph_ds%k_fine%wtq)
 ABI_DEALLOCATE(elph_ds%k_phon%wtq)

 ABI_DEALLOCATE(coskr)
 ABI_DEALLOCATE(sinkr)

 DBG_EXIT("COLL")

end subroutine mka2f
!!***

!!****f* m_elphon/mka2fQgrid
!! NAME
!! mka2fQgrid
!!
!! FUNCTION
!!  Calculate the Eliashberg function only using the phonon linewidths evaluated
!!  in the irreducible q-points of the coarse q-grid.
!!  The obtained results are useful to check the validity of the Fourier interpolation
!!
!! INPUTS
!!  elph_ds = electron-phonon dataset
!!  nunit = integer number for the output file
!!
!! OUTPUT
!!  Only write
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      simpson_int,wrtout
!!
!! SOURCE

subroutine mka2fQgrid(elph_ds,fname)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mka2fQgrid'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: fname
 type(elph_type),intent(in) :: elph_ds

!Local variables -------------------------
!scalars
 integer :: ibranch,iomega,iost,ismear,isppol,nsmear,nunit,qptirred
 real(dp) :: a2f_factor,estep,gaussfactor,gaussprefactor,gaussval,lambda_iso
 real(dp) :: omega,omegalog,omegastep,smear,tc_macmill,weight,xx
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: a2f_1d(:),a2f_1mom(:),a2f_1mom_int(:),a2flogmom(:)
 real(dp),allocatable :: a2flogmom_int(:),eli_smear(:,:,:),tmpa2f(:)

! *********************************************************************

!grid for the representation of alpha^2F (same as mka2f)
!WARNING : supposing that the maximum and minimum value of frequency
!have been defined in mkelph_linwid.

 omegastep = (elph_ds%omega_max-elph_ds%omega_min)/(elph_ds%na2f-one)

 nunit = get_unit()
 open (unit=nunit,file=fname,form='formatted',status='unknown',iostat=iost)
 if (iost /= 0) then
   MSG_ERROR("Opening file: " //trim(fname))
 end if

 write (msg,'(3a)')&
& '# Eliashberg function evaluated using only the irred q-points ',ch10,'#'
 call wrtout(nunit,msg,'COLL')

 write (msg,'(a,i5,2a,es16.8,2a,es16.8,2a,es16.8,2a)')&
& '# number of frequencies = ',elph_ds%na2f,ch10,         &
& '# omega_min = ',elph_ds%omega_min,ch10,                &
& '# omega_max = ',elph_ds%omega_max,ch10,                &
& '# step = ',omegastep,ch10,'#'
 call wrtout(nunit,msg,'COLL')


 nsmear=5
 estep=0.00002_dp !0.54422767 meV

 write (msg,'(a,i5,3a,f10.6,3a,f10.6,3a)')                &
& '# Using ',nsmear,' values for the gaussian smearing ',ch10,&
& '# starint from ',elph_ds%a2fsmear,' (Ha)',ch10,            &
& '# energy step of ',estep,' (Ha)',ch10,'#'
 call wrtout(nunit,msg,'COLL')

!e-ph quantities will be calculated for nsmear gaussian smearing values
!starting from elph_ds%a2fsmearwith an energy step of estep Hartree

 write (msg,'(3a)')'#      Smear(Ha) Lambda_Iso  isppol  <ln w> (K)    Tc_McMill (K) ',ch10,'#'
 call wrtout(nunit,msg,'COLL')

 ABI_ALLOCATE(a2f_1mom,(elph_ds%na2f))
 ABI_ALLOCATE(a2f_1mom_int,(elph_ds%na2f))
 ABI_ALLOCATE(a2flogmom,(elph_ds%na2f))
 ABI_ALLOCATE(a2flogmom_int,(elph_ds%na2f))
 ABI_ALLOCATE(a2f_1d,(elph_ds%na2f))
 ABI_ALLOCATE(tmpa2f,(elph_ds%na2f))
 ABI_ALLOCATE(eli_smear,(nsmear,elph_ds%nsppol,elph_ds%na2f))
 eli_smear(:,:,:)=zero

 do ismear=0,nsmear-1

   smear = elph_ds%a2fsmear+ismear*estep
   gaussprefactor = sqrt(piinv) / smear
   gaussfactor = one / smear

   do isppol=1,elph_ds%nsppol  ! spin pol channels

     a2f_1d(:) = zero
     tmpa2f(:) = zero

     do qptirred=1,elph_ds%nqptirred ! sum over irred qpoints
       do ibranch=1,elph_ds%nbranch

         if (abs(elph_ds%qgrid_data(qptirred,ibranch,isppol,1)) < tol10) cycle
         omega = elph_ds%omega_min
!        MG the weights in elph_ds%wtq(qptirred) are relative to the full grid qpt_full,
!        we need the mapping qirredtofull
         weight=elph_ds%wtq(elph_ds%qirredtofull(qptirred))
         a2f_factor=weight*elph_ds%qgrid_data(qptirred,ibranch,isppol,2)/abs(elph_ds%qgrid_data(qptirred,ibranch,isppol,1))

         do iomega=1,elph_ds%na2f
           xx = (omega-elph_ds%qgrid_data(qptirred,ibranch,isppol,1))*gaussfactor
           gaussval = gaussprefactor*exp(-xx*xx)
           tmpa2f(iomega) = tmpa2f(iomega) + gaussval*a2f_factor
           omega = omega+omegastep
         end do

       end do !end ibranch do
     end do !end qptirred

     a2f_1d(:)= tmpa2f(:)/(2*pi*elph_ds%n0(isppol))
     eli_smear(ismear+1,isppol,:)=a2f_1d(:) !save values

!    Do isotropic calculation of lambda and output lambda, Tc(MacMillan)
     a2f_1mom(:) = zero
     omega = elph_ds%omega_min

     do iomega=1,elph_ds%na2f
       if (abs(omega) > tol10) a2f_1mom(iomega) = two*a2f_1d(iomega)/abs(omega)
       omega=omega+omegastep
     end do

     call simpson_int(elph_ds%na2f,omegastep,a2f_1mom,a2f_1mom_int)
     lambda_iso = a2f_1mom_int(elph_ds%na2f)

!    Get log moment of alpha^2F
     a2flogmom(:) = zero
     omega = elph_ds%omega_min
     do iomega=1,elph_ds%na2f
       if (abs(omega) > tol10) then
         a2flogmom(iomega) = (two/lambda_iso)*a2f_1d(iomega)*log(abs(omega))/abs(omega)
       end if
       omega=omega+omegastep
     end do

     call simpson_int(elph_ds%na2f,omegastep,a2flogmom,a2flogmom_int)
     omegalog = exp(a2flogmom_int(elph_ds%na2f))

     tc_macmill = (omegalog/1.2_dp) * &
&     exp((-1.04_dp*(one+lambda_iso)) / (lambda_iso-elph_ds%mustar*(one+0.62_dp*lambda_iso)))

!    write data
     write(msg,'(a,5x,f10.6,f10.6,i5,2x,f12.7,2x,f12.6,2x,es16.8)')&
&     '# ',smear,lambda_iso,isppol,omegalog/kb_HaK,tc_macmill/kb_HaK
     call wrtout(nunit,msg,'COLL')

   end do !end isppol

 end do !ismear

 ABI_DEALLOCATE(a2f_1mom)
 ABI_DEALLOCATE(a2f_1mom_int)
 ABI_DEALLOCATE(a2flogmom)
 ABI_DEALLOCATE(a2flogmom_int)

!write to file
 write(msg,'(4a)')'#',ch10,'# Eliashberg function calculated for different gaussian smearing values',ch10
 call wrtout(nunit,msg,'COLL')

 do isppol=1,elph_ds%nsppol
   omega = elph_ds%omega_min
   write(nunit,'(a,i5)') '# smeared alpha2F for isppol = ',isppol
   do iomega=1,elph_ds%na2f
     write(nunit,'(6(f17.12,1x))')omega,eli_smear(:,isppol,iomega)
     omega=omega+omegastep
   end do
   write(nunit,*)
 end do

 ABI_DEALLOCATE(eli_smear)
 ABI_DEALLOCATE(a2f_1d)
 ABI_DEALLOCATE(tmpa2f)

 close (nunit)

end subroutine mka2fQgrid
!!***

end module m_elphon
!!***
