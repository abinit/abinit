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
!! Copyright (C) 2004-2019 ABINIT group (MVer, BXu, MG, JPC)
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
 use defs_elphon
 use m_abicore
 use m_krank
 use m_errors
 use m_xmpi
 use m_hdr
 use m_ebands

 use defs_datatypes,    only : ebands_t
 use m_fstrings,        only : int2char4
 use m_io_tools,        only : open_file, is_open, get_unit
 use m_time,            only : timein
 use m_numeric_tools,   only : wrap2_pmhalf, simpson, simpson_int
 use m_pptools,         only : printvtk
 use m_dynmat,          only : ftgam_init, ftgam
 use m_geometry,        only : phdispl_cart2red
 use m_kpts,            only : getkgrid, smpbz
 use m_crystal,         only : crystal_t
 use m_ifc,             only : ifc_type
 use m_nesting,         only : mknesting, bfactor
 use m_anaddb_dataset,  only : anaddb_dataset_type
 use m_eliashberg_1d,   only : eliashberg_1d
 use m_iogkk,           only : read_el_veloc,  read_gkk
 use m_bz_mesh,         only : make_path
 use m_fstab,           only : mkqptequiv
 use m_epweights,       only : d2c_weights, ep_el_weights, ep_fs_weights
 use m_a2ftr,           only : mka2f_tr, mka2f_tr_lova, get_tau_k
 use m_symkpt,          only : symkpt

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
!!      get_nv_fs_temp,get_rank,get_tau_k,get_veloc_tr,hdr_bcast
!!      hdr_fort_read,hdr_free,integrate_gamma,integrate_gamma_tr
!!      integrate_gamma_tr_lova,mka2f,mka2f_tr,mka2f_tr_lova,mka2fqgrid
!!      mkfskgrid,mknesting,mkph_linwid,mkqptequiv,order_fs_kpts,outelph
!!      printvtk,rchkgsheader,read_el_veloc,symkpt,timein,wrap2_pmhalf,wrtout
!!      xmpi_bcast
!!
!! SOURCE

subroutine elphon(anaddb_dtset,Cryst,Ifc,filnam,comm)

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
 integer,allocatable :: pair2red(:,:), red2pair(:,:), bz2ibz_smap(:,:)
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

 if (anaddb_dtset%ep_prt_yambo==1) then
   elph_ds%nFSband = nband
   elph_ds%minFSband = 1
   elph_ds%maxFSband = nband
 end if

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
   ABI_ALLOCATE(bz2ibz_smap, (6, elph_ds%k_phon%nkpt))

   wtk_fullbz(:) = one/dble(elph_ds%k_phon%nkpt) !weights normalized to unity
   call symkpt(0,cryst%gmet,indkpt1,0,elph_ds%k_phon%kpt,elph_ds%k_phon%nkpt,elph_ds%k_phon%new_nkptirr,&
&   Cryst%nsym,Cryst%symrec,timrev,wtk_fullbz,wtk_folded, bz2ibz_smap, xmpi_comm_self)

   ABI_FREE(bz2ibz_smap)

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

   elph_ds%k_fine%krank = elph_ds%k_phon%krank%copy()

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
     symrankkpt = elph_ds%k_phon%krank%get_rank (kpt)
     iFSkpq = elph_ds%k_phon%krank%invrank(symrankkpt)
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
     symrankkpt = elph_ds%k_phon%krank%get_rank (kpt)
     iFSkpq = elph_ds%k_phon%krank%invrank(symrankkpt)
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
   pair2red = 0

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
   red2pair = 0
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
 type(krank_t) :: krank
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

 krank = krank_new(elph_ds%k_phon%nkpt, elph_ds%k_phon%kpt)

 ABI_ALLOCATE(nestfactor,(nqptirred))

!NOTE: weights are not normalised, the normalisation factor in reintroduced in bfactor
 call bfactor(elph_ds%k_phon%nkpt,elph_ds%k_phon%kpt,nqptirred,qirred,krank,&
& elph_ds%k_phon%nkpt,elph_ds%k_phon%wtk,elph_ds%nFSband,nestfactor)

 ABI_DEALLOCATE(qirred)
 call krank%free()


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
!!      destroy_kptrank,get_rank,mkkptrank,sort_int,wrap2_pmhalf,wrtout
!!
!! SOURCE

subroutine mkFSkgrid (elph_k, nsym, symrec, timrev)

 use m_sort

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
 elph_k%krank = krank_new(elph_k%nkptirr, elph_k%kptirr)
 ABI_ALLOCATE(rankallk,(elph_k%krank%max_rank))

!elph_k%krank%invrank is used as a placeholder in the following loop
 rankallk = -1
 elph_k%krank%invrank = -1

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

       symrankkpt = elph_k%krank%get_rank (kpt)

!      is the kpt on the full grid (may have lower symmetry than full spgroup)
!      is kpt among the full FS kpts found already?
       if (elph_k%krank%invrank(symrankkpt) == -1) then
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

         elph_k%krank%invrank(symrankkpt) = elph_k%nkpt
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
 call elph_k%krank%free()

!make proper full rank arrays
 elph_k%krank = krank_new(elph_k%nkpt, elph_k%kpt)

!find correspondence table between irred FS kpoints and a full one
 ABI_ALLOCATE(elph_k%irr2full,(elph_k%nkptirr))
 elph_k%irr2full(:) = 0

 do ikpt1=1,elph_k%nkptirr
   symrankkpt = elph_k%krank%get_rank (elph_k%kptirr(:,ikpt1))
   elph_k%irr2full(ikpt1) = elph_k%krank%invrank(symrankkpt)
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
       symrankkpt = elph_k%krank%get_rank (kpt)
       ikpt2 = elph_k%krank%invrank(symrankkpt)
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
 use m_epweights,      only : d2c_wtq, ep_ph_weights

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
   call ifc%fourq(cryst,elph_ds%k_fine%kpt(:,iFSqpt),phfrq(:,iFSqpt),displ_cart,out_eigvec=pheigvec(:,iFSqpt))
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

     call ifc%fourq(cryst,kpt(:,iFSqpt),phfrq(:,iFSqpt),displ_cart,out_eigvec=pheigvec)

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
!  electron lifetimes eq 4.48 in [[cite:Grimvall1981]] electron phonon coupling in Metals (with T dependency). Also 5.69-5.72, 5.125, section 3.4
!  phonon lifetimes eq 19 in Savrasov PhysRevB.54.16487 [[cite:Savrasov1996]] (T=0)
!  a first T dependent expression in Allen PRB 6 2577 [[cite:Allen1972]] eq 10. Not sure about the units though
!
       do itemp = 1, ntemp
         temp = (itemp-1)*10._dp*kb_HaK
         linewidth_integrand(iomega, itemp) = a2f_1d(iomega) * (fermi_dirac(omega,zero,temp) + bose_einstein(omega,temp))
       end do
     end if
     omega=omega + domega
   end do
!
!  From Allen PRL 59 1460 [[cite:Allen1987]]
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

!!****f* m_elphon/order_fs_kpts
!!
!! NAME
!! order_fs_kpts
!!
!! FUNCTION
!! This routine re-orders the kpoints on the standard grid which belong
!!  to the Fermi surface: put them in increasing z, then y,  then x
!!
!! INPUTS
!!   nkptirr = number of irreducible FS kpoints
!!   nkpt = input nkpt from header
!!   kptns = input kpt from header
!!
!! OUTPUT
!!   FSirredtoGS = mapping of irreducible kpoints to GS set
!!   kptirr = irreducible FS kpoint coordinates
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      destroy_kptrank,mkkptrank,wrap2_pmhalf
!!
!! SOURCE

subroutine order_fs_kpts(kptns, nkpt, kptirr,nkptirr,FSirredtoGS)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkptirr
 integer,intent(in) :: nkpt

!arrays
 integer,intent(out) :: FSirredtoGS(nkptirr)
 real(dp),intent(in) :: kptns(3,nkpt)
 real(dp),intent(out) :: kptirr(3,nkptirr)

!Local variables-------------------------------
!scalars
 integer :: irank,ikpt,jkpt,kkpt,new, ik
 real(dp) :: res
 type(krank_t) :: krank
!arrays
 integer :: kptirrank(nkptirr)

! *************************************************************************

!rank is used to order kpoints
 krank = krank_new(nkpt, kptns)

 ik=1
 do ikpt=1,nkpt
   irank = krank%get_rank(kptns(:,ikpt))
!  add kpt to FS kpts, in order, increasing z, then y, then x !
   new = 1
!  look for position to insert kpt ikpt among irredkpts already found
   do jkpt=1,ik-1
     if (kptirrank(jkpt) > irank) then
!      shift all the others up
       do kkpt=ik-1,jkpt,-1
         kptirr(:,kkpt+1) = kptirr(:,kkpt)
         kptirrank(kkpt+1) = kptirrank(kkpt)
         FSirredtoGS(kkpt+1) = FSirredtoGS(kkpt)
       end do
!      insert kpoint ikpt
       call wrap2_pmhalf(kptns(1,ikpt),kptirr(1,jkpt),res)
       call wrap2_pmhalf(kptns(2,ikpt),kptirr(2,jkpt),res)
       call wrap2_pmhalf(kptns(3,ikpt),kptirr(3,jkpt),res)

       kptirrank(jkpt) = irank
       FSirredtoGS(jkpt) = ikpt
       new=0
       exit
     end if
   end do
!  ikpt not counted yet and higher rank than all previous
   if (new == 1) then
     call wrap2_pmhalf(kptns(1,ikpt),kptirr(1,ikpt),res)
     call wrap2_pmhalf(kptns(2,ikpt),kptirr(2,ikpt),res)
     call wrap2_pmhalf(kptns(3,ikpt),kptirr(3,ikpt),res)
     kptirrank(ik) = irank
     FSirredtoGS(ik) = ikpt
   end if
   ik=ik+1
 end do

 call krank%free()

end subroutine order_fs_kpts
!!***

!!****f* m_elphon/ep_setupqpt
!!
!! NAME
!! ep_setupqpt
!!
!! FUNCTION
!!  set up qpoint grid for elphon.
!!  2 modes, either uniform grid from anaddb input nqpt
!!  or take qpt from anaddb input (explicitly listed)
!!
!! INPUTS
!!   crystal>crystal_t>=data type gathering info on the crystalline structure.
!!   anaddb_dtset=dataset with input variables
!!     %qgrid_type gives type of q grid 1=uniform 2=take from input
!!     %ep_nqpt    number of auxiliary qpoints
!!     %ep_qptlist list of qpoints,
!!
!! OUTPUT
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      getkgrid,smpbz,symkpt,wrap2_pmhalf,wrtout
!!
!! NOTES
!!
!! SOURCE

subroutine ep_setupqpt (elph_ds,crystal,anaddb_dtset,qptrlatt,timrev)

!Arguments -------------------------------
!scalars
 integer, intent(in) :: timrev
 type(crystal_t),intent(in) :: crystal
 type(anaddb_dataset_type), intent(in) :: anaddb_dtset
 type(elph_type), intent(inout) :: elph_ds
!arrays
 integer, intent(out) :: qptrlatt(3,3)

!Local variables -------------------------
!scalars
 integer :: nqshft,option,iqpt, nqpt1
 integer :: iscf,mqpt,iout,berryopt,nqpt_computed
 real(dp) :: qptrlen, res
 character(len=500) :: message
!arrays
 integer :: vacuum(3)
 integer,allocatable :: indqpt1(:)
 real(dp) :: kpt(3)
 integer, allocatable :: bz2ibz_smap(:,:)
 real(dp),allocatable :: wtq_folded(:)
 real(dp), allocatable :: wtq(:),qpt_full(:,:),tmpshifts(:,:)

! *********************************************************************

!default is to expect a uniform grid
 elph_ds%tuniformgrid = 1

!if we use the normal grid way of generating the qpoints:
 if (anaddb_dtset%qgrid_type==1) then
!  qpoint lattice vectors (inverse, like kptrlatt)
   qptrlatt(:,:)=0
   qptrlatt(1,1)=anaddb_dtset%ngqpt(1)
   qptrlatt(2,2)=anaddb_dtset%ngqpt(2)
   qptrlatt(3,3)=anaddb_dtset%ngqpt(3)

   if (anaddb_dtset%nqshft /= 1) then
!    try to reduce the qpoint grid to a single qshift, otherwise stop
!    dummy args for call to getkgrid
     vacuum(:) = 0
     iscf = 3

     mqpt = anaddb_dtset%ngqpt(1)*anaddb_dtset%ngqpt(2)*anaddb_dtset%ngqpt(3)*anaddb_dtset%nqshft
     ABI_ALLOCATE(qpt_full,(3,mqpt))
     ABI_ALLOCATE(wtq,(mqpt))
     ABI_ALLOCATE(tmpshifts,(3,MAX_NSHIFTK))

     wtq(:) = one

     tmpshifts(:,:) = zero
     tmpshifts(:,1:4) = anaddb_dtset%q1shft(:,:)

     iout=6

     berryopt = 1

!    just call with identity, to get full set of kpts in qpt_full, but
!    reduce qshfts

     nqshft=anaddb_dtset%nqshft
     call getkgrid(0,0,iscf,qpt_full,3,qptrlatt,qptrlen, &
&     1,mqpt,nqpt_computed,nqshft,1,crystal%rprimd,tmpshifts,crystal%symafm, &
&     crystal%symrel,vacuum,wtq)
     ABI_DEALLOCATE(qpt_full)
     ABI_DEALLOCATE(wtq)
     ABI_DEALLOCATE(tmpshifts)

     if (anaddb_dtset%nqshft /= 1) then
       write (message,'(a,i0)')&
&       ' multiple qpt shifts not treated yet (should be possible), nqshft= ', anaddb_dtset%nqshft
       MSG_ERROR(message)
     end if
   end if  ! end multiple shifted qgrid


   write(message,'(a,9(i0,1x))')' elphon : enter smpbz with  qptrlatt = ',qptrlatt
   call wrtout(std_out,message,'COLL')

   option=1
!  mqpt=anaddb_dtset%ngqpt(1)*anaddb_dtset%ngqpt(2)*anaddb_dtset%ngqpt(3)*anaddb_dtset%nqshft
   mqpt= qptrlatt(1,1)*qptrlatt(2,2)*qptrlatt(3,3) &
&   +qptrlatt(1,2)*qptrlatt(2,3)*qptrlatt(3,1) &
&   +qptrlatt(1,3)*qptrlatt(2,1)*qptrlatt(3,2) &
&   -qptrlatt(1,2)*qptrlatt(2,1)*qptrlatt(3,3) &
&   -qptrlatt(1,3)*qptrlatt(2,2)*qptrlatt(3,1) &
&   -qptrlatt(1,1)*qptrlatt(2,3)*qptrlatt(3,2)

   ABI_ALLOCATE(qpt_full,(3,mqpt))
   iout = 6
   call smpbz(anaddb_dtset%brav,iout,qptrlatt,mqpt,elph_ds%nqpt_full,anaddb_dtset%nqshft,option,anaddb_dtset%q1shft,qpt_full)


!  save the q-grid for future reference
   ABI_ALLOCATE(elph_ds%qpt_full,(3,elph_ds%nqpt_full))

!  reduce qpt_full to correct zone
   do iqpt=1,elph_ds%nqpt_full
     call wrap2_pmhalf(qpt_full(1,iqpt),kpt(1),res)
     call wrap2_pmhalf(qpt_full(2,iqpt),kpt(2),res)
     call wrap2_pmhalf(qpt_full(3,iqpt),kpt(3),res)
     qpt_full(:,iqpt) = kpt
     elph_ds%qpt_full(:,iqpt)=kpt
   end do
   ABI_DEALLOCATE(qpt_full)

 else if (anaddb_dtset%qgrid_type==2) then ! use explicit list of qpoints from anaddb input
   qptrlatt(:,:)=0
   qptrlatt(1,1)=1
   qptrlatt(2,2)=1
   qptrlatt(3,3)=1

   elph_ds%nqpt_full=anaddb_dtset%ep_nqpt
   ABI_ALLOCATE(elph_ds%qpt_full,(3,elph_ds%nqpt_full))

   elph_ds%qpt_full = anaddb_dtset%ep_qptlist

   elph_ds%tuniformgrid = 0
 end if ! type of qgrid for elphon

!=================================================================
!Calculate weights, needed to estimate lambda using the weighted
!sum of the uninterpolated e-ph matrix elements
!=================================================================
 call wrtout(std_out,' setqgrid : calling symkpt to find irred q points',"COLL")

 ABI_ALLOCATE(indqpt1,(elph_ds%nqpt_full))
 ABI_ALLOCATE(wtq_folded,(elph_ds%nqpt_full))
 ABI_ALLOCATE(wtq,(elph_ds%nqpt_full))
 ABI_ALLOCATE(bz2ibz_smap, (6, elph_ds%nqpt_full))

 wtq(:) = one/dble(elph_ds%nqpt_full) !weights normalized to unity

!
!NOTE: this reduction of irred qpt may not be identical to that in GKK file
!which would be more practical to use.
!
 iout=0 !do not write to ab_out
!should we save indqpt1 for use inside elph_ds?
 call symkpt(0,crystal%gmet,indqpt1,iout,elph_ds%qpt_full,elph_ds%nqpt_full,nqpt1,crystal%nsym,crystal%symrec,&
& timrev,wtq,wtq_folded, bz2ibz_smap, xmpi_comm_self)

 ABI_FREE(bz2ibz_smap)

 write (message,'(2a,i0)')ch10,' Number of irreducible q-points = ',nqpt1
 call wrtout(std_out,message,'COLL')
 elph_ds%nqptirred=nqpt1

 call wrtout(std_out,' === Irreducible q points with weights ==== ','COLL')

 do iqpt=1,elph_ds%nqpt_full
   if (wtq_folded(iqpt) /= zero) then
     write (message,'(1x,i4,a2,4es16.8)')iqpt,') ',elph_ds%qpt_full(:,iqpt),wtq_folded(iqpt)
     call wrtout(std_out,message,'COLL')
   end if
 end do

 call wrtout(std_out,ch10,'COLL')

 ABI_ALLOCATE(elph_ds%wtq,(elph_ds%nqpt_full))

 elph_ds%wtq(:)=wtq_folded(:)
!MEMO indqpt could be useful to test the qgrid read by abinit
 ABI_DEALLOCATE(indqpt1)
 ABI_DEALLOCATE(wtq_folded)
 ABI_DEALLOCATE(wtq)

end subroutine ep_setupqpt
!!***

!!****f* ABINIT/mkph_linwid
!!
!! NAME
!! mkph_linwid
!!
!! FUNCTION
!!  Calculate the phonon linewidths on a trajectory in q space
!!
!! INPUTS
!!  Cryst<crystal_t>=Info on the unit cell and symmetries.
!!  Ifc<ifc_type>=Object containing the interatomic force constants.
!!  elph_ds = datastructure with phonon matrix elements
!!  nqpath = dimension of qpath_vertices
!!  qpath_vertices = vertices of reciprocal space trajectory
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      ftgam,ftgam_init,gam_mult_displ,ifc_fourq,make_path,phdispl_cart2red
!!      wrap2_pmhalf,wrtout,zgemm
!!
!! SOURCE

subroutine mkph_linwid(Cryst,ifc,elph_ds,nqpath,qpath_vertices)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqpath
 type(crystal_t),intent(in) :: Cryst
 type(ifc_type),intent(in) :: ifc
 type(elph_type),intent(inout) :: elph_ds
!arrays
 real(dp),intent(in) :: qpath_vertices(3,nqpath)

!Local variables-------------------------------
!scalars
 integer :: ibranch,natom,ii,indx,ipoint,nbranch,nqbz,nsppol,nrpt
 integer :: isppol,jbranch,qtor,unit_bs,unit_lambda,unit_lwd,npt_tot
 real(dp) :: diagerr,res
 character(len=500) :: msg
 character(len=fnlen) :: fname,base_name
!arrays
 integer :: ndiv(nqpath-1)
 integer, allocatable :: indxprtqpt(:)
 real(dp),parameter :: c0(2)=(/0._dp,0._dp/),c1(2)=(/1._dp,0._dp/)
 real(dp) :: displ_cart(2,3*Cryst%natom,3*Cryst%natom)
 real(dp) :: displ_red(2,3*Cryst%natom,3*Cryst%natom)
 real(dp) :: eigval(3*Cryst%natom)
 real(dp) :: gam_now(2,(3*Cryst%natom)**2)
 real(dp) :: imeigval(3*Cryst%natom)
 real(dp) :: lambda(3*Cryst%natom)
 real(dp) :: pheigvec(2*3*Cryst%natom*3*Cryst%natom),phfrq_tmp(3*Cryst%natom)
 real(dp) :: qpt(3),redkpt(3)
 real(dp) :: tmpgam1(2,3*Cryst%natom,3*Cryst%natom)
 real(dp) :: tmpgam2(2,3*Cryst%natom,3*Cryst%natom)
 real(dp), allocatable :: coskr(:,:), sinkr(:,:),finepath(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

 natom     = Cryst%natom
 nbranch   = elph_ds%nbranch
 nsppol    = elph_ds%nsppol
 base_name = elph_ds%elph_base_name
 nrpt = ifc%nrpt

!===================================================================
!Definition of the q path along which ph linwid will be interpolated
!===================================================================
 call make_path(nqpath,qpath_vertices,Cryst%gmet,'G',20,ndiv,npt_tot,finepath)
 ABI_ALLOCATE(indxprtqpt,(npt_tot))
 indxprtqpt = 0

!==========================================================
!Open _LWD file and write header
!==========================================================
 fname=trim(base_name) // '_LWD'
 if (open_file(fname,msg,newunit=unit_lwd,status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if

 write (unit_lwd,'(a)')       '#'
 write (unit_lwd,'(a)')       '# ABINIT package : Phonon linewidth file'
 write (unit_lwd,'(a)')       '#'
 write (unit_lwd,'(a,i10,a)') '#  Phonon linewidths calculated on ',npt_tot,' points along the qpath'
 write (unit_lwd,'(a)')       '#  Description of the Q-path :'
 write (unit_lwd, '(a,i10)')  '#  Number of line segments = ',nqpath-1
 write (unit_lwd,'(a)')       '#  Vertices of the Q-path and corresponding index = '

 indx=1
 indxprtqpt(1) = 1
 indxprtqpt(npt_tot) = 1

 do ii=1,nqpath
   write (unit_lwd,'(a,3(e16.6,1x),i8)')'#  ',qpath_vertices(:,ii),indx
   if (ii<nqpath) then
     indx=indx+ndiv(ii)
     indxprtqpt(indx) = 1
   end if
 end do

 write (unit_lwd,'(a)')'#'

!==========================================================
!Open _BST file and write header
!==========================================================
 fname=trim(base_name) // '_BST'
 if (open_file(fname,msg,newunit=unit_bs,status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if

 write (unit_bs, '(a)')      '#'
 write (unit_bs, '(a)')      '# ABINIT package : Phonon band structure file'
 write (unit_bs, '(a)')      '#'
 write (unit_bs, '(a,i10,a)')'# Phonon BS calculated on ', npt_tot,' points along the qpath'
 write (unit_bs, '(a,i10)')  '# Number of line segments = ', nqpath-1
 indx=1
 do ii=1,nqpath
   write (unit_bs,'(a,3(E16.6,1x),i8)')'#  ',qpath_vertices(:,ii),indx
   if (ii<nqpath) indx=indx+ndiv(ii)
 end do
 write (unit_bs,'(a)')'#'

!MG20060606
!==========================================================
!open _LAMBDA file and write header
!contains \omega(q,n) and \lambda(q,n) and can be plotted using xmgrace
!==========================================================
 fname=trim(base_name) // '_LAMBDA'
 if (open_file(fname,msg,newunit=unit_lambda,status="unknown") /= 0) then
   MSG_ERROR(msg)
 end if

 write (unit_lambda,'(a)')      '#'
 write (unit_lambda,'(a)')      '# ABINIT package : Lambda file'
 write (unit_lambda,'(a)')      '#'
 write (unit_lambda,'(a,i10,a)')'#  Lambda(q,nu) calculated on ',npt_tot,' Q-points'
 write (unit_lambda,'(a)')      '# Description of the Q-path :'
 write (unit_lambda,'(a,i10)')  '# Number of line segments = ',nqpath-1
 write (unit_lambda,'(a)')      '# Vertices of the Q-path and corresponding index = '

 indx=1
 do ii=1,nqpath
   write (unit_lambda,'(a,3(E16.6,1x),i8)')'#  ',qpath_vertices(:,ii),indx
   if (ii<nqpath) indx=indx+ndiv(ii)
 end do
 write (unit_lambda,'(a)')'#'
 write (unit_lambda,'(a)')'# index frequency lambda(q,n) frequency lambda(q,n) .... lambda_tot'
 write (unit_lambda,'(a)')'#'

!real space to q space
 qtor=0

!initialize the maximum phonon frequency
 elph_ds%omega_min = zero
 elph_ds%omega_max = zero

 ABI_ALLOCATE(coskr, (npt_tot,nrpt))
 ABI_ALLOCATE(sinkr, (npt_tot,nrpt))
 call ftgam_init(ifc%gprim, npt_tot, nrpt, finepath, ifc%rpt, coskr, sinkr)

 write (std_out,*) ' mkph_linwid : shape(elph_ds%gamma_qpt) = ',shape(elph_ds%gamma_qpt)
 nqbz =  SIZE(elph_ds%gamma_qpt,DIM=4)
 write(std_out,*) " nqbz =  SIZE(elph_ds%gamma_qpt,DIM=4) = ",nqbz
!
!Big do loop over spin polarizations
!could put in locally, so phonon stuff is not done twice...
!
 do isppol=1,nsppol
   indx=1

!  Output to the main output file
   write(msg,'(a,a)')ch10,&
&   ' Output of the linewidths for the first point of each segment. Linewidths are given in Hartree.'
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')

   write (std_out,*) ' mkph_linwid : elph_ds%ep_scalprod = ', elph_ds%ep_scalprod

   qtor = 0

!  Interpolation along specified path in q space
   do ipoint=1,npt_tot

!    Get qpoint along the path from qpath_vertices
     qpt(:) = finepath(:,ipoint)

     call wrap2_pmhalf(qpt(1),redkpt(1),res)
     call wrap2_pmhalf(qpt(2),redkpt(2),res)
     call wrap2_pmhalf(qpt(3),redkpt(3),res)
     qpt(:) = redkpt(:)
!
!    This reduced version of ftgkk supposes the kpoints have been integrated
!    in integrate_gamma. Do FT from real-space gamma grid to 1 qpt.
     call ftgam(ifc%wghatm,gam_now,elph_ds%gamma_rpt(:,:,isppol,:),natom,1,ifc%nrpt,qtor, &
&     coskr(ipoint,:), sinkr(ipoint,:))
!
!    get phonon freqs and eigenvectors anyway
!
     call ifc%fourq(cryst,qpt,phfrq_tmp,displ_cart,out_eigvec=pheigvec)
!
!    additional frequency factor for some cases
!
!    If the matrices do not contain the scalar product with the displ_cart vectors yet do it now
     if (elph_ds%ep_scalprod == 0) then

       call phdispl_cart2red(natom,Cryst%gprimd,displ_cart,displ_red)

       tmpgam2 = reshape (gam_now, (/2,nbranch,nbranch/))
       call gam_mult_displ(nbranch, displ_red, tmpgam2, tmpgam1)

       do jbranch=1,nbranch
         eigval(jbranch) = tmpgam1(1, jbranch, jbranch)
         imeigval(jbranch) = tmpgam1(2, jbranch, jbranch)

         if (abs(imeigval(jbranch)) > tol8) then
           write (msg,'(a,i0,a,es16.8)')' imaginary values for branch = ',jbranch,' imeigval = ',imeigval(jbranch)
           MSG_WARNING(msg)
         end if
       end do

     else if (elph_ds%ep_scalprod == 1) then
!
!      Diagonalize gamma matrix at qpoint (complex matrix).
!      MJV NOTE: gam_now is recast implicitly here to matrix
       call ZGEMM ( 'N', 'N', 3*natom, 3*natom, 3*natom, c1, gam_now, 3*natom,&
&       pheigvec, 3*natom, c0, tmpgam1, 3*natom)

       call ZGEMM ( 'C', 'N', 3*natom, 3*natom, 3*natom, c1, pheigvec, 3*natom,&
&       tmpgam1, 3*natom, c0, tmpgam2, 3*natom)

       diagerr = zero
       do ibranch=1,nbranch

         eigval(ibranch) = tmpgam2(1,ibranch,ibranch)

         do jbranch=1,ibranch-1
           diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))+abs(tmpgam2(2,jbranch,ibranch))
         end do
         do jbranch=ibranch+1,nbranch
           diagerr = diagerr + abs(tmpgam2(1,jbranch,ibranch))+abs(tmpgam2(2,jbranch,ibranch))
         end do
         diagerr = diagerr + abs(tmpgam2(2,ibranch,ibranch))
       end do

       if (diagerr > tol12) then
         write (msg,'(a,es14.6)')' Numerical error in diagonalization of gamma with phon eigenvectors: ', diagerr
         MSG_WARNING(msg)
       end if

     else
       write (msg,'(a,i0)')' Wrong value for elph_ds%ep_scalprod = ',elph_ds%ep_scalprod
       MSG_BUG(msg)
     end if ! end elph_ds%ep_scalprod if
!
!    ==========================================================
!    write data to files for each q point
!    ==========================================================
     write (unit_lwd,'(i5)', advance='no') indx
     do ii=1, nbranch
       write (unit_lwd,'(E16.5)',advance='no') eigval(ii)
     end do
     write (unit_lwd,*)

!    only print phonon BS for isppol 1: independent of electron spins
     if (isppol==1) then
       write (unit_bs,'(i5)', advance='no') indx
       do ii=1, nbranch
         write (unit_bs,'(E16.5)',advance='no') phfrq_tmp(ii)
       end do
       write (unit_bs,*)
     end if

     write (unit_lambda,'(i5)', advance='no') indx
     do ii=1,nbranch
       lambda(ii)=zero
       if (abs(phfrq_tmp(ii)) > tol10) lambda(ii)=eigval(ii)/(pi*elph_ds%n0(isppol)*phfrq_tmp(ii)**2)
       write (unit_lambda,'(es16.8)',advance='no')phfrq_tmp(ii),lambda(ii)
     end do
     write (unit_lambda,'(es16.8)',advance='no') sum(lambda)
     write (unit_lambda,*)

!    MG NOTE: I wrote a piece of code to output all these quantities using units
!    chosen by the user, maybe in version 5.2?
!    In this version the output of lambda(q,\nu) has been added

!    Output to the main output file, for first point in segment
     if(indxprtqpt(ipoint)==1)then
       write(msg,'(a,a,3es16.6,a,i4,a,a)')ch10,&
&       ' Q point =',qpt(:),'   isppol = ',isppol,ch10,&
&       ' Mode number    Frequency (Ha)  Linewidth (Ha)  Lambda(q,n)'
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       do ii=1,nbranch
         write(msg,'(i8,es20.6,2es16.6)' )ii,phfrq_tmp(ii),eigval(ii),lambda(ii)
         call wrtout(std_out,msg,'COLL')
         call wrtout(ab_out,msg,'COLL')
       end do
     end if

!    find max/min phonon frequency along path chosen
!    presumed to be representative of full BZ to within 10 percent
     elph_ds%omega_min = min(elph_ds%omega_min,1.1_dp*phfrq_tmp(1))
     elph_ds%omega_max = max(elph_ds%omega_max,1.1_dp*phfrq_tmp(nbranch))

     indx = indx+1
   end do !  end ipoint do

!  add blank lines to output files between sppol
   write(msg,'(a)' ) ''
   call wrtout(unit_lwd,msg,'COLL')
   call wrtout(unit_lambda,msg,'COLL')
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end do ! isppol

 ABI_DEALLOCATE(coskr)
 ABI_DEALLOCATE(sinkr)

 close(unit=unit_lwd)
 close(unit=unit_bs)
 close(unit=unit_lambda)

 ABI_DEALLOCATE(finepath)
 ABI_DEALLOCATE(indxprtqpt)

 write(std_out,*) ' elph_linwid : omega_min, omega_max = ',elph_ds%omega_min, elph_ds%omega_max

 DBG_EXIT("COLL")

end subroutine mkph_linwid
!!***

!!****f* ABINIT/get_fs_bands
!!
!! NAME
!! get_fs_bands
!!
!! FUNCTION
!! This routine determines the bands which contribute to the Fermi surface
!!
!! INPUTS
!!  eigenGS = ground state eigenvalues
!!  hdr = header from input GS file
!!  ep_b_min, ep_b_max=A non-zero value is used to impose certain bands.
!!  fermie=Fermi level.
!!  eigenGS(hdr%nband(1),hdr%nkpt,hdr%nsppol)=Energies.
!!
!! OUTPUT
!!  minFSband,maxFSband=Minimun and maximum index for the bands that cross the Fermi level
!!  nkptirr=Number of irreducible points for which there exist at least one band that crosses the Fermi level.
!!
!! TODO
!!  1) Indeces and dimensions should should be spin dependent.
!!  2) In the present status of the code, all the k-points in the IBZ are used!
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine get_fs_bands(eigenGS,hdr,fermie,ep_b_min,ep_b_max,minFSband,maxFSband,nkptirr)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ep_b_min, ep_b_max
 integer,intent(out) :: minFSband,maxFSband,nkptirr
 real(dp),intent(in) :: fermie
 type(hdr_type),intent(in) :: hdr
!arrays
 real(dp),intent(in) :: eigenGS(hdr%nband(1),hdr%nkpt,hdr%nsppol)

!Local variables-------------------------------
!scalars
 integer :: iband,ikpt,isppol,nband
 real(dp) :: epsFS,gausstol,gaussig
 character(len=500) :: message
 integer :: kpt_phonflag(hdr%nkpt)

! *************************************************************************

!supposes nband is equal for all kpts
 nband = hdr%nband(1)

!gausstol = minimum weight value for integration weights on FS
!should be set to reproduce DOS at Ef (Ref. PRB 34, 5065 [[cite:Lam1986]] p. 5067)
 gausstol = 1.0d-10

!use same band indices in both spin channels
 maxFSband=1
 minFSband=nband

!window of states around fermi Energy is contained in +/- epsFS
!should be adjusted to take into account a minimal but sufficient
!fraction of the kpoints: see the loop below.
!The 1000 is purely empirical!!!
!Should also take into account the density of kpoints.
!gaussig = width of gaussian energy window around fermi energy
!needed to get a good fraction of kpoints contributing to the FS

 gaussig = (maxval(eigenGS)-minval(eigenGS))/1000.0_dp

 write (message,'(a,f11.8,2a)')' get_fs_bands : initial energy window = ',gaussig,ch10,&
& ' The window energy will be increased until the full k-grid is inside the range'
 call wrtout(std_out,message,'COLL')

!NOTE: could loop back to here and change gaussig until we have
!a certain fraction of the kpoints in the FS region...
 nkptirr = 0

!Do not use restricted fermi surface: include all kpts -> one
 do while (nkptirr < hdr%nkpt)
   gaussig = gaussig*1.05_dp

!  we must take into account kpoints with states within epsFS:
   epsFS = gaussig*sqrt(log(one/(gaussig*sqrt(pi)*gausstol)))

!  check if there are eigenvalues close to the Fermi surface
!  (less than epsFS from it)
   kpt_phonflag(:) = 0

!  do for each sppol channel
   do isppol=1,hdr%nsppol
     do ikpt=1,hdr%nkpt
       do iband=1,nband
         if (abs(eigenGS(iband,ikpt,isppol) - fermie) < epsFS) then
           kpt_phonflag(ikpt) = 1
           if (iband > maxFSband) maxFSband = iband
           if (iband < minFSband) minFSband = iband
         end if
       end do
     end do
   end do ! isppol

!  if user imposed certain bands for e-p, make sure they are kept
   if (ep_b_min /= 0 .and. ep_b_min < minFSband) then
     minFSband = ep_b_min
   end if
   if (ep_b_max /= 0 .and. ep_b_max > maxFSband) then
     maxFSband = ep_b_max
   end if

!  number of irreducible kpoints (by all sym) contributing to the Fermi surface (to be completed by symops).
   nkptirr = sum(kpt_phonflag(:))
 end do

 write(std_out,*) ' Energy window around Fermi level= ',epsFS,' nkptirr= ',nkptirr

end subroutine get_fs_bands
!!***

!!****f* ABINIT/get_all_gkk2
!! NAME
!! get_all_gkk2
!!
!! FUNCTION
!! This routine determines where to store gkk2 matrix elements (disk or RAM)
!! and calls interpolate_gkk to calculate them.
!! This is the most time consuming step.
!!
!! INPUTS
!!   acell = lengths of unit cell vectors
!!   amu = masses of atoms
!!   atmfrc = atomic force constants
!!   dielt = dielectric tensor
!!   dipdip = dipole-dipole contribution flag
!!   dyewq0 =
!!   elph_ds = datastructure for elphon data and dimensions
!!   kptirr_phon = irreducible set of fermi-surface kpoints
!!   kpt_phon = full set of fermi-surface kpoints
!!   ftwghtgkk = weights for FT of matrix elements
!!   gmet = metric in reciprocal space
!!   indsym = indirect mapping of atoms under symops
!!   mpert = maximum number of perturbations
!!   msym = maximum number of symmetries (usually nsym)
!!   nsym = number of symmetries
!!   ntypat = number of types of atoms
!!   onegkksize = size of one gkk record, in bytes
!!   rmet = real-space metric
!!   rprim = unit cell lattice vectors (dimensionless)
!!   rprimd = real-space unit-cell lattice vectors
!!   rpt = points in real space for FT, in canonical coordinates
!!   symrel = symmetry operations in reduced real space
!!   trans = Atomic translations : xred = rcan + trans
!!   typat = array of types of atoms
!!   ucvol = unit cell volume
!!   xred = reduced coordinates of atoms
!!   zeff = Born effective charges
!!
!! OUTPUT
!!   elph_ds = calculated |gkk|^2 are in elph_ds%gkk2
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      interpolate_gkk
!!
!! SOURCE

subroutine get_all_gkk2(crystal,ifc,elph_ds,kptirr_phon,kpt_phon)

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: crystal
 type(ifc_type),intent(in) :: ifc
 type(elph_type),intent(inout) :: elph_ds
!arrays
 real(dp),intent(in) :: kpt_phon(3,elph_ds%k_phon%nkpt)
 real(dp),intent(in) :: kptirr_phon(3,elph_ds%k_phon%nkptirr)

!Local variables-------------------------------
!scalars
 integer :: iost,onediaggkksize,sz1,sz2,sz3,sz4
 real(dp) :: realdp_ex
 !character(len=500) :: msg

! *************************************************************************

 if (elph_ds%nsppol /= 1) then
   MSG_ERROR('get_all_gkk2: nsppol>1 not coded yet!')
 end if

 onediaggkksize = elph_ds%nbranch*elph_ds%k_phon%nkpt*kind(realdp_ex)

 elph_ds%unit_gkk2 = 37
 if (elph_ds%gkk2write == 0) then
   write(std_out,*) 'get_all_gkk2 : keep gkk2 in memory. Size = ',&
&   4.0*dble(elph_ds%k_phon%nkpt)*dble(onediaggkksize)/&
&   1024.0_dp/1024.0_dp, " Mb"
   sz1=elph_ds%nbranch
   sz2=elph_ds%ngkkband
   sz3=elph_ds%ngkkband
   sz4=elph_ds%k_phon%nkpt
   ABI_ALLOCATE(elph_ds%gkk2,(sz1,sz2,sz3,sz4,elph_ds%k_phon%nkpt,1))
   elph_ds%gkk2(:,:,:,:,:,:) = zero

 else if (elph_ds%gkk2write == 1) then
   write(std_out,*) 'get_all_gkk2 : About to open gkk2 file : '
   write(std_out,*) elph_ds%unit_gkk2,onediaggkksize
   open (unit=elph_ds%unit_gkk2,file='gkk2file',access='direct',&
&   recl=onediaggkksize,form='unformatted',status='new',iostat=iost)
   if (iost /= 0) then
     MSG_ERROR('error opening gkk2file as new')
   end if
!  rewind (elph_ds%unit_gkk2)
   write(std_out,*) 'get_all_gkk2 : disk file with gkk^2 created'
   write(std_out,*) '  calculate from real space gkk and phonon modes'
   write(std_out,*) '  gkk2write = 1 is forced: can take a lot of time! '
   write(std_out,*) ' size = ', 4.0*dble(onediaggkksize)*dble(elph_ds%k_phon%nkpt)/&
&   1024.0_dp/1024.0_dp, ' Mb'
 else
   MSG_ERROR('bad value of gkk2write')
 end if

!here do the actual calculation of |g_kk|^2
 MSG_ERROR("MGNOTE: interpolate_gkk is broken")
 ABI_UNUSED(kptirr_phon(1,1))
 call interpolate_gkk (crystal,ifc,elph_ds,kpt_phon)

 !MG: This was the old coding in version 7.6.2:

! call interpolate_gkk (elph_ds,kptirr_phon,kpt_phon,natom,nrpt,phon_ds,rcan,wghatm)
!
! and interpolate_gkk had the prototype:
!
!subroutine interpolate_gkk(elph_ds,kpt_phon,gprim,natom,nrpt,phon_ds,rpt,wghatm)

! hence we were associating kpt_phon to gprim!

end subroutine get_all_gkk2
!!***

!!****f* ABINIT/interpolate_gkk
!! NAME
!! interpolate_gkk
!!
!! FUNCTION
!! This routine interpolates the gkk matrices for all q vectors
!! between points on the full kpt_phon grid.
!!
!! INPUTS
!!   elph_ds = elphon datastructure with data and dimensions
!!   kpt_phon = coordinates of all kpoints close to the FS
!!
!! OUTPUT
!!   elph_ds = modified gkq
!!
!! NOTES
!!  inspired to some extent by epcouple.f from the DecAFT package by J. Kay Dewhurst
!!  most inputs taken from mkifc.f
!!  in anaddb set ifcflag 1 such that the IFC are calculated in atmfrc prior to calling elphon
!!
!! PARENTS
!!      get_all_gkk2
!!
!! CHILDREN
!!      ftgkk,ifc_fourq,wrap2_pmhalf,zhpev
!!
!! SOURCE

subroutine interpolate_gkk(crystal,ifc,elph_ds,kpt_phon)

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: crystal
 type(ifc_type),intent(in) :: ifc
 type(elph_type),intent(inout) :: elph_ds
!arrays
 real(dp),intent(in) :: kpt_phon(3,elph_ds%k_phon%nkpt)

!Local variables-------------------------------
  ! output variables for dfpt_phfrq
! variables for zhpev
! variables for phonon interpolation
!scalars
 integer :: i1,i2,ikpt_phon2,iFSqpt,ib1,ib2,ier,ii
 integer :: iost,isppol,qtor,natom
 integer :: sz1,sz2,sz3,sz4,unit_gkkp
 real(dp) :: qphnrm,res
 !character(len=500) :: msg
!arrays
 real(dp) :: gprim(3,3)
 real(dp) :: displ(2,elph_ds%nbranch,elph_ds%nbranch),eigval(3*crystal%natom)
 real(dp) :: eigvec(3*3*crystal%natom*3*crystal%natom)
 real(dp) :: pheigvec(2*elph_ds%nbranch*elph_ds%nbranch)
 real(dp) :: phfrq_tmp(elph_ds%nbranch),qphon(3),redkpt(3)
 real(dp),allocatable :: gkk2_diag_tmp(:,:,:,:),gkk2_tmp(:,:,:,:,:,:,:)
 real(dp),allocatable :: matrx(:,:),zhpev1(:,:)
 real(dp),allocatable :: zhpev2(:)

! *************************************************************************

!
!NOTE: mjv 18/5/2008 reverted to old style of ftgkk with all kpt done together.
!may want to modify this later to use the new cleaner format with 1 FT at a
!time.
!
 write(std_out,*) 'interpolate_gkk : enter'

 natom = crystal%natom
 gprim = ifc%gprim

 if (elph_ds%nsppol /= 1) then
   MSG_ERROR("interpolate_gkk not coded with nsppol>1 yet")
 end if
 isppol = 1


!------------------------------------------------------
!complete dynamical matrices for all qpts between points
!on full kpt grid (interpolation from IFC)
!------------------------------------------------------

 sz1=elph_ds%ngkkband;sz2=elph_ds%nbranch
 sz3=elph_ds%k_phon%nkpt;sz4=elph_ds%nFSband
!allocate (gkk_tmp(2,sz1,sz1,sz2,sz2,1,1))
!DEBUG
!allocate (gkk_tmp_full(2,sz1,sz1,sz2,elph_ds%nFSband,sz3))
!allocate (gkk_tmp_full(2,s2,sz4,sz4,sz3))
!ENDDEBUG
 ABI_ALLOCATE(gkk2_tmp,(2,sz1,sz1,sz2,sz2,sz3,1))
 ABI_ALLOCATE(gkk2_diag_tmp,(sz1,sz1,sz2,sz3))
 ABI_ALLOCATE(zhpev1,(2,2*3*natom-1))
 ABI_ALLOCATE(zhpev2,(3*3*natom-2))
 ABI_ALLOCATE(matrx,(2,(3*natom*(3*natom+1))/2))

 qphnrm = one
!in this part use the inverse Fourier transform to get 1 (arbitrary) qpt at a
!time
 ii = 0
 qtor = 0
 unit_gkkp = 150
 open (unit=unit_gkkp,file='gkkp_file_ascii',form='formatted',status='unknown',iostat=iost)
 if (iost /= 0) then
   MSG_ERROR("error opening gkkpfile as new")
 end if

!loop over all FS pairs.
!do ikpt1=1,elph_ds%k_phon%nkptirr
!do iFSqpt=1,elph_ds%k_phon%nkpt

!
!this should run through the sparse mesh of 2x2x2 kpoints
!
 do iFSqpt=1,elph_ds%k_phon%nkpt
   res = 2.0_dp*(kpt_phon(1,iFSqpt)+one)
   if (abs(res-int(res)) > tol10) cycle
   res = 2.0_dp*(kpt_phon(2,iFSqpt)+one)
   if (abs(res-int(res)) > tol10) cycle
   res = 2.0_dp*(kpt_phon(3,iFSqpt)+one)
   if (abs(res-int(res)) > tol10) cycle

!  do ikpt1=1,1
!
!  NOTE: should be very easy to parallelize!
!
!  write(std_out,*) ' interpolate_gkk : ikpt1 = ',ikpt1, ' / ', elph_ds%k_phon%nkptirr
   write(std_out,*) ' interpolate_gkk : ikpt1 = ',iFSqpt, ' / ', elph_ds%k_phon%nkpt

!  DEBUG
!  write(std_out,*) ' interpolate_gkk : Warning debug version'
!  cycle
!  ENDDEBUG

   gkk2_tmp(:,:,:,:,:,:,:) = zero

!  qphon = 1 - 2    ie.  1 = 2+qphon
   qphon(:) = kpt_phon(:,iFSqpt)

!  shouldnt be necessary here, but oh well
   call wrap2_pmhalf(qphon(1),redkpt(1),res)
   call wrap2_pmhalf(qphon(2),redkpt(2),res)
   call wrap2_pmhalf(qphon(3),redkpt(3),res)

   qphon(:) = redkpt(:)
   redkpt(1) = qphon(1)*gprim(1,1)+qphon(2)*gprim(1,2)+qphon(3)*gprim(1,3)
   redkpt(2) = qphon(1)*gprim(2,1)+qphon(2)*gprim(2,2)+qphon(3)*gprim(2,3)
   redkpt(3) = qphon(1)*gprim(3,1)+qphon(2)*gprim(3,2)+qphon(3)*gprim(3,3)
   write (unit_gkkp,*) 'qp= ', redkpt

   call ifc%fourq(crystal,qphon,phfrq_tmp,displ,out_eigvec=pheigvec)
   write (unit_gkkp,*) phfrq_tmp(:)*Ha_cmm1

   ii = ii+1
!  if(ii > 0 .and. ii < 1000) write(std_out,'(a,i5,3E16.6,2x)') &
!  &   ' wrote phfrq_tmp for time ', ii, phfrq_tmp
!  end if

!  phonon eigenvectors are in eigvec
!  real and imaginary parts
!  phonon displacements = eigvec/sqrt(M_i) are in displ
!  real and imaginary parts

!  DEBUG
!  test: uniform phonon frequency
!  phfrq_tmp(:) = 0.0001_dp
!  ENDDEBUG

!  FT gamma matrices for all kpt_phon points, and
!  for qpoint = qphon(:) = kpt_phon(ikpt_phon)

   call ftgkk(ifc%wghatm,gkk2_tmp,elph_ds%gkk_rpt,elph_ds%gkqwrite,&
&   elph_ds%gkk_rptwrite,gprim,1,&
&   natom,elph_ds%k_phon%nkpt,elph_ds%ngkkband,elph_ds%k_phon%nkpt,1,ifc%nrpt,elph_ds%nsppol,&
&   qtor,ifc%rpt,qphon,elph_ds%unit_gkk_rpt,elph_ds%unitgkq)

!  NOTE: Normally the eigenvectors of the gkk2_tmp should be the same as eigvec

!  Diagonalize gamma matrices at qpoint (complex matrix) for all kpt_phon.
!  Copied from dfpt_phfrq
   do ikpt_phon2=1,elph_ds%k_phon%nkpt
     res = 8.0_dp*(kpt_phon(1,ikpt_phon2)+one)
     if (abs(res-int(res)) > tol10) cycle
     res = 8.0_dp*(kpt_phon(2,ikpt_phon2)+one)
     if (abs(res-int(res)) > tol10) cycle
     res = 8.0_dp*(kpt_phon(3,ikpt_phon2)+one)
     if (abs(res-int(res)) > tol10) cycle

     write (unit_gkkp,*) 'kp= ', kpt_phon(:,ikpt_phon2)

     do ib1=1,elph_ds%ngkkband
       do ib2=1,elph_ds%ngkkband
         ier=0
         ii=1
         do i2=1,3*natom
           do i1=1,i2
             matrx(1,ii)=gkk2_tmp(1,ib1,ib2,i1,i2,ikpt_phon2,1)
             matrx(2,ii)=gkk2_tmp(2,ib1,ib2,i1,i2,ikpt_phon2,1)
             ii=ii+1
           end do
         end do
         call ZHPEV ('N','U',3*natom,matrx,eigval,eigvec,3*natom,zhpev1,&
&         zhpev2,ier)

         gkk2_diag_tmp(ib2,ib1,:,ikpt_phon2) = eigval(:)
         do i1=1,3*natom
           write (unit_gkkp,*) elph_ds%minFSband-1+ib1,elph_ds%minFSband-1+ib2,i1,&
&           eigval(i1)
         end do
       end do
     end do
   end do

   if (elph_ds%gkk2write == 1) then
     write(std_out,*) 'WARNING COMMENTED WRITE TO BINARY FILE!!!'
!    write (elph_ds%unit_gkk2,REC=iFSqpt) gkk2_diag_tmp(:,:,:,:)
     write(std_out,'(a,i4,4(2E16.6,2x))') ' gkk2 loop ', &
&     iFSqpt,gkk2_diag_tmp(1,1,:,1:2),gkk2_diag_tmp(1,1,:,elph_ds%k_phon%nkpt-1:elph_ds%k_phon%nkpt)
!    &    ikpt1,gkk2_tmp(:,1,1,1,1,1:2),gkk2_tmp(:,1,1,elph_ds%k_phon%nkpt-1:elph_ds%k_phon%nkpt)
   else if (elph_ds%gkk2write == 0) then
     elph_ds%gkk2(:,:,:,:,iFSqpt,isppol) = gkk2_diag_tmp(:,:,:,:)
!    elph_ds%gkk2(:,:,:,:,ikpt1) = gkk2_tmp
     write(std_out,*) ' interpolate_gkk : gkk2(b=1,b=1,:,kpt=1,iFSqpt) = '
     write(std_out,*) gkk2_diag_tmp(1,1,:,1)
   end if

 end do
!end do on iFSqpt

 ABI_DEALLOCATE(matrx)
 ABI_DEALLOCATE(zhpev1)
 ABI_DEALLOCATE(zhpev2)

end subroutine interpolate_gkk
!!***

!!****f* ABINIT/get_all_gkq
!!
!! NAME
!! get_all_gkq
!!
!! FUNCTION
!! This routine determines what to do with the initial qspace
!!   matrix elements of the electron phonon coupling (to disk or in memory),
!!   then reads those given in the gkk file and completes them
!!   (for kpts, then perturbations)
!!   01/2010: removed completion on qpoints here (MJV)
!!
!! INPUTS
!!   elph_ds = elphon datastructure with data and dimensions
!!   Cryst<crystal_t>=Info on the unit cell and on its symmetries.
!!   Ifc<ifc_type>=Object containing the interatomic force constants.
!!   Bst<ebands_t>=GS energies, occupancies and Fermi level.
!!   FSfullpqtofull = mapping of k+q to another k
!!   kphon_full2full = mapping of FS kpoints under symops
!!   kpt_phon = fermi surface kpoints
!!   %k_phon%wtk = integration weights for bands and kpoints near the FS
!!   gkk_flag = flag to
!!   nband = number of bands
!!   n1wf = number of file headers from perturbation calculations
!!      which are present in the initial gkk input file.
!!   onegkksize = size of one record of the new gkk output file, in bytes
!!   qpttoqpt = mapping of qpoints onto each other under symmetries
!!   unitgkk = fortran unit for initial gkk input file
!!   xred = reduced coordinates of atoms
!!
!! OUTPUT
!!   elph_ds%gkq = recip space elphon matrix elements.
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      complete_gkk,int2char4,read_gkk,wrtout
!!
!! SOURCE

subroutine get_all_gkq (elph_ds,Cryst,ifc,Bst,FSfullpqtofull,nband,n1wf,onegkksize,&
&    qpttoqpt,ep_prt_yambo,unitgkk,ifltransport)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1wf,nband,onegkksize,unitgkk,ep_prt_yambo,ifltransport
 type(crystal_t),intent(in) :: Cryst
 type(ifc_type),intent(in) :: ifc
 type(ebands_t),intent(in) :: Bst
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 integer,intent(in) :: qpttoqpt(2,Cryst%nsym,elph_ds%nqpt_full)

!Local variables-------------------------------
!scalars
 integer :: iost,ierr,me,sz2,sz3,sz4,sz5,sz6
 character(len=10) :: procnum
 character(len=500) :: message
 character(len=fnlen) :: fname
!arrays
 integer,allocatable :: gkk_flag(:,:,:,:,:)

! *************************************************************************

!attribute file unit number
 elph_ds%unitgkq = get_unit()

!============================================
!save gkk for all qpts in memory or to disk
!============================================

!DEBUG
!write(std_out,*) ' 4 bytes / ??'
!write(std_out,*) ' kind(real) = ', kind(one)
!write(std_out,*) ' elph_ds%ngkkband = ', elph_ds%ngkkband, '^2'
!write(std_out,*) ' elph_ds%nbranch = ', elph_ds%nbranch, '^2'
!write(std_out,*) ' elph_ds%k_phon%nkpt = ', elph_ds%k_phon%nkpt
!write(std_out,*) ' elph_ds%nsppol = ', elph_ds%nsppol
!write(std_out,*) ' elph_ds%nqptirred ', elph_ds%nqptirred
!ENDDEBUG

 write(message,'(a,f14.4,a)')&
& ' get_all_gkq : gkq file/array size = ',&
 4.0*dble(onegkksize)*dble(elph_ds%k_phon%my_nkpt)*dble(elph_ds%nqptirred)/1024.0_dp/1024.0_dp/1024.0_dp,' Gb'
 call wrtout(std_out,message,'COLL')

 if (elph_ds%gkqwrite == 0) then !calculate gkk(q) keeping all in memory

   call wrtout(std_out,' get_all_gkq : keep gkk(q) in memory ','COLL')

   sz2=elph_ds%ngkkband*elph_ds%ngkkband
   sz3=elph_ds%nbranch*elph_ds%nbranch
   sz4=elph_ds%k_phon%my_nkpt
   sz5=elph_ds%nsppol
   if (ifltransport == 3) then
     sz6=elph_ds%nqpt_full
   else
     sz6=elph_ds%nqptirred
   end if
   ABI_MALLOC_OR_DIE(elph_ds%gkk_qpt,(2,sz2,sz3,sz4,sz5,sz6), ierr)

   elph_ds%gkk_qpt = zero

 else if (elph_ds%gkqwrite == 1) then !calculate gkk(q) and write to file
   me = xmpi_comm_rank(xmpi_world)
   call int2char4(me,procnum)
   ABI_CHECK((procnum(1:1)/='#'),'Bug: string length too short!')
   fname=trim(elph_ds%elph_base_name) // "_P" // trim(procnum) // '_GKKQ'

   iost=open_file(file=fname,iomsg=message,newunit=elph_ds%unitgkq,access='direct',&
&   recl=onegkksize,form='unformatted')
   if (iost /= 0) then
     write (message,'(2a)')' get_all_gkq : ERROR- opening file ',trim(fname)
     MSG_ERROR(message)
   end if

   write (message,'(5a)')&
&   ' get_all_gkq : gkq matrix elements  will be written to file : ',trim(fname),ch10,&
&   ' Nothing is in files yet',ch10
   call wrtout(std_out,message,'COLL')

 else
   write(message,'(a,i0)')' gkqwrite must be 0 or 1 while it is : ',elph_ds%gkqwrite
   MSG_BUG(message)
 end if !if gkqwrite

!=====================================================
!read in g_kk matrix elements for all bands, kpoints,
!and calculated qpoints
!=====================================================
 call wrtout(std_out,' get_all_gkq : calling read_gkk to read in the g_kk matrix elements',"COLL")

 sz2=elph_ds%nbranch;sz3=elph_ds%k_phon%my_nkpt
 sz4=elph_ds%nsppol;sz5=elph_ds%nqpt_full
 ABI_MALLOC_OR_DIE(gkk_flag,(sz2,sz2,sz3,sz4,sz5), ierr)

 call read_gkk(elph_ds,Cryst,ifc,Bst,FSfullpqtofull,gkk_flag,n1wf,nband,ep_prt_yambo,unitgkk)

!if (elph_ds%symgkq ==1) then
!MJV 01/2010 removed the completion on qpt here: it should be done after FS integration
!so that everything is lighter in memory etc... (only irred qpt)
! if (0==1) then
 if (ifltransport == 3) then !  bxu, complete gkk is necessary

!  ==============================================================
!  complete gkk matrices for other qpoints on the full grid qpt_full
!  inspired and cannibalized from symdm9.f
!  FIXME: should add the possibility to copy over to other qpoints,
!  without full symmetrization, for testing purposes.
!  ==============================================================

   write(message,'(4a)')ch10,&
&   ' get_all_gkq : calling complete_gkk to complete ',ch10,&
&   ' gkk matrices for other qpoints on the full grid'
   call wrtout(std_out,message,'COLL')

   call complete_gkk(elph_ds,gkk_flag,Cryst%gprimd,Cryst%indsym,&
&   Cryst%natom,Cryst%nsym,qpttoqpt,Cryst%rprimd,Cryst%symrec,Cryst%symrel)

   call wrtout(std_out,' get_all_gkq : out of complete_gkk','COLL')

 end if !symgkq

!TODO Do we need gkk_flag in elphon?
 ABI_DEALLOCATE(gkk_flag)

end subroutine get_all_gkq
!!***

!!****f* ABINIT/get_all_gkr
!! NAME
!! get_all_gkr
!!
!! FUNCTION
!! This routine determines what to do with the rspace
!! matrix elements of the el phon coupling (to disk or in memory),
!! then reads those given in the gkq file and Fourier Transforms them
!!
!! INPUTS
!!   elph_ds = elphon datastructure with data and dimensions
!!   gprim = reciprocal space lattice vectors
!!   natom = number of atoms
!!   nrpt = number of real-space points used for FT
!!   onegkksize = size of one record of the new gkk output file, in bytes
!!   rpt = positions of real-space points for FT
!!   qpt_full = qpoint coordinates
!!   wghatm = weights for real-space rpt in FT
!!
!! OUTPUT
!!   elph_ds%gkr = real space elphon matrix elements.
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      ftgkk
!!
!! SOURCE

subroutine get_all_gkr (elph_ds,gprim,natom,nrpt,onegkksize,rpt,qpt_full,wghatm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrpt,onegkksize
 type(elph_type),intent(inout) :: elph_ds
!arrays
 real(dp),intent(in) :: gprim(3,3),rpt(3,nrpt),qpt_full(3,elph_ds%nqpt_full)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon0,iost,qtor,sz2,sz3,sz4,sz5

! *************************************************************************

!
!WARNING : disk file used for large arrays gkk_rpt and
!(eventually) gkk2
!
!allocate (gkk_rpt(2,elph_ds%nbranch,elph_ds%nFSband,elph_ds%nFSband,&
!&  elph_ds%k_phon%nkpt,nrpt))
 elph_ds%unit_gkk_rpt = 36
!see if the gkk_rpt should be written to a file (only available option now)
 if (elph_ds%gkk_rptwrite == 1) then
!  file is not present : we need to do the FT
   open (unit=elph_ds%unit_gkk_rpt,file='gkk_rpt_file',access='direct',&
&   recl=onegkksize,form='unformatted',&
&   status='new',iostat=iost)
   if (iost /= 0) then
     MSG_ERROR('get_all_gkr : error opening gkk_rpt_file as new')
   end if
   write(std_out,*) ' get_all_gkr : will write real space gkk to a disk file.'
   write(std_out,*) ' size = ', 4.0*dble(onegkksize)*dble(nrpt)/&
&   1024.0_dp/1024.0_dp, ' Mb'

!  else if (elph_ds%gkk_rptwrite  == 0) then
 else
   write(std_out,*) ' get_all_gkr : will keep real space gkk in memory.'
   write(std_out,*) ' size = ', 4.0*dble(onegkksize)*dble(nrpt)/&
&   1024.0_dp/1024.0_dp, ' Mb'
   sz2=elph_ds%ngkkband*elph_ds%ngkkband
   sz3=elph_ds%nbranch*elph_ds%nbranch
   sz4=elph_ds%k_phon%nkpt
   sz5=elph_ds%nsppol
   ABI_ALLOCATE(elph_ds%gkk_rpt,(2,sz2,sz3,sz4,sz5,nrpt))
!  write(std_out,*) ' get_all_gkr: invalid value for gkk_rptwrite'
!  stop
 end if
 write(std_out,*) '    about to FT the recip space gkk to real space '
 qtor = 1

!
!NOTE: should be very easy to parallelize!
!
 ikpt_phon0 = 1
 call ftgkk (wghatm,elph_ds%gkk_qpt,elph_ds%gkk_rpt,&
& elph_ds%gkqwrite,elph_ds%gkk_rptwrite,gprim,1,natom,&
& elph_ds%k_phon%nkpt,elph_ds%ngkkband,elph_ds%k_phon%nkpt,elph_ds%nqpt_full,&
& nrpt,elph_ds%nsppol,qtor,rpt,qpt_full,elph_ds%unit_gkk_rpt,elph_ds%unitgkq)

!call ftgkk (elph_ds,gprim,ikpt_phon0,natom,nrpt,qtor,rpt,qpt_full,wghatm)
 write(std_out,*) ' get_all_gkr : done with FT of gkk to real space'

!No longer need the gkk_qpt?
!if (elph_ds%gkqwrite == 0) deallocate (elph_ds%gkk_qpt)

!!DEBUG
!Test the FT of the gkk elements.
!call test_ftgkk(elph_ds,gprim,natom,nrpt,rpt,qpt_full,wghatm)
!!ENDDEBUG

!DEBUG
!do irpt=1,nrpt
!do ipert1=1,elph_ds%nbranch
!write(std_out,'(6(F16.5,1x))') elph_ds%gkk_rpt(:,ipert1,1,1,1,irpt)
!end do
!end do
!ENDDEBUG

end subroutine get_all_gkr
!!***

!!****f* ABINIT/complete_gkk
!!
!! NAME
!! complete_gkk
!!
!! FUNCTION
!! Use the set of special q points calculated by the Monkhorst &
!! Pack Technique.
!! Check if all the information for the q points are present in
!! the DDB to determine the elphon interaction matrices
!! Generate the gkk matrices of the set of q points which
!! samples homogeneously the entire Brillouin zone.
!!
!! INPUTS
!! elph_ds = datastructure for elphon information (mainly
!!      matrix elements and dimensions)
!!   elph_ds%k_phon%full2full = kpt_phon index mapping under symops
!! gkk_flag = flag for existence of matrix element
!! gprimd(3,3)=dimensionful primitive translations in reciprocal space
!! indsym = map of atoms by inverses of symrels
!! natom=number of atoms in unit cell
!! nsym=number of space group symmetries
!! qpttoqpt = qpoint index mapping under symops
!! rprimd(3,3)=dimensionful primitive translations in real space
!! symrec(3,3,nsym)=3x3 matrices of the group symmetries (recip space)
!! symrel(3,3,nsym)=3x3 matrices of the group symmetries (real space)
!! tnons(3,nsym)=nonsymmorphic translations associated to symrel
!!
!! OUTPUT
!! elph_ds%gkk_qpt = gkk matrices for all qpts on a full mesh
!!
!! PARENTS
!!      get_all_gkq
!!
!! CHILDREN
!!      xmpi_sum,zgemm
!!
!! SOURCE

subroutine complete_gkk(elph_ds,gkk_flag,gprimd,indsym,natom,nsym,qpttoqpt,rprimd,symrec,symrel)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: indsym(4,nsym,natom)
 integer,intent(in) :: qpttoqpt(2,nsym,elph_ds%nqpt_full),symrec(3,3,nsym)
 integer,intent(in) :: symrel(3,3,nsym)
 integer,intent(inout) :: gkk_flag(elph_ds%nbranch,elph_ds%nbranch,elph_ds%k_phon%my_nkpt,elph_ds%nsppol,elph_ds%nqpt_full)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,ib1,ibranch,ieqqpt,ii, ierr,comm
 integer :: iqpt,isppol,isym
 integer :: itim,jbranch,jj,kk,ll
 integer :: neqqpt,symikpt_phon
 integer :: iatom,ancestor_iatom
 integer :: ik_this_proc, me,sz1,sz2

 real(dp),parameter :: tol=2.d-8
!arrays
 integer :: symmetrized_qpt(elph_ds%nqpt_full)
 real(dp) :: ss(3,3)
 real(dp) :: tmp_mat(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_mat2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),allocatable :: gkk_qpt_new(:,:,:,:,:),gkk_qpt_tmp(:,:,:,:,:)

 real(dp) :: ss_allatoms(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: c_one(2), c_zero(2)


! *********************************************************************

 c_one = (/one,zero/)
 c_zero = (/zero,zero/)

!Generation of the gkk matrices relative to the q points
!of the set which samples the entire Brillouin zone

 comm = xmpi_world
 me = xmpi_comm_rank(comm)

 symmetrized_qpt(:) = -1

!FIXME bxu, why set it to 1?
!isppol=1

 sz1=elph_ds%ngkkband*elph_ds%ngkkband
 sz2=elph_ds%nbranch*elph_ds%nbranch

!these arrays are not parallelized, to enable symmetrization: syms swap k-points.
 ABI_ALLOCATE(gkk_qpt_new,(2,sz1,sz2,elph_ds%k_phon%nkpt,elph_ds%nsppol))
 ABI_ALLOCATE(gkk_qpt_tmp,(2,sz1,sz2,elph_ds%k_phon%nkpt,elph_ds%nsppol))

 do iqpt=1,elph_ds%nqpt_full

!  Already symmetrized?
   if (symmetrized_qpt(iqpt) == 1) cycle

   gkk_qpt_new(:,:,:,:,:) = zero
!   gkk_qpt_tmp(:,:,:,:,:) = zero

!  loop over qpoints equivalent to iqpt
   neqqpt=0
!  do not use time reversal symmetry to complete the qpoints:
!  do not know what happens to the gamma matrices
!  itim=1

   do itim=1,2
     do isym=1,nsym
!      ieqqpt is sent onto iqpt by itim/isym
       ieqqpt = qpttoqpt(itim,isym,iqpt)
       gkk_qpt_tmp(:,:,:,:,:) = zero


       if (gkk_flag(1,1,1,1,ieqqpt) == -1) cycle
!      if we have information on this qpt
!      iqpt is equivalent to ieqqpt: get it from file or memory
       do ik_this_proc =1,elph_ds%k_phon%my_nkpt
         ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

         if (elph_ds%gkqwrite == 0) then
           gkk_qpt_tmp(:,:,:,ikpt_phon,:) = elph_ds%gkk_qpt(:,:,:,ik_this_proc,:,ieqqpt)
         else if (elph_ds%gkqwrite == 1) then
           read(elph_ds%unitgkq,REC=((ieqqpt-1)*elph_ds%k_phon%my_nkpt+ik_this_proc)) gkk_qpt_tmp(:,:,:,ikpt_phon,:)
         end if
       end do

!      condense everything
       call xmpi_sum (gkk_qpt_tmp, comm, ierr)

       neqqpt=neqqpt+1

       if (elph_ds%ep_scalprod==1) then
         do ii=1,3
           do jj=1,3
             ss(ii,jj)=0.0_dp
             do kk=1,3
               do ll=1,3
                 ss(ii,jj)=ss(ii,jj)+rprimd(ii,kk)*symrel(kk,ll,isym)*gprimd(ll,jj)
               end do
             end do
           end do
         end do
       else
         do ii=1,3
           do jj=1,3
             ss(ii,jj) = symrec(jj,ii,isym)
           end do
         end do
       end if

       ss_allatoms(:,:,:) = zero
       do iatom=1,natom
         ancestor_iatom = indsym(4,isym,iatom)
!        do jatom=1,natom
!        ancestor_jatom = indsym(4,isym,jatom)
         ss_allatoms(1,(ancestor_iatom-1)*3+1:(ancestor_iatom-1)*3+3,&
&         (iatom-1)*3+1:         (iatom-1)*3+3) = ss(1:3,1:3)
!        end do
       end do


!      NOTE   ssinv(ii,jj)=ssinv(ii,jj)+gprimd(ii,kk)*rprimd(jj,ll)*symrec(ll,kk,isym)

       do isppol=1,elph_ds%nsppol
         do ikpt_phon=1,elph_ds%k_phon%nkpt
!          symikpt_phon is sent onto ikpt_phon by itim/isym
           symikpt_phon=elph_ds%k_phon%full2full(itim,isym,ikpt_phon)

!          Do each element band1, band2 separately...
           do ib1=1,elph_ds%ngkkband*elph_ds%ngkkband

!            multiply by the ss matrices
             tmp_mat2(:,:,:) = zero
             tmp_mat(:,:,:) = reshape(gkk_qpt_tmp(:,ib1,:,ikpt_phon,isppol),&
&             (/2,elph_ds%nbranch,elph_ds%nbranch/))
             call ZGEMM ('N','N',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,&
&             c_one,ss_allatoms,elph_ds%nbranch,tmp_mat,elph_ds%nbranch,c_zero,&
&             tmp_mat2,elph_ds%nbranch)
             call ZGEMM ('N','T',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,&
&             c_one,tmp_mat2,elph_ds%nbranch,ss_allatoms,elph_ds%nbranch,c_zero,&
&             tmp_mat,elph_ds%nbranch)

!            add to gkk_qpt_new
             do ibranch =1,elph_ds%nbranch
               do jbranch =1,elph_ds%nbranch
                 gkk_qpt_new(:,ib1,(jbranch-1)*elph_ds%nbranch+ibranch,symikpt_phon,isppol) = &
&                 gkk_qpt_new(:,ib1,(jbranch-1)*elph_ds%nbranch+ibranch,symikpt_phon,isppol) + &
&                 tmp_mat(:,jbranch,ibranch)
               end do
             end do

           end do ! end ib1 do
         end do ! end ikpt_phon do
       end do ! end isppol do

     end do ! end isym do
   end do ! itim

   if (neqqpt > 1) then
     write(std_out,*) ' found several equiv qpts and am symmetrizing them ', neqqpt
   end if

!  divide by number of equivalent qpts found
   gkk_qpt_new(:,:,:,:,:) = gkk_qpt_new(:,:,:,:,:)/neqqpt

!  copy the symmetrized version into all the equivalent qpoints, appropriately transformed
!  See above
!  itim=1
   do itim=1,2
     do isym=1,nsym
!      ieqqpt is sent onto iqpt by itim/isym
       ieqqpt = qpttoqpt(itim,isym,iqpt)

       if (symmetrized_qpt(ieqqpt) /= -1) cycle
       gkk_qpt_tmp(:,:,:,:,:) = zero

!      use symrec matrices to get inverse transform from isym^{-1}
       if (elph_ds%ep_scalprod==1) then
         do ii=1,3
           do jj=1,3
             ss(ii,jj)=0.0_dp
             do kk=1,3
               do ll=1,3
!                Use inverse of symop matrix here to get back to ieqqpt (inv+transpose is in symrec and in gprimd)
                 ss(ii,jj)=ss(ii,jj)+rprimd(ii,kk)*symrec(ll,kk,isym)*gprimd(ll,jj)
               end do
             end do
           end do
         end do
       else
         do ii=1,3
           do jj=1,3
             ss(ii,jj) = symrel(ii,jj,isym)
           end do
         end do
       end if

       ss_allatoms(:,:,:) = zero
       do iatom=1,natom
         ancestor_iatom = indsym(4,isym,iatom)
!        do jatom=1,natom
!        ancestor_jatom = indsym(4,isym,jatom)
         ss_allatoms(1,(ancestor_iatom-1)*3+1:(ancestor_iatom-1)*3+3,&
&         (iatom-1)*3+1:          (iatom-1)*3+3) = ss(1:3,1:3)
!        end do
       end do

!      ! Use inverse of symop matrix here to get back to ieqqpt
!      ssinv(ii,jj)=ssinv(ii,jj)+gprimd(ii,kk)*rprimd(jj,ll)*symrel(kk,ll,isym)

       do isppol=1,elph_ds%nsppol
         do ikpt_phon=1,elph_ds%k_phon%nkpt
!          symikpt_phon is sent onto ikpt_phon by itim/isym
           symikpt_phon=elph_ds%k_phon%full2full(itim,isym,ikpt_phon)

           do ib1=1,elph_ds%ngkkband*elph_ds%ngkkband

!            multiply by the ss^{-1} matrices
             tmp_mat2(:,:,:) = zero
             tmp_mat(:,:,:) = reshape(gkk_qpt_new(:,ib1,:,ikpt_phon,isppol),&
&             (/2,elph_ds%nbranch,elph_ds%nbranch/))
             call ZGEMM ('N','N',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,&
&             c_one,ss_allatoms,elph_ds%nbranch,tmp_mat,elph_ds%nbranch,c_zero,&
&             tmp_mat2,elph_ds%nbranch)
             call ZGEMM ('N','T',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,&
&             c_one,tmp_mat2,elph_ds%nbranch,ss_allatoms,elph_ds%nbranch,c_zero,&
&             tmp_mat,elph_ds%nbranch)

             do ibranch =1,elph_ds%nbranch
               do jbranch =1,elph_ds%nbranch
                 gkk_qpt_tmp(:,ib1,(jbranch-1)*elph_ds%nbranch+ibranch,symikpt_phon,isppol) =&
&                 tmp_mat(:,jbranch,ibranch)
               end do
             end do

             do ik_this_proc =1,elph_ds%k_phon%my_nkpt
               if (elph_ds%k_phon%my_ikpt(ik_this_proc) == symikpt_phon) then
                 if (gkk_flag (1,1,ik_this_proc,isppol,ieqqpt) == -1) gkk_flag (:,:,ik_this_proc,isppol,ieqqpt) = 0
                 exit
               end if
             end do
!             if (gkk_flag (1,1,symikpt_phon,isppol,ieqqpt) == -1) then
!               gkk_flag (:,:,symikpt_phon,isppol,ieqqpt) = 0
!             end if

           end do ! end ib1 do
         end do ! end ikpt_phon do
       end do ! end isppol do


!      save symmetrized matrices for qpt ieqqpt
       do ik_this_proc =1,elph_ds%k_phon%my_nkpt
         ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

         if (elph_ds%gkqwrite == 0) then
           elph_ds%gkk_qpt(:,:,:,ik_this_proc,:,ieqqpt) = gkk_qpt_tmp(:,:,:,ikpt_phon,:)
         else if (elph_ds%gkqwrite == 1) then
           write(elph_ds%unitgkq,REC=((ieqqpt-1)*elph_ds%k_phon%my_nkpt+ik_this_proc)) gkk_qpt_tmp(:,:,:,ikpt_phon,:)
         end if
       end do

       symmetrized_qpt(ieqqpt) = 1

     end do ! end isym do
   end do ! end itim do

 end do
!end iqpt do

 ABI_DEALLOCATE(gkk_qpt_new)
 ABI_DEALLOCATE(gkk_qpt_tmp)

end subroutine complete_gkk
!!***

!!****f* ABINIT/get_nv_fs_en
!! NAME
!!  get_nv_fs_en
!!
!! FUNCTION
!! This routine finds the energy grids for the integration on epsilon
!! and epsilon prime. It then calculates the DOS and FS averaged velocity_sq at
!! these energies. Metals and semiconductors are treated differently, to deal
!! correctly with the gap.
!!
!! INPUTS
!! crystal<crystal_t>=data type gathering info on the crystalline structure.
!! Ifc<ifc_type>=Object containing the interatomic force constants.
!!  elph_ds
!!    elph_ds%nband = number of bands in ABINIT
!!    elph_ds%k_fine%nkptirr = Number of irreducible points for which there exist at least one band that crosses the Fermi level.
!!    elph_ds%nbranch = number of phonon branches = 3*natom
!!    elph_ds%k_phon%nkpt = number of k points
!!    elph_ds%k_fine%irredtoGS = mapping of elph k-points to ground state grid
!!    elph_ds%minFSband = lowest band included in the FS integration
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%fermie = fermi energy
!!    elph_ds%tempermin = minimum temperature at which resistivity etc are calculated (in K)
!!    elph_ds%temperinc = interval temperature grid on which resistivity etc are calculated (in K)
!!    elph_ds%ep_b_min= first band taken into account in FS integration (if telphint==2)
!!    elph_ds%ep_b_max= last band taken into account in FS integration (if telphint==2)
!!    elph_ds%telphint = flag for integration over the FS with 0=tetrahedra 1=gaussians
!!    elph_ds%elphsmear = smearing width for gaussian integration
!!           or buffer in energy for calculations with tetrahedra (telphint=0)
!!
!!  elph_tr_ds
!!    elph_tr_ds%el_veloc = electronic velocities from the fine k-grid
!!
!!  eigenGS = Ground State eigenvalues
!!  kptrlatt_fine = k-point grid vectors (if divided by determinant of present matrix)
!!  max_occ = maximal occupancy for a band
!!
!! OUTPUT
!!  elph_ds%nenergy = number of energy points for integration on epsilon
!!  elph_tr_ds%en_all = energy points
!!  elph_tr_ds%de_all = differences between energy points
!!  elph_tr_ds%dos_n = DOS at selected energy points
!!  elph_tr_ds%veloc_sq = FS averaged velocity square at selected energy points
!!  elph_tr_ds%tmp_gkk_intweight = integration weights at coarse k grid
!!  elph_tr_ds%tmp_velocwtk = velocity times integration weights at coarse k grid
!!  elph_tr_ds%tmp_vvelocwtk = velocity square times integration weights at coarse k grid
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      d2c_weights,ep_el_weights,ep_fs_weights,get_veloc_tr,ifc_fourq,wrtout
!!
!! SOURCE

subroutine get_nv_fs_en(crystal,ifc,elph_ds,eigenGS,max_occ,elph_tr_ds,omega_max)

!Arguments ------------------------------------
!Scalars
 real(dp), intent(in)  :: max_occ
 real(dp), intent(out) :: omega_max
 type(ifc_type),intent(in) :: ifc
 type(crystal_t),intent(in) :: crystal
 type(elph_type),intent(inout) :: elph_ds
 type(elph_tr_type),intent(inout) :: elph_tr_ds
!Arrays

 real(dp), intent(in)  :: eigenGS(elph_ds%nband,elph_ds%k_fine%nkptirr,elph_ds%nsppol)

!Local variables-------------------------------
!scalars
 integer ::  iFSqpt,isppol,ie1,ierr
 integer ::  i_metal,low_T
 integer ::  in_nenergy, out_nenergy
 integer ::  n_edge1, n_edge2, edge
 integer ::  ie_all, ne_all
 integer ::  sz1, sz2, sz3, sz4
  real(dp) :: e_vb_max, e_cb_min,ucvol
 real(dp) :: e1,max_e,fine_range
 real(dp) :: enemin,enemax
 real(dp) :: Temp,e_tiny,de0
 real(dp) :: eff_mass1, eff_mass2, tmp_dos
 character(len=500) :: message
!arrays
 real(dp) :: gprimd(3,3)
 real(dp) :: kpt_2nd(3), e_cb_2nd(2), en1(2)
 real(dp),allocatable :: dos_e1(:,:),tmp_wtk(:,:,:,:)
 real(dp),allocatable :: phfrq(:,:)
 real(dp),allocatable :: displ(:,:,:,:)

! *************************************************************************

 gprimd = crystal%gprimd
 ucvol = crystal%ucvol

 Temp             = elph_ds%tempermin+elph_ds%temperinc
 elph_ds%delta_e  = kb_HaK*Temp ! about 1000 cm^-1/100, no need to be omega_max
 max_e            = elph_ds%nenergy*kb_HaK*Temp
 e_tiny           = kb_HaK*0.00001_dp ! this is the min. delta_e
 de0              = kb_HaK*Temp ! Kb*T

 in_nenergy = elph_ds%nenergy

 ABI_ALLOCATE(tmp_wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt,elph_ds%nsppol,4))
 ABI_ALLOCATE(dos_e1,(elph_ds%nsppol,3))

 ABI_ALLOCATE(phfrq,(elph_ds%nbranch, elph_ds%k_phon%nkpt))
 ABI_ALLOCATE(displ,(2, elph_ds%nbranch, elph_ds%nbranch, elph_ds%k_phon%nkpt))

 do iFSqpt=1,elph_ds%k_phon%nkpt
   call ifc%fourq(crystal,elph_ds%k_phon%kpt(:,iFSqpt),phfrq(:,iFSqpt),displ(:,:,:,iFSqpt))
 end do

 omega_max = maxval(phfrq)*1.1_dp
 ABI_DEALLOCATE(phfrq)
 ABI_DEALLOCATE(displ)

 write(message,'(a,E20.12)')' The max phonon energy is  ', omega_max
 call wrtout(std_out,message,'COLL')

 enemin = elph_ds%fermie - max_e*2
 enemax = elph_ds%fermie + max_e
 call ep_el_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
& enemin, enemax, 4, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
& elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
& elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine, tmp_wtk)

 do isppol=1,elph_ds%nsppol
   dos_e1(isppol,1) = sum(tmp_wtk(:,:,isppol,2))/elph_ds%k_fine%nkpt
   dos_e1(isppol,2) = sum(tmp_wtk(:,:,isppol,3))/elph_ds%k_fine%nkpt
   dos_e1(isppol,3) = sum(tmp_wtk(:,:,isppol,4))/elph_ds%k_fine%nkpt

!  ! BXU, only treat metallic case at this moment, as variational method may not
!  ! apply to insulators
!  i_metal = -1
   i_metal = 1
!  if (dos_e1(isppol,1) .gt. 0.1_dp .and. dos_e1(isppol,2) .gt. 0.1_dp .and. &
!  &   dos_e1(isppol,3) .gt. 0.1_dp) then ! metal
!  i_metal = 1
   if (i_metal == 1) then
     write(message,'(a)')' This is a metal.'
     call wrtout(std_out,message,'COLL')

     fine_range = 1.5_dp
     e1 = elph_ds%fermie + omega_max*fine_range
     out_nenergy = 0
     low_T = 1
     if (omega_max*fine_range .lt. max_e) then
       low_T = 0
       de0 = omega_max*fine_range/in_nenergy ! energy spacing within Ef +/- omega_max
       do while ((e1-elph_ds%fermie) .lt. max_e)
         e1 = e1 + elph_ds%delta_e
         out_nenergy = out_nenergy + 1
       end do
     end if

     if (low_T == 0) max_e = e1 - elph_ds%fermie
     elph_ds%nenergy = in_nenergy*2 + 1 + out_nenergy*2

   else ! semiconductor/insulator, need careful consideration later
     i_metal = 0
!    between CB min and the next k point, use free electron to replace
!    The weights will be proportional to the DOS, relative to the weights
!    calculated with ep_fs_weights, tetrahedron method prefered

!    output VB and CB edges for semiconductor/insulator
     e_vb_max = maxval(eigenGS(elph_ds%minFSband+elph_ds%nFSband/2-1,:,isppol))
     e_cb_min = minval(eigenGS(elph_ds%minFSband+elph_ds%nFSband/2,:,isppol))
     e_cb_2nd(1) = eigenGS(elph_ds%minFSband+elph_ds%nFSband/2,2,isppol)
     e_cb_2nd(2) = eigenGS(elph_ds%minFSband+elph_ds%nFSband/2+1,2,isppol)
     write(message,'(a,E20.12,2x,E20.12)')' elphon : top of VB, bottom of CB = ',&
&     e_vb_max, e_cb_min
     call wrtout(std_out,message,'COLL')
     write(message,'(a,E20.12)')' elphon : energy at the neighbor kpt = ',e_cb_2nd(1)
     call wrtout(std_out,message,'COLL')

     n_edge1 = 4 ! at the very edge
     n_edge2 = 8  ! sparse to the end of free-electron part

     kpt_2nd(:) = gprimd(:,1)*elph_ds%k_fine%kptirr(1,2) + &
&     gprimd(:,2)*elph_ds%k_fine%kptirr(2,2) + &
&     gprimd(:,3)*elph_ds%k_fine%kptirr(3,2)
     write(message,'(a,3E20.12)')' The neighbor k point is:  ', elph_ds%k_fine%kptirr(:,2)
     call wrtout(std_out,message,'COLL')

     if (dabs(elph_ds%fermie-e_cb_min) .lt. dabs(elph_ds%fermie-e_vb_max)) then
       e1 = e_cb_2nd(1)
     else
       e1 = e_vb_max
     end if
     call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&     e1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&     elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&     elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)

     elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt

     eff_mass1 = (kpt_2nd(1)*kpt_2nd(1) + kpt_2nd(2)*kpt_2nd(2) + kpt_2nd(3)*kpt_2nd(3)) / &
&     (2.0_dp*(e_cb_2nd(1)-e_cb_min))
     write(message,'(a,E20.12)')' The eff. mass from band1 is: ', eff_mass1
     call wrtout(std_out,message,'COLL')
     eff_mass2 = (kpt_2nd(1)*kpt_2nd(1) + kpt_2nd(2)*kpt_2nd(2) + kpt_2nd(3)*kpt_2nd(3)) / &
&     (2.0_dp*(e_cb_2nd(2)-e_cb_min))
     write(message,'(a,E20.12)')' The eff. mass from band2 is: ', eff_mass2
     call wrtout(std_out,message,'COLL')

!    bxu, but the eff. mass estimated in this way is too small
!    The following is obtained by roughly fitting to the DOS of 48x48x48
     eff_mass1 = 0.91036
     write(message,'(a,E20.12)')' The eff. mass we are using is: ', eff_mass1
     call wrtout(std_out,message,'COLL')

     tmp_dos = (ucvol/2.0_dp/pi**2.0_dp)*(2.0_dp*eff_mass1)**1.5_dp*(e1-e_cb_min)**0.5_dp + &
&     2.0_dp*(ucvol/2.0_dp/pi**2.0_dp)*(2.0_dp*eff_mass2)**1.5_dp*(e1-e_cb_min)**0.5_dp
     write(message,'(a,E20.12)')' The fake DOS at kpt1 =   ', tmp_dos
     call wrtout(std_out,message,'COLL')
     write(message,'(a,E20.12)')' The calculated DOS at kpt1 =   ', elph_ds%n0(isppol)
     call wrtout(std_out,message,'COLL')


     e1 = elph_ds%fermie - max_e
     ie_all = 1
     ne_all = 0
     edge = 0

     call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&     e1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&     elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&     elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)

     elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt
     do while ((e1-elph_ds%fermie) .lt. max_e)
       if (e1 .lt. e_cb_min .and. elph_ds%n0(isppol) .lt. tol9) then
         e1 = e_cb_2nd(1)
         edge = 1
         e1 = e1 + de0
       end if

       if (e1 .lt. e_cb_2nd(1)) then
         e1 = e_cb_2nd(1)
         edge = 1
         e1 = e1 + de0
       end if

       if (e1 .gt. e_cb_2nd(1)) then
         call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&         e1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&         elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&         elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)

         elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt

         e1 = e1 + de0
         ie_all = ie_all + 1
       end if
     end do ! e_all
     ne_all = ie_all - 1 + (n_edge1 + n_edge2 - 1)*edge ! energy levels in the free-electron range
     write(message,'(a,i3,a,i3,a)')' For spin', isppol, '  there are ', &
&     ne_all, '  energy levels considered '
     call wrtout(std_out,message,'COLL')

     elph_ds%nenergy = ne_all
   end if ! metal or insulator
 end do ! isppol

 ABI_DEALLOCATE(tmp_wtk)

 if (elph_ds%nenergy .lt. 2) then
   MSG_ERROR('There are too few energy levels for non-LOVA')
 end if

 sz1=elph_ds%ngkkband;sz2=elph_ds%k_phon%nkpt
 sz3=elph_ds%nsppol;sz4=elph_ds%nenergy+1
 ABI_ALLOCATE(elph_tr_ds%dos_n,(sz4,sz3))
 ABI_ALLOCATE(elph_tr_ds%veloc_sq,(3,sz3,sz4))
 ABI_ALLOCATE(elph_tr_ds%en_all,(sz3,sz4))
 ABI_ALLOCATE(elph_tr_ds%de_all,(sz3,sz4+1))
 ABI_ALLOCATE(elph_tr_ds%tmp_gkk_intweight,(sz1,sz2,sz3,sz4))
 ABI_ALLOCATE(elph_tr_ds%tmp_velocwtk,(sz1,sz2,3,sz3,sz4))
 ABI_ALLOCATE(elph_tr_ds%tmp_vvelocwtk,(sz1,sz2,3,3,sz3,sz4))

 elph_tr_ds%dos_n = zero
 elph_tr_ds%veloc_sq = zero
 elph_tr_ds%tmp_gkk_intweight = zero
 elph_tr_ds%tmp_velocwtk = zero
 elph_tr_ds%tmp_vvelocwtk = zero

 ABI_MALLOC_OR_DIE(elph_ds%k_phon%velocwtk,(elph_ds%nFSband,elph_ds%k_phon%nkpt,3,elph_ds%nsppol), ierr)

 ABI_MALLOC_OR_DIE(elph_ds%k_phon%vvelocwtk,(elph_ds%nFSband,elph_ds%k_phon%nkpt,3,3,elph_ds%nsppol), ierr)

 elph_ds%k_phon%velocwtk = zero
 elph_ds%k_phon%vvelocwtk = zero

!metal
 if (i_metal .eq. 1) then
   e1 = elph_ds%fermie - max_e
   en1(:) = elph_ds%fermie - max_e
   if (low_T .eq. 1) then
     enemin = elph_ds%fermie - max_e - elph_ds%delta_e
     enemax = elph_ds%fermie + max_e

     ABI_ALLOCATE(tmp_wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt,elph_ds%nsppol,elph_ds%nenergy+1))
     call ep_el_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&     enemin, enemax, elph_ds%nenergy+1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&     elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&     elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine, tmp_wtk)

     do isppol=1,elph_ds%nsppol
       do ie1 = 1, elph_ds%nenergy
         elph_tr_ds%en_all(isppol,ie1) = en1(isppol)
         elph_tr_ds%de_all(isppol,ie1) = elph_ds%delta_e

         elph_ds%k_fine%wtk(:,:,isppol) = tmp_wtk(:,:,isppol,ie1+1)
         elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt ! for get_veloc_tr
         elph_tr_ds%dos_n(ie1,isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt

         call get_veloc_tr(elph_ds,elph_tr_ds)
         elph_tr_ds%veloc_sq(:,isppol,ie1)=elph_tr_ds%FSelecveloc_sq(:,isppol)

         call d2c_weights(elph_ds,elph_tr_ds)

         elph_tr_ds%tmp_gkk_intweight(:,:,isppol,ie1) = elph_ds%k_phon%wtk(:,:,isppol)
         elph_tr_ds%tmp_velocwtk(:,:,:,isppol,ie1) = elph_ds%k_phon%velocwtk(:,:,:,isppol)
         elph_tr_ds%tmp_vvelocwtk(:,:,:,:,isppol,ie1) = elph_ds%k_phon%vvelocwtk(:,:,:,:,isppol)
         en1(isppol) = en1(isppol) + elph_ds%delta_e
       end do
     end do
     ABI_DEALLOCATE(tmp_wtk)

   else ! low_T = 0
     enemin = e1 - elph_ds%delta_e
     enemax = e1 + (out_nenergy-1)*elph_ds%delta_e

     ABI_ALLOCATE(tmp_wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt,elph_ds%nsppol,out_nenergy+1))
     call ep_el_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&     enemin, enemax, out_nenergy+1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&     elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&     elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine, tmp_wtk)
     do isppol=1,elph_ds%nsppol
       do ie1 = 1, out_nenergy
         elph_tr_ds%en_all(isppol,ie1) = en1(isppol)
         elph_tr_ds%de_all(isppol,ie1) = elph_ds%delta_e

         elph_ds%k_fine%wtk(:,:,isppol) = tmp_wtk(:,:,isppol,ie1+1)
         elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt ! for get_veloc_tr
         elph_tr_ds%dos_n(ie1,isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt

         call get_veloc_tr(elph_ds,elph_tr_ds)
         elph_tr_ds%veloc_sq(:,isppol,ie1)=elph_tr_ds%FSelecveloc_sq(:,isppol)

         call d2c_weights(elph_ds,elph_tr_ds)

         elph_tr_ds%tmp_gkk_intweight(:,:,isppol,ie1) = elph_ds%k_phon%wtk(:,:,isppol)
         elph_tr_ds%tmp_velocwtk(:,:,:,isppol,ie1) = elph_ds%k_phon%velocwtk(:,:,:,isppol)
         elph_tr_ds%tmp_vvelocwtk(:,:,:,:,isppol,ie1) = elph_ds%k_phon%vvelocwtk(:,:,:,:,isppol)

         en1(isppol) = en1(isppol) + elph_ds%delta_e
       end do
     end do
     ABI_DEALLOCATE(tmp_wtk)

     e1 = en1(1)
     enemin = e1 - de0
     enemax = e1 + in_nenergy*2*de0

     ABI_ALLOCATE(tmp_wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt,elph_ds%nsppol,in_nenergy*2+2))
     call ep_el_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&     enemin, enemax, in_nenergy*2+2, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&     elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&     elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine, tmp_wtk)

     do isppol=1,elph_ds%nsppol
       do ie1 = out_nenergy+1, out_nenergy+in_nenergy*2+1
         elph_tr_ds%en_all(isppol,ie1) = en1(isppol)
         elph_tr_ds%de_all(isppol,ie1) = de0

         elph_ds%k_fine%wtk(:,:,isppol) = tmp_wtk(:,:,isppol,ie1-out_nenergy+1)
         elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt ! for get_veloc_tr
         elph_tr_ds%dos_n(ie1,isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt

         call get_veloc_tr(elph_ds,elph_tr_ds)
         elph_tr_ds%veloc_sq(:,isppol,ie1)=elph_tr_ds%FSelecveloc_sq(:,isppol)

         call d2c_weights(elph_ds,elph_tr_ds)

         elph_tr_ds%tmp_gkk_intweight(:,:,isppol,ie1) = elph_ds%k_phon%wtk(:,:,isppol)
         elph_tr_ds%tmp_velocwtk(:,:,:,isppol,ie1) = elph_ds%k_phon%velocwtk(:,:,:,isppol)
         elph_tr_ds%tmp_vvelocwtk(:,:,:,:,isppol,ie1) = elph_ds%k_phon%vvelocwtk(:,:,:,:,isppol)

         en1(isppol) = en1(isppol) + de0
       end do
     end do
     ABI_DEALLOCATE(tmp_wtk)

     e1 = en1(1)
     enemin = e1 - elph_ds%delta_e
     enemax = e1 + (out_nenergy-1)*elph_ds%delta_e

     ABI_ALLOCATE(tmp_wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt,elph_ds%nsppol,out_nenergy+1))
     call ep_el_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&     enemin, enemax, out_nenergy+1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&     elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&     elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine, tmp_wtk)

     en1(:) = en1(:) - de0 + elph_ds%delta_e ! adjust to make the points symmetric around Ef
     do isppol=1,elph_ds%nsppol
       do ie1 = out_nenergy+in_nenergy*2+2, in_nenergy*2+1+out_nenergy*2
         elph_tr_ds%en_all(isppol,ie1) = en1(isppol)
         elph_tr_ds%de_all(isppol,ie1) = elph_ds%delta_e

         elph_ds%k_fine%wtk(:,:,isppol) = tmp_wtk(:,:,isppol,ie1-out_nenergy-in_nenergy*2)
         elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt ! for get_veloc_tr
         elph_tr_ds%dos_n(ie1,isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt

         call get_veloc_tr(elph_ds,elph_tr_ds)
         elph_tr_ds%veloc_sq(:,isppol,ie1)=elph_tr_ds%FSelecveloc_sq(:,isppol)

         call d2c_weights(elph_ds,elph_tr_ds)

         elph_tr_ds%tmp_gkk_intweight(:,:,isppol,ie1) = elph_ds%k_phon%wtk(:,:,isppol)
         elph_tr_ds%tmp_velocwtk(:,:,:,isppol,ie1) = elph_ds%k_phon%velocwtk(:,:,:,isppol)
         elph_tr_ds%tmp_vvelocwtk(:,:,:,:,isppol,ie1) = elph_ds%k_phon%vvelocwtk(:,:,:,:,isppol)

         en1(isppol) = en1(isppol) + elph_ds%delta_e
       end do
     end do
     ABI_DEALLOCATE(tmp_wtk)
   end if

!semiconductor
 else if (i_metal .eq. 0) then
   e1 = elph_ds%fermie - max_e
   ie_all = 1

   call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&   e1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&   elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&   elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)

   elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt
   do while ((e1-elph_ds%fermie) .lt. max_e)
     if (e1 .lt. e_cb_min .and. elph_ds%n0(isppol) .lt. tol9) then
       e1 = e_cb_min
     end if

     if (ie_all .ge. n_edge1+n_edge2) then
       if (ie_all .eq. n_edge1+n_edge2) e1 = e1 + de0
       call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&       e1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&       elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&       elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)

       elph_tr_ds%tmp_gkk_intweight(:,:,isppol,ie_all) = elph_ds%k_fine%wtk(:,:,isppol)
       elph_tr_ds%dos_n(ie_all,isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt
       elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt ! for get_veloc_tr

       elph_tr_ds%en_all(isppol,ie_all) = e1
       call get_veloc_tr(elph_ds,elph_tr_ds)
       elph_tr_ds%veloc_sq(:,isppol,ie_all)=elph_tr_ds%FSelecveloc_sq(:,isppol)
!      bxu
!      veloc_sq(1,isppol,ie_all) is "1" good and general??

       elph_tr_ds%de_all(isppol,ie_all) = de0
       e1 = e1 + elph_tr_ds%de_all(isppol,ie_all)
       ie_all = ie_all + 1
     else ! divided according to the 1/DOS (evenly)
       if (ie_all .lt. n_edge1) then
         elph_tr_ds%en_all(isppol,ie_all) = e_cb_min + &
&         (e_tiny**(-0.5_dp) - ie_all*(e_tiny**(-0.5_dp)-(e_cb_2nd(1)-e_cb_min)**(-0.5_dp))/ &
&         dble(n_edge1))**(-2.0_dp)
         if (ie_all .gt. 1) then
           elph_tr_ds%de_all(isppol,ie_all) = elph_tr_ds%en_all(isppol,ie_all) - elph_tr_ds%en_all(isppol,ie_all-1)
         else
           elph_tr_ds%de_all(isppol,ie_all) = elph_tr_ds%en_all(isppol,ie_all) - e_cb_min - e_tiny
         end if
         e1 = elph_tr_ds%en_all(isppol,ie_all)
       else
         elph_tr_ds%en_all(isppol,ie_all) = e_cb_min + &
&         ((ie_all-n_edge1+1)/dble(n_edge2))**2.0_dp*(e_cb_2nd(1)-e_cb_min)
         if (ie_all .gt. 1) then
           elph_tr_ds%de_all(isppol,ie_all) = elph_tr_ds%en_all(isppol,ie_all) - elph_tr_ds%en_all(isppol,ie_all-1)
         else
           elph_tr_ds%de_all(isppol,ie_all) = (e_cb_2nd(1)-e_cb_min)/(dble(n_edge2)**2.0_dp)
         end if
         e1 = elph_tr_ds%en_all(isppol,ie_all)
       end if

       call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
&       e1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, max_occ, &
&       elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&       elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)

       elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt ! for get_veloc_tr

       tmp_dos = (ucvol/2.0_dp/pi**2.0_dp)*(2.0_dp*eff_mass1)**1.5_dp*(e1-e_cb_min)**0.5_dp + &
&       2.0_dp*(ucvol/2.0_dp/pi**2.0_dp)*(2.0_dp*eff_mass2)**1.5_dp*(e1-e_cb_min)**0.5_dp
       elph_tr_ds%dos_n(ie_all,isppol) = tmp_dos
       elph_tr_ds%tmp_gkk_intweight(:,:,isppol,ie_all) = elph_ds%k_fine%wtk(:,:,isppol)*tmp_dos/elph_ds%n0(isppol)

       call get_veloc_tr(elph_ds,elph_tr_ds)
       elph_tr_ds%veloc_sq(:,isppol,ie_all)=elph_tr_ds%FSelecveloc_sq(:,isppol)

       if (ie_all .eq. (n_edge1+n_edge2)) e1 = e_cb_2nd(1) + de0
       ie_all = ie_all + 1
     end if
   end do ! ie_all
 else
   MSG_BUG('check i_metal!')
 end if ! metal or insulator

 ABI_DEALLOCATE(dos_e1)

end subroutine get_nv_fs_en
!!***

!!****f* ABINIT/get_nv_fs_temp
!! NAME
!!  get_nv_fs_temp
!!
!! FUNCTION
!! This routine calculates the fermi energy, FD smeared DOS(Ef) and
!! Veloc_sq(Ef) at looped temperatures.
!!
!! INPUTS
!!  elph_ds
!!    elph_ds%nband = number of bands in ABINIT
!!    elph_ds%k_fine%nkptirr = Number of irreducible points for which there exist at least one band that crosses the Fermi level.
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%k_fine%nkpt = number of k points for fine k-grid
!!    elph_ds%k_phon%nkpt = number of k points for coarse k-grid
!!    elph_ds%tempermin = minimum temperature at which resistivity etc are calculated (in K)
!!    elph_ds%temperinc = interval temperature grid on which resistivity etc are calculated (in K)
!!    elph_ds%ep_b_min= first band taken into account in FS integration (if telphint==2)
!!    elph_ds%ep_b_max= last band taken into account in FS integration (if telphint==2)
!!    elph_ds%telphint = flag for integration over the FS with 0=tetrahedra 1=gaussians
!!    elph_ds%elphsmear = smearing width for gaussian integration
!!           or buffer in energy for calculations with tetrahedra (telphint=0)
!!
!!  eigenGS = Ground State eigenvalues
!!  gprimd = reciprocal lattice vectors (dimensionful)
!!  kptrlatt_fine = k-point grid vectors (if divided by determinant of present matrix)
!!  max_occ = maximal occupancy for a band
!!
!! OUTPUT
!!  elph_ds%fermie=Fermi level at input temperature
!!  elph_tr_ds%dos_n0=DOS(Ef) at looped temperatures
!!  elph_tr_ds%veloc_sq0=FS averaged velocity at Ef at looped temperatures
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      ebands_update_occ,ep_fs_weights,get_veloc_tr,wrtout
!!
!! SOURCE

subroutine get_nv_fs_temp(elph_ds,BSt,eigenGS,gprimd,max_occ,elph_tr_ds)

!Arguments ------------------------------------
 type(elph_type),intent(inout) :: elph_ds
 type(ebands_t),intent(inout)   :: BSt
 type(elph_tr_type),intent(inout) :: elph_tr_ds

!Scalars
 real(dp), intent(in) :: max_occ

! arrays
 real(dp), intent(in) :: gprimd(3,3)
 real(dp), intent(in) :: eigenGS(elph_ds%nband,elph_ds%k_fine%nkptirr,elph_ds%nsppol)

!Local variables-------------------------------

 integer :: isppol!, ie1
 integer :: itemp, tmp_nenergy

 character(len=500) :: message

 real(dp) :: Temp, tmp_elphsmear, tmp_delta_e
! real(dp) :: xtr, e1
! real(dp),allocatable :: tmp_wtk(:,:)

! *************************************************************************

 ABI_ALLOCATE(elph_tr_ds%dos_n0,(elph_ds%ntemper,elph_ds%nsppol))
 ABI_ALLOCATE(elph_tr_ds%veloc_sq0,(elph_ds%ntemper,3,elph_ds%nsppol))
!if (elph_ds%use_k_fine == 1) then
!ABI_ALLOCATE(tmp_wtk,(elph_ds%nFSband,elph_ds%k_fine%nkpt))
!else
!ABI_ALLOCATE(tmp_wtk,(elph_ds%nFSband,elph_ds%k_phon%nkpt))
!end if

 elph_tr_ds%dos_n0 = zero
 elph_tr_ds%veloc_sq0 = zero

 tmp_nenergy = 8
 do itemp=1,elph_ds%ntemper  ! runs over temperature in K
   Temp=elph_ds%tempermin + elph_ds%temperinc*dble(itemp)
   tmp_delta_e = kb_HaK*Temp
   Bst%occopt = 3
   Bst%tsmear = Temp*kb_HaK
   tmp_elphsmear = Temp*kb_HaK
   call ebands_update_occ(Bst,-99.99_dp)
   write(message,'(a,f12.6,a,E20.12)')'At T=',Temp,' Fermi level is:',Bst%fermie
   call wrtout(std_out,message,'COLL')
   if (abs(elph_ds%fermie) < tol10) then
     elph_ds%fermie = BSt%fermie
   end if

!  FD smeared DOS and veloc

   call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, tmp_elphsmear, &
&   elph_ds%fermie, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine,&
&   max_occ, elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
&   elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)

   do isppol=1,elph_ds%nsppol
     elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt
     write(message,'(a,f12.6,a,f12.6)')'At T=',Temp,' The DOS at Ef is:', elph_ds%n0(isppol)
     call wrtout(std_out,message,'COLL')

!    For the non-LOVA case, N(Ef) is not that important (canceled out eventually).
!    Should not be important for metal, comment out for now
!    tmp_wtk = zero
!    do ie1=-tmp_nenergy,tmp_nenergy ! use ie1 here, hope there is no confusion
!    e1=Bst%fermie+ie1*tmp_delta_e
!    xtr=(e1-Bst%fermie)/(2.0_dp*kb_HaK*Temp)
!
!    call ep_fs_weights(elph_ds%ep_b_min, elph_ds%ep_b_max, eigenGS, elph_ds%elphsmear, &
!    &       e1, gprimd, elph_ds%k_fine%irredtoGS, elph_ds%kptrlatt_fine, &
!    &       max_occ, elph_ds%minFSband, elph_ds%nband, elph_ds%nFSband, &
!    &       elph_ds%nsppol, elph_ds%telphint, elph_ds%k_fine)
!
!    tmp_wtk(:,:) = tmp_wtk(:,:) + elph_ds%k_fine%wtk(:,:,isppol)* &
!    &       tmp_delta_e/(4.0d0*kb_HaK*Temp)/(COSH(xtr)**2.0d0)
!    end do ! ie1

!    elph_ds%k_fine%wtk(:,:,isppol) = tmp_wtk(:,:)
     elph_tr_ds%dos_n0(itemp,isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt
!    elph_ds%n0(isppol) = sum(elph_ds%k_fine%wtk(:,:,isppol))/elph_ds%k_fine%nkpt ! for get_veloc_tr
!    write(message,'(a,f12.6,a,f12.6)')'At T=',Temp,' The eff. DOS at Ef is:', elph_tr_ds%dos_n0(itemp,isppol)
!    call wrtout(std_out,message,'COLL')
   end do ! isppol
   call get_veloc_tr(elph_ds,elph_tr_ds)
   elph_tr_ds%veloc_sq0(itemp,:,:) = elph_tr_ds%FSelecveloc_sq(:,:)

 end do ! temperature

end subroutine get_nv_fs_temp
!!***

!!****f* ABINIT/get_veloc_tr
!!
!! NAME
!! get_veloc_tr
!!
!! FUNCTION
!!  calculate the (in) and (out) velocity factors for transport
!!
!! INPUTS
!!  elph_ds
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%k_fine%nkpt = number of kpts included in the FS integration
!!    elph_ds%nFSband = number of bands included in the FS integration
!!    elph_ds%minFSband = index of the lowest FS band
!!    elph_ds%nqpt_full  = number of Q pts
!!    elph_ds%nqptirred  = number of irreducible Q pts
!!  to index the GS electronic states :
!!  kphon_full2irr = mapping of full FS kpts to irreducible ones
!!   FSfullpqtofull = mapping of k + q to k
!!   FSirredtoGS = mapping of irreducible kpoints to GS set
!!
!! OUTPUT
!! elph_tr_ds%FSelecveloc_sq = avergae FS electronic velocity
!!
!! PARENTS
!!      elphon,get_nv_fs_en,get_nv_fs_temp
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_veloc_tr(elph_ds,elph_tr_ds)

!Arguments ------------------------------------
!arrays
  type(elph_type),intent(in) :: elph_ds
  type(elph_tr_type),intent(inout) :: elph_tr_ds

!Local variables-------------------------------
  !scalars
  integer :: ikpt_fine
  integer :: ib1,fib1,isppol, ii
  real(dp) :: eta2
  !arrays
  real(dp) :: elvelock(3)

! *********************************************************************

 ABI_CHECK(allocated(elph_tr_ds%FSelecveloc_sq),"FSele not associated")


!precalculate the Fermi speed modulus squared
 elph_tr_ds%FSelecveloc_sq = zero
 do isppol=1,elph_ds%nsppol
   do ikpt_fine=1,elph_ds%k_fine%nkpt
     do ib1=1,elph_ds%nFSband
       fib1=ib1+elph_ds%minFSband-1
       elvelock(:)=elph_tr_ds%el_veloc(ikpt_fine,fib1,:,isppol)
       do ii=1, 3
         eta2=elvelock(ii)*elvelock(ii)
         elph_tr_ds%FSelecveloc_sq(ii, isppol)=elph_tr_ds%FSelecveloc_sq(ii, isppol)&
&         +eta2*elph_ds%k_fine%wtk(ib1,ikpt_fine,isppol)
       end do
     end do
   end do
   elph_tr_ds%FSelecveloc_sq(:,isppol) = elph_tr_ds%FSelecveloc_sq(:,isppol)/elph_ds%k_fine%nkpt/elph_ds%n0(isppol)
!  for factor 1/elph_ds%n0(isppol) see eq 12 of Allen prb 17 3725 [[cite:Allen1978]] : sum of v**2 over all k gives n0 times FSelecveloc_sq
 end do ! end isppol
 write (std_out,*) '  get_veloc_tr: FSelecveloc_sq ', elph_tr_ds%FSelecveloc_sq

 write (std_out,*) 'out of get_veloc_tr'

end subroutine get_veloc_tr
!!***

!!****f* ABINIT/integrate_gamma
!!
!! NAME
!! integrate_gamma
!!
!! FUNCTION
!! This routine integrates the electron phonon coupling matrix
!! over the kpoints on the fermi surface. A dependency on qpoint
!! remains for gamma_qpt
!!
!! INPUTS
!!   elph_ds = elphon datastructure with data and dimensions
!!      elph_ds%qpt_full = qpoint coordinates
!!      elph_ds%nqptirred = number of irred qpoints
!!      elph_ds%qirredtofull = indexing of the GKK qpoints found
!!   FSfullpqtofull = mapping of k+q to k
!!
!! OUTPUT
!!   elph_ds = modified elph_ds%gamma_qpt and created elph_ds%gamma_rpt
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      get_rank,wrtout,xmpi_sum
!!
!! SOURCE

subroutine integrate_gamma(elph_ds,FSfullpqtofull)

!Arguments ------------------------------------
!scalars
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)

!Local variables-------------------------------
!scalars
 integer :: comm,ikpt_phon,ikpt_phonq,ib1,ib2,ibeff,iqpt,iqpt_fullbz,isppol,ierr
 integer :: irec, symrankkpt_phon,nbranch,nsppol,ngkkband, ik_this_proc
 character(len=500) :: message
 character(len=fnlen) :: fname
!arrays
 real(dp),allocatable :: tmp_gkk(:,:,:,:)

! *************************************************************************

 comm = xmpi_world

 write (message,'(3a)')ch10,' entering integrate_gamma ',ch10
 call wrtout(std_out,message,'COLL')

 nsppol   = elph_ds%nsppol
 nbranch  = elph_ds%nbranch
 ngkkband = elph_ds%ngkkband

 ABI_ALLOCATE(elph_ds%gamma_qpt,(2,nbranch**2,nsppol,elph_ds%nqpt_full))
 elph_ds%gamma_qpt = zero

 ABI_ALLOCATE(tmp_gkk ,(2,ngkkband**2,nbranch**2,nsppol))

 if (elph_ds%gkqwrite == 0) then
   call wrtout(std_out,' integrate_gamma : keeping gamma matrices in memory','COLL')
 else if (elph_ds%gkqwrite == 1) then
   fname=trim(elph_ds%elph_base_name) // '_GKKQ'
   write (message,'(2a)')' integrate_gamma : reading gamma matrices from file ',trim(fname)
   call wrtout(std_out,message,'COLL')
 else
   write (message,'(a,i0)')' Wrong value for gkqwrite = ',elph_ds%gkqwrite
   MSG_BUG(message)
 end if



 do iqpt=1,elph_ds%nqptirred
   iqpt_fullbz = elph_ds%qirredtofull(iqpt)
   symrankkpt_phon = elph_ds%k_phon%krank%get_rank (elph_ds%k_phon%kpt(:,iqpt_fullbz))
   write (std_out,*) ' iqpt_fullbz in qpt grid only,  rank ', iqpt_fullbz, symrankkpt_phon

   do ik_this_proc =1,elph_ds%k_phon%my_nkpt
     ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

     if (elph_ds%gkqwrite == 0) then
       tmp_gkk = elph_ds%gkk_qpt(:,:,:,ik_this_proc,:,iqpt)
     else if (elph_ds%gkqwrite == 1) then
       irec = (iqpt-1)*elph_ds%k_phon%my_nkpt+ik_this_proc
       if (ikpt_phon == 1) then
         write (std_out,*) ' integrate_gamma  read record ', irec
       end if
       read (elph_ds%unitgkq,REC=irec) tmp_gkk(:,:,:,:)
     end if

     do isppol=1,nsppol
       ikpt_phonq = FSfullpqtofull(ikpt_phon,iqpt_fullbz)
!
       do ib1=1,ngkkband
         do ib2=1,ngkkband
           ibeff = ib2+(ib1-1)*ngkkband
           elph_ds%gamma_qpt(:,:,isppol,iqpt_fullbz) = elph_ds%gamma_qpt(:,:,isppol,iqpt_fullbz) + &
&           tmp_gkk(:,ibeff,:,isppol)&
&           *elph_ds%gkk_intweight(ib1,ikpt_phon,isppol)*elph_ds%gkk_intweight(ib2,ikpt_phonq,isppol)
!          NOTE: if ngkkband==1 we are using trivial weights since average
!          over bands was done in normsq_gkk (nmsq_gam_sumFS or nmsq_pure_gkk)
         end do ! ib2
       end do ! ib1
     end do ! isppol
   end do ! ikpt_phon
 end do ! iqpt

 call xmpi_sum (elph_ds%gamma_qpt, comm, ierr)

 ABI_DEALLOCATE(tmp_gkk)

!need prefactor of 1/nkpt for each integration over 1 kpoint index. NOT INCLUDED IN elph_ds%gkk_intweight
 do iqpt=1,elph_ds%nqptirred
   iqpt_fullbz = elph_ds%qirredtofull(iqpt)
!  elph_ds%gamma_qpt(:,:,:,iqpt_fullbz) = elph_ds%gamma_qpt(:,:,:,iqpt_fullbz) / elph_ds%k_phon%nkpt / n0(1) / n0(1)
!  elph_ds%gamma_qpt(:,:,:,iqpt_fullbz) = elph_ds%gamma_qpt(:,:,:,iqpt_fullbz) / elph_ds%k_phon%nkpt / elph_ds%k_phon%nkpt
   elph_ds%gamma_qpt(:,:,:,iqpt_fullbz) = elph_ds%gamma_qpt(:,:,:,iqpt_fullbz) * elph_ds%occ_factor / elph_ds%k_phon%nkpt
 end do

 call wrtout(std_out,' integrate_gamma: gamma matrices have been calculated for recip space and irred qpoints ',"COLL")

end subroutine integrate_gamma
!!***

!!****f* ABINIT/integrate_gamma_tr
!!
!! NAME
!! integrate_gamma_tr
!!
!! FUNCTION
!! This routine integrates the TRANSPORT electron phonon coupling matrices
!! over the kpoints on the fermi surface. A dependency on qpoint
!! remains for gamma_qpt_in/out
!! Copied from integrate_gamma
!!
!! INPUTS
!!   elph_ds = elphon datastructure with data and dimensions
!!      elph_ds%qpt_full = qpoint coordinates
!!   FSfullpqtofull = mapping of k+q to k
!!   veloc_sq1 = mean square electronic velocity on constant energy surface
!!   veloc_sq2 = mean square electronic velocity on constant energy surface
!!
!! OUTPUT
!!   elph_tr_ds%gamma_qpt_tr and created elph_tr_ds%gamma_rpt_tr
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      wrtout,xmpi_sum
!!
!! SOURCE

subroutine integrate_gamma_tr(elph_ds,FSfullpqtofull,s1,s2, veloc_sq1,veloc_sq2,elph_tr_ds)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: s1,s2
 type(elph_tr_type), intent(inout) :: elph_tr_ds
 type(elph_type),intent(in) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)
 real(dp),intent(in) :: veloc_sq1(3,elph_ds%nsppol), veloc_sq2(3,elph_ds%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,ikpt_phonq,ib1,ib2,ibeff,ierr,iqpt,iqpt_fullbz,isppol
 integer :: itensor, icomp, jcomp,comm
 integer :: fib1, fib2
 integer :: ik_this_proc
! integer :: ikpttemp
 character(len=500) :: message
 real(dp) :: wtk, wtkpq, interm
 real(dp) :: veloc1_i, veloc1_j, veloc2_i, veloc2_j
!arrays
 real(dp) :: elvelock(3), elvelockpq(3)
 real(dp) :: velocwtk(3), velocwtkpq(3)
 real(dp) :: vvelocwtk(3,3), vvelocwtkpq(3,3)
 real(dp),allocatable :: tmp_gkk(:,:,:,:)

! *************************************************************************

 comm = xmpi_world

!information
 if (elph_ds%gkqwrite == 0) then
   write (message,'(a)')' integrate_gamma_tr : keeping gamma matrices in memory'
   call wrtout(std_out,message,'COLL')
 else if (elph_ds%gkqwrite == 1) then
   write (message,'(a)')' integrate_gamma_tr : reading gamma matrices from disk'
   call wrtout(std_out,message,'COLL')
 else
   write (message,'(3a,i3)')' integrate_gamma_tr : BUG-',ch10,&
&   ' Wrong value for gkqwrite = ',elph_ds%gkqwrite
   MSG_BUG(message)
 end if

!allocate temp variables
 ABI_MALLOC_OR_DIE(tmp_gkk,(2,elph_ds%ngkkband**2,elph_ds%nbranch**2,elph_ds%nsppol), ierr)

 do iqpt=1,elph_ds%nqptirred
   iqpt_fullbz = elph_ds%qirredtofull(iqpt)
!  write(std_out,*)'iqpt, iqptfullbz  ',iqpt, iqpt_fullbz

   do ik_this_proc =1,elph_ds%k_phon%my_nkpt
     ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

     if (elph_ds%gkqwrite == 0) then
       tmp_gkk = elph_ds%gkk_qpt(:,:,:,ik_this_proc,:,iqpt)
     else if (elph_ds%gkqwrite == 1) then
       read(elph_ds%unitgkq,REC=((iqpt-1)*elph_ds%k_phon%my_nkpt+ik_this_proc)) tmp_gkk
     end if

     ikpt_phonq = FSfullpqtofull(ikpt_phon,iqpt_fullbz)

     do isppol=1,elph_ds%nsppol
       do ib1=1,elph_ds%ngkkband !FS bands
         fib1=ib1+elph_ds%minFSband-1 ! full bands
         elvelock(:)=elph_tr_ds%el_veloc(ikpt_phon,fib1,:,isppol)
         wtk=elph_tr_ds%tmp_gkk_intweight1(ib1,ikpt_phon,isppol)
         velocwtk(:)=elph_tr_ds%tmp_velocwtk1(ib1,ikpt_phon,:,isppol)
         vvelocwtk(:,:)=elph_tr_ds%tmp_vvelocwtk1(ib1,ikpt_phon,:,:,isppol)

         do ib2=1,elph_ds%ngkkband ! FS bands
           ibeff=ib2+(ib1-1)*elph_ds%ngkkband ! full bands
           fib2=ib2+elph_ds%minFSband-1
           elvelockpq(:)= elph_tr_ds%el_veloc(ikpt_phonq,fib2,:,isppol)
           wtkpq=elph_tr_ds%tmp_gkk_intweight2(ib2,ikpt_phonq,isppol)
           velocwtkpq(:)=elph_tr_ds%tmp_velocwtk2(ib2,ikpt_phonq,:,isppol)
           vvelocwtkpq(:,:)=elph_tr_ds%tmp_vvelocwtk2(ib2,ikpt_phonq,:,:,isppol)

!          MJV 31/03/2009: Note that the following is valid for any geometry, not just cubic!
!          see eq 5 and 6 of prb 36 4103 (Al-Lehaibi et al 1987) [[cite:Al-Lehaibi1987]],
!          see also Allen PRB 17 3725 [[cite:Allen1978]]
!          generalization to tensorial quantities is simple, by keeping the directional
!          references of velock and velockpq as indices.
           do icomp = 1, 3
             do jcomp = 1, 3
               itensor = (icomp-1)*3+jcomp
!              FIXME: could use symmetry i <-> j

               veloc1_i = sqrt(veloc_sq1(icomp,isppol))
               veloc1_j = sqrt(veloc_sq1(jcomp,isppol))
               veloc2_i = sqrt(veloc_sq2(icomp,isppol))
               veloc2_j = sqrt(veloc_sq2(jcomp,isppol))
               if (elph_ds%use_k_fine == 1) then
                 interm = vvelocwtk(icomp,jcomp)*wtkpq/veloc1_i/veloc1_j + &
&                 s1*s2*vvelocwtkpq(icomp,jcomp)*wtk/veloc2_i/veloc2_j - &
&                 s1*velocwtk(jcomp)*velocwtkpq(icomp)/veloc1_j/veloc2_i - &
&                 s2*velocwtk(icomp)*velocwtkpq(jcomp)/veloc1_i/veloc2_j

                 elph_tr_ds%gamma_qpt_tr(:,itensor,:,isppol,iqpt_fullbz) = &
&                 elph_tr_ds%gamma_qpt_tr(:,itensor,:,isppol,iqpt_fullbz) + &
&                 tmp_gkk(:,ibeff,:,isppol)*interm
               else
                 elph_tr_ds%gamma_qpt_tr(:,itensor,:,isppol,iqpt_fullbz) = &
&                 elph_tr_ds%gamma_qpt_tr(:,itensor,:,isppol,iqpt_fullbz) + &
&                 tmp_gkk(:,ibeff,:,isppol) &
&                 *(elvelock(icomp)/veloc1_i - s1*elvelockpq(icomp)/veloc2_i) &
&                 *(elvelock(jcomp)/veloc1_j - s2*elvelockpq(jcomp)/veloc2_j) &
&                 *wtk*wtkpq
               end if
             end do
           end do

         end do
       end do
     end do ! isppol

   end do ! ik
 end do ! iq

 call xmpi_sum (elph_tr_ds%gamma_qpt_tr, comm, ierr)

 ABI_DEALLOCATE(tmp_gkk)


!need prefactor of 1/nkpt for each integration over 1 kpoint index.
!NOT INCLUDED IN elph_ds%gkk_intweight
!Add a factor of 1/2 for the cross terms of (v-v')(v-v')
 elph_tr_ds%gamma_qpt_tr = elph_tr_ds%gamma_qpt_tr* elph_ds%occ_factor*0.5_dp / elph_ds%k_phon%nkpt

 write (message,'(2a)')' integrate_gamma_tr : transport gamma matrices are calculated ',&
& ' in recip space and for irred qpoints'
!call wrtout(std_out,message,'COLL')

end subroutine integrate_gamma_tr
!!***

!!****f* ABINIT/integrate_gamma_tr_lova
!!
!! NAME
!! integrate_gamma_tr_lova
!!
!! FUNCTION
!! This routine integrates the TRANSPORT electron phonon coupling matrices
!! over the kpoints on the fermi surface. A dependency on qpoint
!! remains for gamma_qpt_in/out
!! Copied from integrate_gamma
!!
!! INPUTS
!!   elph_ds = elphon datastructure with data and dimensions
!!      elph_ds%qpt_full = qpoint coordinates
!!   FSfullpqtofull = mapping of k+q to k
!!
!! OUTPUT
!!   elph_tr_ds%gamma_qpt_trout
!!   elph_tr_ds%gamma_qpt_trin
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      wrtout,xmpi_sum
!!
!! SOURCE

subroutine integrate_gamma_tr_lova(elph_ds,FSfullpqtofull,elph_tr_ds)

!Arguments ------------------------------------
!scalars
 type(elph_tr_type), intent(inout) :: elph_tr_ds
 type(elph_type),intent(in) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%k_phon%nkpt,elph_ds%nqpt_full)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,ikpt_phonq,ib1,ib2,ibeff,ierr,iqpt,iqpt_fullbz,isppol
 integer :: itensor, icomp, jcomp,comm
 integer :: fib1, fib2
 integer :: ik_this_proc
 real(dp) :: etain, etaout
 character(len=500) :: message
!arrays
 real(dp) :: elvelock(3), elvelockpq(3)
 real(dp),allocatable :: tmp_gkk(:,:,:,:)

! *************************************************************************

 comm = xmpi_world

 ib1=elph_ds%nbranch*elph_ds%nbranch ; ib2=elph_ds%nqpt_full
 ABI_MALLOC_OR_DIE(elph_tr_ds%gamma_qpt_trin,(2,9,ib1,elph_ds%nsppol,ib2), ierr)
 elph_tr_ds%gamma_qpt_trin = zero

 ABI_MALLOC_OR_DIE(elph_tr_ds%gamma_qpt_trout,(2,9,ib1,elph_ds%nsppol,ib2), ierr)
 elph_tr_ds%gamma_qpt_trout = zero

!information
 if (elph_ds%gkqwrite == 0) then
   write (message,'(a)')' integrate_gamma_tr : keeping gamma matrices in memory'
   call wrtout(std_out,message,'COLL')
 else if (elph_ds%gkqwrite == 1) then
   write (message,'(a)')' integrate_gamma_tr : reading gamma matrices from disk'
   call wrtout(std_out,message,'COLL')
 else
   write (message,'(3a,i3)')' integrate_gamma_tr : BUG-',ch10,&
&   ' Wrong value for gkqwrite = ',elph_ds%gkqwrite
   MSG_ERROR(message)
 end if

!allocate temp variables
 ABI_MALLOC_OR_DIE(tmp_gkk,(2,elph_ds%ngkkband**2,elph_ds%nbranch**2,elph_ds%nsppol), ierr)

 do iqpt=1,elph_ds%nqptirred
   iqpt_fullbz = elph_ds%qirredtofull(iqpt)
   write(std_out,*)'iqpt, iqptfullbz  ',iqpt, iqpt_fullbz

   do ik_this_proc =1,elph_ds%k_phon%my_nkpt
     ikpt_phon = elph_ds%k_phon%my_ikpt(ik_this_proc)

     if (elph_ds%gkqwrite == 0) then
       tmp_gkk = elph_ds%gkk_qpt(:,:,:,ik_this_proc,:,iqpt)
     else if (elph_ds%gkqwrite == 1) then
       read(elph_ds%unitgkq,REC=((iqpt-1)*elph_ds%k_phon%my_nkpt+ik_this_proc)) tmp_gkk
     end if

     ikpt_phonq = FSfullpqtofull(ikpt_phon,iqpt_fullbz)

     do isppol=1,elph_ds%nsppol
       do ib1=1,elph_ds%ngkkband
         fib1=ib1+elph_ds%minFSband-1
         elvelock(:)=elph_tr_ds%el_veloc(ikpt_phon,fib1,:,isppol)

         do ib2=1,elph_ds%ngkkband
           ibeff=ib2+(ib1-1)*elph_ds%ngkkband
           fib2=ib2+elph_ds%minFSband-1
           elvelockpq(:)= elph_tr_ds%el_veloc(ikpt_phonq,fib2,:,isppol)


!          MJV 31/03/2009: Note that the following is valid for any geometry, not just cubic!
!          see eq 5 and 6 of prb 36 4103 (Al-Lehaibi et al 1987) [[cite:Al-Lehaibi1987]]
!          see also Allen PRB 17 3725 [[cite:Allen1978]]
!          generalization to tensorial quantities is simple, by keeping the directional
!          references of velock and velockpq as indices.
           do icomp = 1, 3
             do jcomp = 1, 3
               itensor = (icomp-1)*3+jcomp
!              FIXME: could use symmetry i <-> j

               etain  = elvelock(icomp)*elvelockpq(jcomp)
               etaout = elvelock(icomp)*elvelock(jcomp)


               elph_tr_ds%gamma_qpt_trin(:,itensor,:,isppol,iqpt_fullbz) = &
&               elph_tr_ds%gamma_qpt_trin(:,itensor,:,isppol,iqpt_fullbz) + &
&               tmp_gkk(:,ibeff,:,isppol) &
&               *etain &
&               *elph_ds%gkk_intweight(ib1,ikpt_phon,isppol)*elph_ds%gkk_intweight(ib2,ikpt_phonq,isppol)

               elph_tr_ds%gamma_qpt_trout(:,itensor,:,isppol,iqpt_fullbz) = &
&               elph_tr_ds%gamma_qpt_trout(:,itensor,:,isppol,iqpt_fullbz) + &
&               tmp_gkk(:,ibeff,:,isppol) &
&               *etaout &
&               *elph_ds%gkk_intweight(ib1,ikpt_phon,isppol)*elph_ds%gkk_intweight(ib2,ikpt_phonq,isppol)

             end do
           end do
         end do
       end do

     end do ! isppol
   end do ! ik

 end do ! iq

 ABI_DEALLOCATE(tmp_gkk)

 call xmpi_sum (elph_tr_ds%gamma_qpt_trout, comm, ierr)
 call xmpi_sum (elph_tr_ds%gamma_qpt_trin, comm, ierr)


!
!normalize tensor with 1/sqrt(v_x**2 * v_y**2)
!
!move the veloc into mka2f_tr_lova, where T dependence is dealt with
!This will cause some slight difference to the results
 if (.true.) then
   do isppol=1, elph_ds%nsppol
     do icomp = 1, 3
       do jcomp = 1, 3
         itensor = (icomp-1)*3+jcomp
         if(abs(elph_tr_ds%FSelecveloc_sq(icomp,isppol))>tol14**2 .and. abs(elph_tr_ds%FSelecveloc_sq(jcomp,isppol))>tol14**2)then
           elph_tr_ds%gamma_qpt_trin(:,itensor,:,isppol,:) = elph_tr_ds%gamma_qpt_trin(:,itensor,:,isppol,:) / &
&           sqrt(elph_tr_ds%FSelecveloc_sq(icomp,isppol)*elph_tr_ds%FSelecveloc_sq(jcomp,isppol))
           elph_tr_ds%gamma_qpt_trout(:,itensor,:,isppol,:) = elph_tr_ds%gamma_qpt_trout(:,itensor,:,isppol,:) / &
&           sqrt(elph_tr_ds%FSelecveloc_sq(icomp,isppol)*elph_tr_ds%FSelecveloc_sq(jcomp,isppol))
         else
!          XG120528 Fixed problem with zero velocity
           elph_tr_ds%gamma_qpt_trin(:,itensor,:,isppol,:)=zero
           elph_tr_ds%gamma_qpt_trout(:,itensor,:,isppol,:)=zero
         end if
       end do
     end do
   end do ! isppol
 end if

!need prefactor of 1/nkpt for each integration over 1 kpoint index.
!NOT INCLUDED IN elph_ds%gkk_intweight
 elph_tr_ds%gamma_qpt_trout = elph_tr_ds%gamma_qpt_trout* elph_ds%occ_factor / elph_ds%k_phon%nkpt
 elph_tr_ds%gamma_qpt_trin  = elph_tr_ds%gamma_qpt_trin * elph_ds%occ_factor / elph_ds%k_phon%nkpt

 write (message,'(2a)')' integrate_gamma_tr : transport gamma matrices are calculated ',&
& ' in recip space and for irred qpoints'
 call wrtout(std_out,message,'COLL')

!DEBUG
!write(std_out,*)' integrate_gamma_tr_lova: end  elph_tr_ds%gamma_qpt_trin(1,9,1,1,1)=',elph_tr_ds%gamma_qpt_trin(1,9,1,1,1)
!ENDDEBUG

end subroutine integrate_gamma_tr_lova
!!***

!!****f* ABINIT/ftgkk
!!
!! NAME
!! ftgkk
!!
!! FUNCTION
!! If qtor=1 (q->r):
!! Generates the Fourier transform of the recip space gkk matrices
!! to obtain the real space ones.
!! If qtor=0 (r->q):
!! Generates the Fourier transform of the real space gkk matrices
!! to obtain the reciprocal space ones.
!!
!! INPUTS
!! gkqwrite = flag to write recip space matrix elements to disk
!! gkrwrite = flag to write real space matrix elements to disk
!! gprim(3,3)= Normalized coordinates in reciprocal space
!! ikpt_phon0 = starting kpt number for forward FT.
!! natom= Number of atoms in the unit cell
!! nkpt_phon= Number of kpoints used for the FS
!! ngkkband = number of bands kept in gkq and gkr matrix elements (=1 or nband)
!! nkpt_used= number of FS kpoints used, starting at ikpt_phon0
!! nqpt= Number of q points in the Brillouin zone
!!           if qtor=0 this number is read in the input file
!! nrpt= Number of R points in the Big Box
!! qtor= ( q to r : see above )
!! rpt(3,nprt)= Canonical coordinates of the R points in the unit cell
!!           These coordinates are normalized (=> * acell(3)!!)
!! qpt_full(3,nqpt)= Reduced coordinates of the q vectors in reciprocal space
!!           if qtor=0 these vectors are read in the input file
!! unit_gkk_rpt = fortran unit for writing real-space matrix elements
!! unitgkq = fortran unit for writing reciprocal-space matrix elements
!! wghatm(natom,natom,nrpt)
!!         = Weights associated to a pair of atoms and to a R vector
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/output
!! gkk_qpt(2,3*natom,nFSband,nFSband,nkpt_used,nqpt)
!!  = gkk matrices in recip space coming from the Derivative Data Base
!! gkk_rpt(2,3*natom,nFSband,nFSband,nkpt_phon,nqpt)
!!  = gkk matrices in real space stored in file unit_gkk_rpt
!!
!! PARENTS
!!      get_all_gkr,interpolate_gkk,test_ftgkk
!!
!! CHILDREN
!!
!! NOTES
!!   copied from ftiaf9.f
!!   recip to real space: real space is forced to disk file unit_gkk_rpt
!!                        recip space depends on gkqwrite and unitgkq
!!   real to recip space: real space is forced to disk file unit_gkk_rpt
!!                        recip space is necessarily in memory in gkk_qpt
!!
!!    real space elements are complex, but could be reduced, as (-r) = (+r)*
!!
!! SOURCE

subroutine ftgkk (wghatm,gkk_qpt,gkk_rpt,gkqwrite,gkrwrite,gprim,ikpt_phon0,&
&                  natom,nkpt_phon,ngkkband,nkpt_used,nqpt,nrpt,nsppol,&
&                  qtor,rpt,qpt_full,unit_gkk_rpt,unitgkq)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: gkqwrite,gkrwrite,ikpt_phon0,nkpt_phon,natom,ngkkband
 integer,intent(in) :: nkpt_used,nqpt,nrpt,nsppol,qtor,unit_gkk_rpt,unitgkq
!arrays
 real(dp),intent(in) :: gprim(3,3),rpt(3,nrpt),qpt_full(3,nqpt)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 real(dp),intent(inout) :: gkk_qpt(2,ngkkband*ngkkband,3*natom*3*natom,nkpt_used,nsppol,nqpt)
 real(dp),intent(inout) :: gkk_rpt(2,ngkkband*ngkkband,3*natom*3*natom,nkpt_used,nsppol,nrpt)

!Local variables -------------------------
!scalars
 integer :: ikpt_phon,iatom,ib1,ieffkpt_phon,ip,iqpt,irpt,isppol
 integer :: jatom
 real(dp) :: im,kr,re
 character(len=500) :: message
!arrays
 real(dp) :: coskr(nqpt,nrpt),ftwght(2,3*natom*3*natom)
 real(dp) :: gkk_qpt_tmp(2,ngkkband*ngkkband,3*natom*3*natom,nkpt_used,nsppol)
 real(dp) :: gkk_rpt_tmp(2,ngkkband*ngkkband,3*natom*3*natom,nkpt_phon,nsppol)
 real(dp) :: kk(3),sinkr(nqpt,nrpt)

! *********************************************************************

!rewind (unit_gkk_rpt)

!prepare the phase factors
 do iqpt=1,nqpt
!  Calculation of the k coordinates in Normalized Reciprocal
!  coordinates
   kk(1)=   qpt_full(1,iqpt)*gprim(1,1)+&
&   qpt_full(2,iqpt)*gprim(1,2)+&
&   qpt_full(3,iqpt)*gprim(1,3)
   kk(2)=   qpt_full(1,iqpt)*gprim(2,1)+&
&   qpt_full(2,iqpt)*gprim(2,2)+&
&   qpt_full(3,iqpt)*gprim(2,3)
   kk(3)=   qpt_full(1,iqpt)*gprim(3,1)+&
&   qpt_full(2,iqpt)*gprim(3,2)+&
&   qpt_full(3,iqpt)*gprim(3,3)
   do irpt=1,nrpt
!    Product of k and r
     kr =        kk(1)*rpt(1,irpt)+&
&     kk(2)*rpt(2,irpt)+&
&     kk(3)*rpt(3,irpt)
     coskr(iqpt,irpt)=cos(two_pi*kr)
     sinkr(iqpt,irpt)=sin(two_pi*kr)
!    DEBUG
!    if (iqpt < 1000 .and. (irpt == 101 .or. irpt == 901)) then
!    write(std_out,*) iqpt,irpt,kk,rpt(:,irpt),coskr(iqpt,irpt), sinkr(iqpt,irpt)
!    end if
!    ENDDEBUG
   end do
 end do



!Recip to real space
 if (qtor==1) then
!
   if (nkpt_used /= nkpt_phon) write(std_out,*) 'ftgkk: strange usage of nkpt_used for back FT!'
   do irpt=1,nrpt
!    DEBUG
!    write(std_out,*) ' ftgkk : G->R irpt = ',irpt,' / ',nrpt
!    ENDDEBUG
     gkk_rpt_tmp(:,:,:,:,:) = zero

     do iqpt=1,nqpt

!      write(std_out,*) iqpt

       if (gkqwrite == 0) then
         gkk_qpt_tmp(:,:,:,:,:) = gkk_qpt(:,:,:,:,:,iqpt)
       else
         do ikpt_phon=1, nkpt_phon
           read(unitgkq,REC=((iqpt-1)*nkpt_phon+ikpt_phon)) gkk_qpt_tmp(:,:,:,ikpt_phon,:)
         end do
       end if
!      Get the phase factor with normalization!
       re=coskr(iqpt,irpt)/nqpt
       im=sinkr(iqpt,irpt)/nqpt
       do isppol=1,nsppol
         do ikpt_phon=1,nkpt_used
!          DEBUG
!          write(std_out,*) ' ftgkk : G->R ikpt_phon = ',ikpt_phon,' / ',nkpt_used
!          ENDDEBUG
           do ip=1,3*natom*3*natom
!            Real and imaginary part of the real-space gkk matrices -> exp(-i k.r)
             do ib1=1,ngkkband*ngkkband
               gkk_rpt_tmp(1,ib1,ip,ikpt_phon,isppol) = gkk_rpt_tmp(1,ib1,ip,ikpt_phon,isppol)&
&               +re*gkk_qpt_tmp(1,ib1,ip,ikpt_phon,isppol) &
&               +im*gkk_qpt_tmp(2,ib1,ip,ikpt_phon,isppol)
               gkk_rpt_tmp(2,ib1,ip,ikpt_phon,isppol) = gkk_rpt_tmp(2,ib1,ip,ikpt_phon,isppol)&
&               +re*gkk_qpt_tmp(2,ib1,ip,ikpt_phon,isppol) &
&               -im*gkk_qpt_tmp(1,ib1,ip,ikpt_phon,isppol)
             end do
           end do
         end do
       end do
     end do
     if (gkrwrite == 0) then
       gkk_rpt(:,:,:,:,:,irpt) = gkk_rpt_tmp(:,:,:,:,:)
     else
       write (unit_gkk_rpt,REC=irpt) gkk_rpt_tmp
     end if
   end do

!  Real space to recip space
 else if (qtor==0) then

!  write(std_out,*) 'ftgkk : shape(gkk_qpt) = ', shape(gkk_qpt)
   gkk_qpt(:,:,:,:,:,:)=zero

!  rewind (unit_gkk_rpt)
   do irpt=1,nrpt
     if (gkrwrite == 0) then
       gkk_rpt_tmp(:,:,:,:,:) = gkk_rpt(:,:,:,:,:,irpt)
     else
       read(unit_gkk_rpt,REC=irpt) gkk_rpt_tmp
     end if


     do iqpt=1,nqpt

!      Avoid recalculating weights nkpt_used*9 times
       do iatom=1,natom
         do jatom=1,natom
           ip = 3*((iatom-1)*natom+jatom-1)
!          copy same weight for all 3 directions
           ftwght(1,ip+1:ip+3)=coskr(iqpt,irpt)*wghatm(iatom,jatom,irpt)
           ftwght(2,ip+1:ip+3)=sinkr(iqpt,irpt)*wghatm(iatom,jatom,irpt)
         end do
       end do



       do ip=1,3*natom*3*natom
!        Get phase factor
         re = ftwght(1,ip)
         im = ftwght(2,ip)

         do isppol=1,nsppol
           do ikpt_phon=1,nkpt_used


!            DEBUG
!            write(std_out,*) ' ftgkk : R->G ikpt_phon = ',ikpt_phon,' / ',nkpt_used
!            ENDDEBUG
!            effective FS kpt in real space array is ikpt_phon+ikpt_phon0-1 to allow for offset
             ieffkpt_phon = ikpt_phon+ikpt_phon0-1
!            write(std_out,*) 'ftgkk :ikpt_phon,iqpt,ieffkpt_phon ', ikpt_phon,iqpt,ieffkpt_phon

             do ib1=1,ngkkband*ngkkband
!              Real and imaginary part of the gamma matrices
               gkk_qpt(1,ib1,ip,ikpt_phon,isppol,iqpt)=&
&               gkk_qpt(1,ib1,ip,ikpt_phon,isppol,iqpt)&
&               +re*gkk_rpt_tmp(1,ib1,ip,ieffkpt_phon,isppol)&
&               -im*gkk_rpt_tmp(2,ib1,ip,ieffkpt_phon,isppol)
!              !DEBUG
               gkk_qpt(2,ib1,ip,ikpt_phon,isppol,iqpt)=&
&               gkk_qpt(2,ib1,ip,ikpt_phon,isppol,iqpt)&
&               +im*gkk_rpt_tmp(1,ib1,ip,ieffkpt_phon,isppol)&
&               +re*gkk_rpt_tmp(2,ib1,ip,ieffkpt_phon,isppol)
!              !ENDDEBUG

!              if (iqpt < 100 .and. irpt < 100 .and. &
!              &   tmpgkkrim(irpt)**2+tmpgkkrre(irpt)**2 > tol6) then
!              write(std_out,'(2I4,2E16.8,x,2E16.8)') &
!              &   iqpt,irpt,re,im,tmpgkkrre(irpt),tmpgkkrim(irpt)
!              end if

             end do
           end do
!          end ikpt_phon
         end do
!        end isppol
!        write(std_out,'(a)') ' ftgkk :gkk_qpt :'
!        write(std_out,'(4E16.5)') gkk_qpt(:,1,1,,ikpt_phon,1:nqpt)
       end do
!      end ip
     end do
!    end iqpt
   end do
!  end irpt


!  There is no other space to Fourier transform from ??
 else
   write(message,'(a,a,a,i0,a)' )&
&   'The only allowed values for qtor are 0 or 1, while',ch10,&
&   'qtor=',qtor,' has been required.'
   MSG_BUG(message)
 end if

end subroutine ftgkk
!!***

!!****f* ABINIT/test_ftgkk
!! NAME
!! test_ftgkk
!!
!! FUNCTION
!!  Test the fourier transform routine ftgkk for the el-phon matrix elements
!!
!! INPUTS
!!   elph_ds = elphon datastructure with matrix elements
!!   gprim = reciprocal lattice vectors
!!   natom = number of atoms
!!   nrpt = number of real space points for FT interpolation
!!   rpt = coordinates of real space points for FT interpolation
!!   qpt_full = qpoint coordinates
!!   wghatm = weights for pairs of atoms in FT interpolation
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!!  MJV 18/5/2008 reverted to old syntax/use for ftgkk, with all ft being done
!!   in a batch. Might come back to 5.5 version with atomic FT in ftgkk, but later.
!!
!! PARENTS
!!
!! CHILDREN
!!      ftgkk
!!
!! SOURCE

subroutine test_ftgkk(elph_ds,gprim,natom,nrpt,rpt,qpt_full,wghatm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrpt
 type(elph_type),intent(inout) :: elph_ds
!arrays
 real(dp),intent(in) :: gprim(3,3),rpt(3,nrpt),qpt_full(3,elph_ds%nqpt_full)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt)

!Local variables-------------------------------
!scalars
 integer :: ikpt_phon,iqpt,isppol,qtor,sz1,sz2
!arrays
 real(dp),allocatable :: gkq_disk(:,:,:,:,:),tmp_gkq(:,:,:,:,:)

! *************************************************************************

!for each qpt do FT to recuperate original values

 isppol = 1
 qtor = 0
 sz1=elph_ds%ngkkband*elph_ds%ngkkband
 sz2=elph_ds%nbranch*elph_ds%nbranch
 ABI_ALLOCATE(gkq_disk,(2,sz1,sz2,elph_ds%k_phon%nkpt,elph_ds%nsppol))
 ABI_ALLOCATE(tmp_gkq,(2,sz1,sz2,elph_ds%k_phon%nkpt,elph_ds%nsppol))

 do iqpt=1,elph_ds%nqpt_full
   tmp_gkq(:,:,:,:,:) = zero

   call ftgkk (wghatm,tmp_gkq,elph_ds%gkk_rpt,elph_ds%gkqwrite,&
&   elph_ds%gkk_rptwrite,gprim,1,natom,&
&   elph_ds%k_phon%nkpt,elph_ds%ngkkband,elph_ds%k_phon%nkpt,1,&
&   nrpt,elph_ds%nsppol,qtor,rpt,qpt_full,elph_ds%unit_gkk_rpt,elph_ds%unitgkq)

   if (elph_ds%gkqwrite == 0) then
     do ikpt_phon=1,10
       write (93,*) tmp_gkq(:,:,:,ikpt_phon,isppol)-elph_ds%gkk_qpt(:,:,:,ikpt_phon,isppol,iqpt)
     end do
   else
     do ikpt_phon=1, elph_ds%k_phon%nkpt
       read (elph_ds%unitgkq,REC=((iqpt-1)*elph_ds%k_phon%nkpt+ikpt_phon)) gkq_disk(:,:,:,ikpt_phon,:)
     end do
     do ikpt_phon=1,10
       write (93,*) tmp_gkq(:,:,:,ikpt_phon,isppol)-gkq_disk(:,:,:,ikpt_phon,isppol)
     end do
   end if
 end do

end subroutine test_ftgkk
!!***

end module m_elphon
!!***
