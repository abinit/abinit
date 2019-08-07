!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_epjdos
!! NAME
!!  m_epjdos
!!
!! FUNCTION
!!  Tools for the computiation of electronic PJDOSes
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (MVer, XG, SM, MT, BAmadon, MG, MB)
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

module m_epjdos

 use defs_basis
 use defs_abitypes
 use m_abicore
 use m_xmpi
 use m_errors
 use m_htetrahedron
 use m_splines
 use m_cgtools
 use m_atomdata
 use m_crystal
 use m_ebands
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr
 use m_mpinfo
 use m_sort

 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_occ,            only : dos_hdr_write
 use m_time,           only : cwtime, timab
 use m_io_tools,       only : open_file
 use m_numeric_tools,  only : simpson, simpson_int
 use m_fstrings,       only : int2char4, strcat
 use m_special_funcs,  only : jlspline_t, jlspline_new, jlspline_free, jlspline_integral
 use m_kpts,           only : tetra_from_kptrlatt
 use m_kg,             only : ph1d3d, getph
 use m_gsphere,        only : getkpgnorm
 use m_fftcore,        only : sphereboundary
 use m_fft,            only : fftpac, fourwf, fourdp
 use m_pawrad,         only : pawrad_type, simp_gen
 use m_pawtab,         only : pawtab_type
 use m_pawcprj,        only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_free
 use m_initylmg,       only : initylmg

 implicit none

 private
!!***

 public :: dos_calcnwrite      ! Calculate DOS and write results to file(s).
 public :: recip_ylm           ! Project input wavefunctions (real space) on to Ylm
 public :: dens_in_sph         ! Calculate integrated density in sphere around each atom
 public :: partial_dos_fractions   ! Calculate partial DOS fractions to feed to the tetrahedron method (PW part)
 public :: partial_dos_fractions_paw ! Calculate PAW contributions to the partial DOS fractions.
!!***

!----------------------------------------------------------------------

!!****t* m_epjdos/epjdos_t
!! NAME
!! epjdos_t
!!
!! FUNCTION
!!  Stores different contributions to the electronic DOS.
!!
!! NOTES
!!  Please contact gmatteo if you plan to change the internal implementation
!!  or add new DOSes. These results are saved in a netcdf file (see fatbands_ncwrite)
!!  so that one can read it with python and plot fatbands and PJDOSEs.
!!  The python version is able to handle the different cases (L, LM, Spin ...) but
!!  any change in the internal Abinit implementation is likely to break the python interface.
!!
!! SOURCE

 type,public :: epjdos_t

   integer :: mbesslang
   ! Max L+1 used in LM-DOS  (Bessel function expansion)

   integer :: ndosfraction
   ! Defines the last dimension of the dos arrays.
   ! Actual value depends on the other variables.

   integer :: prtdos
   ! 2 --> Standard DOS with tetra.
   ! 3 --> L-DOS with tetra (prtdosm>0 if LM is wanted in Ylm/Slm basis).
   ! 4 --> L-DOS with gaussian (prtdosm if LM is wanted in Ylm/Slm basis).
   ! 5 --> Spin-DOS

   integer :: prtdosm
   ! Option for the m-contributions to the partial DOS
   ! 1 if LM-projection is done onto complex Ylm
   ! 2 if LM-projection is done onto real Slm

   integer :: partial_dos_flag

   integer :: paw_dos_flag
   ! 1 if both PAW contributions are evaluated AND stored

   !integer :: pawfatbnd
   integer :: fatbands_flag

   integer :: nkpt, mband, nsppol
   ! Used to dimension arrays

   integer,allocatable :: mlang_type(:)
   ! mlang_type(ntypat + natsph_extra)
   ! Max L+1 used in LM-DOS for each atom type

   real(dp),allocatable :: fractions(:,:,:,:)
   ! fractions(nkpt,mband,nsppol,ndosfraction))
   ! TODO: replace nsppol with nspden = 1, 2 (nsppol==2) or 4 (nspinor==2)

   real(dp),allocatable :: fractions_m(:,:,:,:)
   ! fractions_m(nkpt,mband,nsppol,ndosfraction*mbesslang)

   real(dp),allocatable :: fractions_paw1(:,:,:,:)
   ! fractions_paw1(nkpt,mband,nsppol,ndosfraction)

   real(dp),allocatable :: fractions_pawt1(:,:,:,:)
   ! fractions_pawt1(nkpt,mband,nsppol,ndosfraction))

 end type epjdos_t

 public :: epjdos_new         ! Create new object
 public :: epjdos_free        ! Free dynamic memory

 public :: prtfatbands        ! Print PJDOS contributions in xmgrace format.
 public :: fatbands_ncwrite   ! Write PJDOS contributions to netcdf file.

!----------------------------------------------------------------------

contains  !============================================================
!!***

!!****f* m_epjdos/epjdos_new
!! NAME
!!  epjdos_new
!!
!! FUNCTION
!!  Create new object from dataset input variables.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! PARENTS
!!
!! SOURCE

type(epjdos_t) function epjdos_new(dtset, psps, pawtab) result(new)

!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: ierr,itypat,iat

! *********************************************************************

 new%nkpt = dtset%nkpt; new%mband = dtset%mband; new%nsppol = dtset%nsppol

 new%prtdos = dtset%prtdos
 new%partial_dos_flag = 0
 if (new%prtdos==2) new%partial_dos_flag = 0 ! Standard DOS with tetra.
 if (new%prtdos==3) new%partial_dos_flag = 1 ! L-DOS with tetra (prtdosm>0 if LM is wanted in Ylm/Slm basis).
 if (new%prtdos==4) new%partial_dos_flag = 1 ! L-DOS with gaussian (prtdosm if LM is wanted in Ylm/Slm basis).
 if (new%prtdos==5) new%partial_dos_flag = 2 ! Spin DOS

 new%prtdosm=0
 if (new%partial_dos_flag==1) new%prtdosm=dtset%prtdosm
 ! paw_dos_flag= 1 if both PAW contributions are evaluated AND stored
 new%paw_dos_flag=0
 if (dtset%usepaw==1 .and. new%partial_dos_flag==1 .and. dtset%pawprtdos==1) new%paw_dos_flag=1

 new%fatbands_flag=0
 if (dtset%pawfatbnd>0 .and. new%prtdosm==0) new%fatbands_flag=1
 if (new%prtdosm==1.and.dtset%pawfatbnd>0)then
!  because they compute quantities in real and complex harmonics respectively
   MSG_ERROR('pawfatbnd>0  and prtdosm=1 are not compatible')
 end if

 ! mjv : initialization is needed as mbesslang is used for allocation below
 ! NOTE: 10/5/2010 the whole of this could be looped over ndosfraction,
 ! to store much less in memory. The DOS is accumulated in an array
 ! and then printed to file at the end.
 new%mbesslang = 1
 if (new%partial_dos_flag==1 .or. new%fatbands_flag==1) then

   ABI_MALLOC(new%mlang_type, (dtset%ntypat + dtset%natsph_extra))
   new%mlang_type = 0

   ! TODO: Could use mbesslang = 4 or compute it from psps/pawtab
   ! Increment by one (could underestimate if vloc = vlmax)
   if (dtset%usepaw == 0) then
     do iat=1,dtset%natsph
       itypat = dtset%typat(dtset%iatsph(iat))
       new%mlang_type(itypat) = 1 + maxval(psps%indlmn(1, :, itypat))
     end do
   else
     do iat=1,dtset%natsph
       itypat= dtset%typat(dtset%iatsph(iat))
       new%mlang_type(itypat) = 1 + (pawtab(itypat)%l_size - 1) / 2
     end do
   end if

   ! Up to l=g if we have natsph_extra.
   if (dtset%natsph_extra > 0) new%mlang_type(dtset%ntypat+1:) = 5

   new%mlang_type = 5  ! This is to preserve the old implementation
   new%mbesslang = maxval(new%mlang_type)
   new%ndosfraction = (dtset%natsph + dtset%natsph_extra) * new%mbesslang

 else if (new%partial_dos_flag == 2) then
   new%ndosfraction = 7

 else
   new%ndosfraction = 1
   new%mbesslang = 0
 end if

 ! Check allocations status as these arrays are not distributed and the wavefunctions are still in memory.
 ABI_MALLOC_OR_DIE(new%fractions, (dtset%nkpt,dtset%mband,dtset%nsppol,new%ndosfraction), ierr)
 new%fractions = zero

 if (new%prtdosm>=1 .or. new%fatbands_flag==1) then
   ABI_MALLOC_OR_DIE(new%fractions_m,(dtset%nkpt,dtset%mband,dtset%nsppol,new%ndosfraction*new%mbesslang), ierr)
   new%fractions_m = zero
 end if

 if (dtset%usepaw==1 .and. new%partial_dos_flag==1) then
   ABI_MALLOC_OR_DIE(new%fractions_paw1,(dtset%nkpt,dtset%mband,dtset%nsppol,new%ndosfraction), ierr)
   ABI_MALLOC_OR_DIE(new%fractions_pawt1,(dtset%nkpt,dtset%mband,dtset%nsppol,new%ndosfraction), ierr)
   new%fractions_paw1 = zero; new%fractions_pawt1 = zero
 end if

end function epjdos_new
!!***

!!****f* m_epjdos/epjdos_free
!! NAME
!!  epjdos_free
!!
!! FUNCTION
!!  deallocate memory
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      cwtime,wrtout
!!
!! SOURCE

subroutine epjdos_free(self)

!Arguments ------------------------------------
 type(epjdos_t),intent(inout) :: self

! *********************************************************************

 ! integer
 if (allocated(self%mlang_type)) then
   ABI_FREE(self%mlang_type)
 end if

 ! real
 if (allocated(self%fractions)) then
   ABI_FREE(self%fractions)
 end if
 if (allocated(self%fractions_m)) then
   ABI_FREE(self%fractions_m)
 end if
 if (allocated(self%fractions_paw1)) then
   ABI_FREE(self%fractions_paw1)
 end if
 if (allocated(self%fractions_pawt1)) then
   ABI_FREE(self%fractions_pawt1)
 end if

end subroutine epjdos_free
!!***

!!****f* m_epjdos/dos_calcnwrite
!! NAME
!! dos_calcnwrite
!!
!! FUNCTION
!! calculate DOS and write results to file(s)
!!
!! INPUTS
!!  dos_fractions= projections of wavefunctions on each angular momentum Ylm
!!     which is the weight going into the DOS for an l-decomposed dos
!!  dos_fractions_m= same as dos_fractions, but m-decomposed not just l-
!!  dos_fractions_paw1= contribution to dos fractions from the PAW partial waves (phi)
!!  dos_fractions_pawt1= contribution to dos fractions from the PAW pseudo partial waves (phi_tild)
!!  dtset     structured datatype, in particular one uses :
!!   kptrlatt(3,3)=lattice vectors for full kpoint grid
!!   nshiftk      =number of kpoint grid shifts
!!   pawprtdos    =option to output the individual contributions to the partial DOS (0, 1 or 2)
!!   shiftk(3,nshiftk)=kpoint shifts
!!   usepaw       =option for PAW
!!  crystal<crystal_t>=Object defining the unit cell and its symmetries.
!!  ebands<ebands_t>=Band structure data.
!!  fermie=Fermi energy
!!  fildata=name of the DOS output file
!!  mbesslang=maximum angular momentum for Bessel function expansion
!!  prtdosm=option for the m-contributions to the partial DOS
!!  ndosfraction= number of types of DOS we are calculating, e.g. the number
!!    of l channels. Could be much more general, for other types of partial DOS
!!  paw_dos_flag= option for partial dos in PAW
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  (no explicit output)
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      cwtime,wrtout
!!
!! SOURCE

subroutine dos_calcnwrite(dos,dtset,crystal,ebands,fildata,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 character(len=*),intent(in) :: fildata
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: crystal
 type(ebands_t),intent(in) :: ebands
 type(epjdos_t),intent(in) :: dos

!Local variables-------------------------------
!scalars
 integer,parameter :: bcorr0=0,master=0
 integer :: iat,iband,iene,ikpt,isppol,natsph,natsph_extra,nkpt,nsppol,i1,i2
 integer :: nene,prtdos,unitdos,ierr,prtdosm,paw_dos_flag,mbesslang,ndosfraction
 integer :: my_rank,nprocs,cnt,ifrac,ii
 real(dp),parameter :: dos_max=9999.9999_dp
 real(dp) :: buffer,deltaene,enemax,enemin,integral_DOS,max_occ
 real(dp) :: cpu,wall,gflops
 logical :: bigDOS,iam_master
 character(len=10) :: tag
 character(len=500) :: frmt,frmt_extra,msg
 type(t_htetrahedron) :: tetra
!arrays
 integer,allocatable :: unt_atsph(:)
 real(dp) :: list_dp(3)
 real(dp),allocatable :: tmp_eigen(:),total_dos(:,:,:),eig_dos(:,:)
 real(dp),allocatable :: dos_m(:,:,:),dos_paw1(:,:,:),dos_pawt1(:,:,:)
 real(dp),allocatable :: wdt(:,:)

! *********************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm); iam_master = (my_rank == master)

 prtdosm = dos%prtdosm; paw_dos_flag = dos%paw_dos_flag
 mbesslang = dos%mbesslang; ndosfraction = dos%ndosfraction

 nkpt = dtset%nkpt; nsppol = dtset%nsppol

!m-decomposed DOS not compatible with PAW-decomposed DOS
 if (prtdosm>=1.and.paw_dos_flag==1) then
   msg = 'm-decomposed DOS (prtdosm>=1) not compatible with PAW-decomposed DOS (pawprtdos=1) !'
   MSG_ERROR(msg)
 end if

!Refuse nband different for different kpoints
!Note: This means we can pass ebands%eig(:,:,:) instead of eigen(mband*nkpt*nsppol) in packed form
 do isppol=1,nsppol
   do ikpt=1,nkpt
     if ( dtset%nband(nkpt*(isppol-1) + ikpt) /= dtset%nband(1) ) then
       write(std_out,*) 'tetrahedron: skip subroutine.'
       write(std_out,*) 'nband must be the same for all kpoints'
       write(std_out,*) 'nband=', dtset%nband
       MSG_WARNING('tetrahedron: skip subroutine. See message above')
       return
     end if
   end do
 end do

 call cwtime(cpu, wall, gflops, "start")

 tetra = tetra_from_kptrlatt(crystal, dtset%kptopt, dtset%kptrlatt, dtset%nshiftk, &
   dtset%shiftk, dtset%nkpt, dtset%kpt, comm, msg, ierr)
 if (ierr /= 0) then
   call tetra%free()
   MSG_WARNING(msg)
   return
 end if

 natsph=dtset%natsph; natsph_extra=dtset%natsph_extra

 ! Master opens the DOS files.
 if (iam_master) then
   if (any(dtset%prtdos == [2, 5])) then
     if (open_file(fildata, msg, newunit=unitdos, status='unknown', form='formatted', action="write") /= 0) then
       MSG_ERROR(msg)
     end if

   else if (dtset%prtdos == 3) then
     ! unt_atsph(0) is used for the total DOS.
     ABI_MALLOC(unt_atsph,(0:natsph+natsph_extra))

     ! Open file for total DOS as well.
     if (open_file(strcat(fildata, '_TOTAL'), msg, newunit=unt_atsph(0), &
         status='unknown', form='formatted', action="write") /= 0) then
       MSG_ERROR(msg)
     end if

     do iat=1,natsph
       call int2char4(dtset%iatsph(iat),tag)
       ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
       if (open_file(strcat(fildata, '_AT', tag), msg, newunit=unt_atsph(iat), &
           status='unknown', form='formatted', action="write") /= 0) then
         MSG_ERROR(msg)
       end if
     end do
     ! do extra spheres in vacuum too. Use _ATEXTRA[NUM] suffix
     do iat=1,natsph_extra
       call int2char4(iat,tag)
       ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
       if (open_file(strcat(fildata, '_ATEXTRA', tag), msg, newunit=unt_atsph(natsph+iat), &
           status='unknown', form='formatted', action="write") /= 0) then
         MSG_ERROR(msg)
       end if
     end do
   end if
 end if

 ! Write the header of the DOS file, and determine the energy range and spacing
 prtdos=dtset%prtdos
 buffer=0.01_dp ! Size of the buffer around the min and max ranges

 ! A Similar section is present is getnel. Should move all DOS stuff to m_ebands
 ! Choose the lower and upper energies
 enemax = maxval(ebands%eig) + buffer
 enemin = minval(ebands%eig) - buffer

 ! Extend the range to a nicer value
 enemax=0.1_dp*ceiling(enemax*10._dp)
 enemin=0.1_dp*floor(enemin*10._dp)

 ! Choose the energy increment
 if(abs(dtset%dosdeltae)<tol10)then
   deltaene=0.001_dp
   if(dtset%prtdos>=2)deltaene=0.0005_dp ! Higher resolution possible (and wanted) for tetrahedron
 else
   deltaene=dtset%dosdeltae
 end if
 nene=nint((enemax-enemin)/deltaene)+1

 call xmpi_bcast(nene, master, comm, ierr)
 if (iam_master) list_dp(1:3) = [deltaene, enemin, enemax]
 call xmpi_bcast(list_dp, master, comm, ierr)
 deltaene = list_dp(1); enemin = list_dp(2); enemax = list_dp(3)

 if (iam_master) then
   if (any(dtset%prtdos == [2, 5])) then
     call dos_hdr_write(deltaene,ebands%eig,enemax,enemin,ebands%fermie,dtset%mband,&
     dtset%nband,nene,nkpt,nsppol,dtset%occopt,prtdos,&
     dtset%tphysel,dtset%tsmear,unitdos)
   else if (dtset%prtdos == 3) then
     do iat=0,natsph+natsph_extra
       call dos_hdr_write(deltaene,ebands%eig,enemax,enemin,ebands%fermie,dtset%mband,&
       dtset%nband,nene,nkpt,nsppol,dtset%occopt,prtdos,&
       dtset%tphysel,dtset%tsmear,unt_atsph(iat))
     end do
   end if
 end if

 ! Tetra weights
 ABI_MALLOC(wdt, (nene, 2))

 ! Allocate arrays to store DOSes and fill with zeros.
 ! 1--> DOS , 2--> IDOS
 ABI_MALLOC(total_dos,(nene,ndosfraction,2))
 ABI_MALLOC(eig_dos, (nene, 2))

 if (paw_dos_flag==1) then
   ABI_MALLOC(dos_paw1,(nene,ndosfraction,2))
   ABI_MALLOC(dos_pawt1,(nene,ndosfraction,2))
 end if
 if (prtdosm>=1) then
   ABI_MALLOC(dos_m, (nene,ndosfraction*mbesslang,2))
 end if

!Get maximum occupation value (2 or 1)
 max_occ = one; if (dtset%nspinor == 1 .and. nsppol == 1) max_occ = two

!-------------------------------------------------------------------
!For each spin polarisation and band, interpolate band over kpoints
!calculate integration weights and DOS contib from
!-------------------------------------------------------------------

 ! Workspace arrays.
 ABI_MALLOC(tmp_eigen,(nkpt))

 cnt = 0
 do isppol=1,nsppol

   total_dos = zero; eig_dos = zero
   if (prtdosm>=1) dos_m = zero
   if (paw_dos_flag==1) then
     dos_paw1 = zero; dos_pawt1 = zero
   end if

   do ikpt=1,nkpt
      do iband=1,ebands%nband(ikpt+(isppol-1)*ebands%nkpt)
        cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! Mpi parallelism.

        ! Accumulate total DOS from eigenvalues (this is the **exact** total DOS)
        tmp_eigen(:) = ebands%eig(iband, :, isppol)
        call tetra%get_onewk(ikpt,bcorr0,nene,nkpt,tmp_eigen,enemin,enemax,max_occ,wdt)
        wdt = wdt*ebands%wtk(ikpt)
        eig_dos = eig_dos + wdt

        ! Accumulate L-DOS.
        do ii=1,2
          do ifrac=1,ndosfraction
            total_dos(:,ifrac,ii) = total_dos(:,ifrac,ii) + wdt(:,ii) * dos%fractions(ikpt,iband,isppol,ifrac)
          end do
        end do

        if (paw_dos_flag==1) then
          ! Accumulate L-DOS (on-site terms).
          do ii=1,2
            do ifrac=1,ndosfraction
              dos_paw1(:,ifrac,ii) = dos_paw1(:,ifrac,ii) + wdt(:,ii) * dos%fractions_paw1(ikpt,iband,isppol,ifrac)
              dos_pawt1(:,ifrac,ii) = dos_pawt1(:,ifrac,ii) + wdt(:,ii) * dos%fractions_pawt1(ikpt,iband,isppol,ifrac)
            end do
          end do
        end if

        if (prtdosm>=1) then
         ! Accumulate LM-DOS.
         do ii=1,2
           do ifrac=1,ndosfraction*mbesslang
             dos_m(:,ifrac,ii) = dos_m(:,ifrac,ii) + wdt(:,ii) * dos%fractions_m(ikpt, iband, isppol, ifrac)
           end do
         end do
        end if

      end do ! ikpt
   end do ! iband

   ! Collect results on master
   call xmpi_sum_master(eig_dos, master, comm, ierr)
   call xmpi_sum_master(total_dos, master, comm, ierr)
   bigDOS=(maxval(total_dos(:,:,1))>999._dp)

   if (paw_dos_flag == 1) then
     call xmpi_sum_master(dos_paw1, master, comm, ierr)
     call xmpi_sum_master(dos_pawt1, master, comm, ierr)
   end if
   if (prtdosm >= 1) call xmpi_sum_master(dos_m, master, comm, ierr)

   ! Write the DOS value in the DOS file
   ! Print the data for this energy. Note the upper limit (dos_max), to be consistent with the format.
   ! The use of "E" format is not adequate, for portability of the self-testing procedure.
   ! header lines depend on the type of DOS (projected etc...) which is output

   if (.not. iam_master) goto 10
   call write_extra_headers()

   if (prtdos==2) then
     ! E, DOS, IDOS
     do iene=1,nene
       write(unitdos, '(f11.5,1x,2(f10.4,1x))') &
        enemin + (iene-1)*deltaene, min(total_dos(iene,:,1), dos_max), total_dos(iene,:,2)
     end do

   else if (prtdos==3) then

     ! Write E, DOS, IDOS
     do iene=1,nene
       write(unt_atsph(0), '(f11.5,1x,2(f10.4,1x))') &
        enemin + (iene-1)*deltaene, min(eig_dos(iene,1), dos_max), eig_dos(iene,2)
     end do

     ! E, DOS(L=1,LMAX), IDOS(L=1,LMAX)
     ! Here we assume mpsang = 5 in the format.
     if (paw_dos_flag/=1.or.dtset%pawprtdos==2) then
       frmt = '(f11.5,1x,5(f9.4,1x),10x,5(f8.2,1x),10x,25(f8.2,1x))'
       if (bigDOS) frmt = '(f11.5,1x,5(f10.4,1x),10x,5(f8.2,1x),10x,25(f8.2,1x))'
       ! for extra atoms in vacuum need more precision
       frmt_extra = '(f11.5,1x,5(f20.16,1x),10x,5(f20.16,1x),10x,25(f20.16,1x))'

       do iat=1,natsph
         i1 = (iat-1)*mbesslang+1; i2 = iat*mbesslang
         if (prtdosm==0) then
           do iene=1,nene
             write(unt_atsph(iat), fmt=frmt) enemin + (iene-1)*deltaene, &
&              min(total_dos(iene, i1:i2, 1), dos_max), total_dos(iene, i1:i2,2)
           end do
         else
           do iene=1,nene
             write(unt_atsph(iat), fmt=frmt) enemin + (iene-1)*deltaene, &
&              min(total_dos(iene, i1:i2, 1), dos_max),&
&              total_dos(iene, i1:i2, 2),&
&              min(dos_m(iene,(iat-1)*mbesslang**2+1:iat*mbesslang**2,1), dos_max)
           end do
         end if
       end do

       ! Extra spheres.
       do iat=natsph+1,natsph+natsph_extra
         i1 = (iat-1)*mbesslang+1; i2 = iat*mbesslang
         if (prtdosm==0) then
           do iene=1,nene
             write(unt_atsph(iat), fmt=frmt_extra) enemin + (iene-1)*deltaene, &
&             total_dos(iene, i1:i2, 1), &
&             total_dos(iene, i1:i2, 2)
           end do
         else
           do iene=1,nene
             write(unt_atsph(iat), fmt=frmt_extra) enemin + (iene-1)*deltaene, &
&             total_dos(iene, i1:i2, 1),&
&             total_dos(iene, i1:i2, 2),&
&             dos_m(iene,(iat-1)*mbesslang**2+1:iat*mbesslang**2, 1)
           end do
         end if
       end do

     else
       frmt = '(f11.5,1x,5(f9.4,1x),3(6x,5f9.4))'
       if (bigDOS) frmt = '(f11.5,1x,5(f10.4,1x),3(6x,5f10.4))'
       ! for extra atom spheres in vacuum need more precision
       frmt_extra = '(f11.5,1x,5(f20.16,1x),3(6x,5f20.16))'

       do iat=1,natsph
         i1 = iat*5-4; i2 = iat*5
         do iene=1,nene
           write(unt_atsph(iat), fmt=frmt) enemin + (iene-1)*deltaene, &
&           min(total_dos(iene,i1:i2,1), dos_max),&
&           min(total_dos(iene,i1:i2,1) - dos_paw1(iene,i1:i2,1) + dos_pawt1(iene,i1:i2,1), dos_max),&
&           min(dos_paw1(iene,i1:i2,1), dos_max),&
&           min(dos_pawt1(iene,i1:i2,1), dos_max)
         end do
       end do

       ! Extra spheres.
       do iat=natsph+1,natsph+natsph_extra
         i1 = iat*5-4; i2 = iat*5
         do iene=1,nene
           write(unt_atsph(iat), fmt=frmt_extra) enemin + (iene-1)*deltaene, &
&            min(total_dos(iene,i1:i2,1), dos_max),&
&            min(total_dos(iene,i1:i2,1) - dos_paw1(iene,i1:i2,1) + dos_pawt1(iene,i1:i2,1), dos_max),&
&            min(dos_paw1(iene,i1:i2,1), dos_max),&
&            min(dos_pawt1(iene,i1:i2,1), dos_max)
         end do
       end do
     end if

   else if (prtdos==5)then
     ! E, SPIN-DOS
     frmt = '(f11.5,1x,7(f9.4,1x),10x,7(f8.2,1x))'
     if (bigDOS) frmt = '(f11.5,1x,7(f10.4,1x),10x,7(f8.2,1x))'
     do iene=1,nene
       write(unitdos, fmt=frmt) enemin + (iene-1)*deltaene, min(total_dos(iene,1:7,1), dos_max), total_dos(iene,1:7,2)
     end do
   end if

10 continue
   integral_DOS=sum(total_dos(nene,:,2))
   write(msg, '(a,es16.8)' ) ' tetrahedron : integrate to',integral_DOS
   call wrtout(std_out,msg,'COLL')
 end do ! isppol

 ! Close files.
 if (iam_master) then
   if (any(prtdos == [2, 5])) then
     close(unitdos)
   else if (prtdos == 3) then
     do iat=0,natsph+natsph_extra
       close(unt_atsph(iat))
     end do
     ABI_FREE(unt_atsph)
   end if
 end if

 ABI_FREE(tmp_eigen)
 ABI_FREE(total_dos)
 ABI_FREE(wdt)
 ABI_FREE(eig_dos)

 if (prtdosm>=1)  then
   ABI_FREE(dos_m)
 end if

 if (paw_dos_flag==1)  then
   ABI_FREE(dos_paw1)
   ABI_FREE(dos_pawt1)
 end if

 call tetra%free()

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2),a)')" tetrahedron: cpu_time: ",cpu,"[s], walltime: ",wall," [s]"
 call wrtout(std_out,msg,"PERS")

contains

subroutine write_extra_headers()

 if (nsppol==2) then
   if(isppol==1) write(msg,'(a,16x,a)')  '#','Spin-up DOS'
   if(isppol==2) write(msg,'(2a,16x,a)')  ch10,'#','Spin-dn DOS'
   ! NB: dtset%prtdos == 5 should not happen for nsppol==2

   if (any(dtset%prtdos == [2, 5])) then
     write(unitdos, "(a)")trim(msg)

   else if (dtset%prtdos == 3) then
     do iat=0,natsph+natsph_extra
       write(unt_atsph(iat), "(a)")trim(msg)
     end do
   end if

 end if

 if (prtdos==2) then
   write(unitdos, '(a)' )'# energy(Ha)     DOS  integrated DOS'

 else if (prtdos==3) then

   write(unt_atsph(0), '(a)' )'# energy(Ha)     DOS  integrated DOS'

   if (paw_dos_flag/=1.or.dtset%pawprtdos==2) then
     do iat=1,natsph
       write(unt_atsph(iat), '(3a,i5,a,i5,a,a,es16.6,3a)' ) &
&       '# Local DOS (columns 2-6) and integrated local DOS (columns 7-11),',ch10,&
&       '# for atom number iat=',iat,'  iatom=',dtset%iatsph(iat),ch10,&
&       '# inside sphere of radius ratsph=',dtset%ratsph(dtset%typat(dtset%iatsph(iat))),' Bohr.',ch10,"#"

       if (dtset%usepaw==1.and.dtset%pawprtdos==2) then
         write(unt_atsph(iat), '(3a)' ) &
&         '# PAW: note that only all-electron on-site part has been used to compute DOS !',ch10,"#"
       end if
       if (bigDOS) then
         write(msg, '(a,a)' ) &
&         '# energy(Ha)   l=0       l=1       l=2       l=3       l=4',&
&         '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
       else
         write(msg, '(a,a)' ) &
&         '# energy(Ha)  l=0      l=1      l=2      l=3      l=4',&
&         '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
       end if
       if (prtdosm>=1) then
         write(msg, '(7a)' ) trim(msg),'          ',&
&         '  lm=0 0',&
&         '  lm=1-1  lm=1 0  lm=1 1',&
&         '  lm=2-2  lm=2-1  lm=2 0  lm=2 1  lm=2 2',&
&         '  lm=3-3  lm=3-2  lm=3-1  lm=3 0  lm=3 1  lm=3 2  lm=3 3',&
&         '  lm=4-4  lm=4-3  lm=4-2  lm=4-1  lm=4 0  lm=4 1  lm=4 2  lm=4 3  lm=4 4'
       end if
       write(unt_atsph(iat), "(a)")trim(msg)
     end do
   else
     do iat=1,natsph
       write(unt_atsph(iat), '(9a,i5,a,i5,a,a,es16.6,3a)' ) &
&       '# Local DOS (columns 2-6),',ch10,&
&       '#  plane-waves contrib. to DOS (columns 7-11),',ch10,&
&       '#  AE on-site  contrib. to DOS (columns 12-16),',ch10,&
&       '# -PS on-site  contrib. to DOS (columns 17-21),',ch10,&
&       '# for atom number iat=',iat,'  iatom=',dtset%iatsph(iat),ch10,&
&       '# inside sphere of radius ratsph=',dtset%ratsph(dtset%typat(dtset%iatsph(iat))),' Bohr.',ch10,"#"
       if (bigDOS) then
         write(msg, '(4a)' ) &
&         '#energy(Ha)   l=0       l=1       l=2       l=3       l=4',&
&         '       (PW)  l=0       l=1       l=2       l=3       l=4',&
&         '      (Phi)  l=0       l=1       l=2       l=3       l=4',&
&         '     (tPhi)  l=0       l=1       l=2       l=3       l=4'
       else
         write(msg, '(4a)' ) &
&         '#energy(Ha)  l=0      l=1      l=2      l=3      l=4',&
&         '       (PW) l=0      l=1      l=2      l=3      l=4',&
&         '      (Phi) l=0      l=1      l=2      l=3      l=4',&
&         '     (tPhi) l=0      l=1      l=2      l=3      l=4'
       end if
       write(unt_atsph(iat), "(a)")trim(msg)
     end do
   end if
   do iat=1,natsph_extra
     write(unt_atsph(natsph+iat), '(3a,i5,2a,es16.6,3a)' ) &
&     '# Local DOS (columns 2-6) and integrated local DOS (columns 7-11),',ch10,&
&     '# for non-atomic sphere number iat=',iat,ch10,&
&     '# of radius ratsph=',dtset%ratsph_extra,' Bohr.',ch10,"#"
     if (bigDOS) then
       write(msg, '(a,a)' ) &
&       '# energy(Ha)   l=0       l=1       l=2       l=3       l=4',&
&       '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
     else
       write(msg, '(a,a)' ) &
&       '# energy(Ha)  l=0      l=1      l=2      l=3      l=4',&
&       '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
     end if
     if (prtdosm>=1) then
       write(msg, '(7a)' ) trim(msg),'          ',&
&       '  lm=0 0',&
&       '  lm=1-1  lm=1 0  lm=1 1',&
&       '  lm=2-2  lm=2-1  lm=2 0  lm=2 1  lm=2 2',&
&       '  lm=3-3  lm=3-2  lm=3-1  lm=3 0  lm=3 1  lm=3 2  lm=3 3',&
&       '  lm=4-4  lm=4-3  lm=4-2  lm=4-1  lm=4 0  lm=4 1  lm=4 2  lm=4 3  lm=4 4'
     end if
     write(unt_atsph(natsph+iat), "(a)")trim(msg)
   end do

 else if (prtdos==5) then
   write(unitdos, '(a)' )&
&    '# energy(Ha)     DOS up,up  up,dn  dn,up  dn,dn  sigma_x sigma_y sigma_z  and integrated DOS components'
 end if ! prtdos value

end subroutine write_extra_headers

end subroutine dos_calcnwrite
!!***

!!****f* m_epjdos/recip_ylm
!! NAME
!! recip_ylm
!!
!! FUNCTION
!! Project input wavefunctions in reciprocal space on to Ylm
!! (real or complex harmonics depending on rc_ylm).
!!
!! INPUTS
!!  bess_fit(mpw,nradintmax,ll) = Bessel functions for L, splined
!!   with arguments $2 \pi |k+G| \Delta r$, for all G vectors in sphere
!!   and all points on radial grid.
!!  cg_1band(2,npw_k)=wavefunction in recip space (note that nspinor is missing, see Notes).
!!  istwfk= storage mode of cg_1band
!!  nradint(natsph)=number of points on radial real-space grid for a given atom.
!!  nradintmax=dimension of rint array.
!!  me_g0=1 if this processor has G=0, 0 otherwise
!!  mlang=maximum angular momentum in Bessel functions.
!!  mpw=Maximum number of planewaves. Used to dimension bess_fit
!!  natsph=number of atoms around which ang mom projection has to be done
!!  typat_extra(natsph)=Type of each atom. ntypat + 1 if empty sphere
!!  mlang_type(ntypat + natsph_extra)=Max L+1 for each atom type
!!  npw_k=number of plane waves for this kpt
!!  nspinor=number of spinor components
!!  ph3d(2,npw_k,natsph)=3-dim structure factors, for each atom and plane wave.
!!  prtsphere= if 1, print a complete analysis of the angular momenta in atomic spheres
!!  rint(nradintmax) = points on radial real-space grid for integration
!!  rmax(natsph)=maximum radius for real space integration sphere
!!  rc_ylm= 1 for real spherical harmonics. 2 for complex spherical harmonics,
!!  ucvol=unit cell volume in bohr**3.
!!  ylm_k(npw_k,mlang**2)=real spherical harmonics for each G and LM.
!!  znucl_sph(natsph)=gives the nuclear number for each type of atom
!!
!! OUTPUT
!!  sum_1ll_1atom(mlang,natsph)= projected scalars for each atom and ang. mom.
!!  sum_1lm_1atom(mlang*mlang,natsph)= projected scalars for each atom and LM component.
!!  cplx_1lm_1atom(2,dtset%nspinor**2,dos%mbesslang**2,natsph_tot) = complex projection of wave function on atomic like orbital
!!
!! NOTES
!!  * ph3d atoms are ordered with natsph and must be provided by the caller in the correct order!
!!
!!  * spinor components are not treated here. This facilitates the implementation of spinor parallelism
!!    because the caller can easily call the routine inside a loop over spinors and then sum the
!!    different contributions outside the loop thus reducing the number of MPI calls.
!!
!! PARENTS
!!      m_cut3d,partial_dos_fractions
!!
!! CHILDREN
!!      cwtime,wrtout
!!
!! SOURCE

subroutine recip_ylm (bess_fit, cg_1band, istwfk, mpi_enreg, nradint, nradintmax, mlang,&
&  mpw, natsph, typat_extra, mlang_type, npw_k, nspinor, ph3d, prtsphere, rint, rmax,&
&  rc_ylm, sum_1ll_1atom, sum_1lm_1atom, cplx_1lm_1atom, ucvol, ylm_k, znucl_sph)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwfk,mlang,mpw,natsph,npw_k,nradintmax
 integer,intent(in) :: nspinor
 integer,intent(in) :: prtsphere,rc_ylm
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: nradint(natsph),typat_extra(natsph),mlang_type(:)
 real(dp),intent(in) :: bess_fit(mpw,nradintmax,mlang),cg_1band(:,:) !(2,my_nspinor*npw_k)
 real(dp),intent(in) :: ph3d(2,npw_k,natsph),rint(nradintmax)
 real(dp),intent(in) :: rmax(natsph),ylm_k(npw_k,mlang*mlang)
 real(dp),intent(in) :: znucl_sph(natsph)
 type(MPI_type),intent(in) :: mpi_enreg
 real(dp),intent(out) :: sum_1ll_1atom(nspinor**2, mlang, natsph)
 real(dp),intent(out) :: sum_1lm_1atom(nspinor**2, mlang*mlang,natsph)
 real(dp),intent(out) :: cplx_1lm_1atom(2,nspinor, mlang*mlang,natsph)

!Local variables-------------------------------
!scalars
 integer :: ilm,iat,ipw,ixint,ll,mm,il,jlm,ierr,lm_size,itypat
 integer :: ispinor
 integer :: ipauli, is, isp
 integer :: my_nspinor
 real(dp),parameter :: invsqrt2=one/sqrt2
 real(dp) :: sum_all, dr, fact
 type(atomdata_t) :: atom
 character(len=500) :: msg
!arrays
 integer :: ilang(mlang**2)
 real(dp) :: c1(2),c2(2)
 real(dp) :: sum_1atom(natsph),sum_1ll(mlang),sum_1lm(mlang**2)
 real(dp) :: func(nradintmax)
 real(dp) :: func_cplx(nradintmax,2)
 complex(dpc) :: vect(npw_k)
 complex(dpc),allocatable :: tmppsia(:,:),tmppsim(:,:),dotc(:)
 integer, allocatable :: ispinors(:)
 complex(dpc),allocatable :: values(:,:,:,:)

! *************************************************************************

 ! Workspace array (used to reduce the number of MPI communications)
 ! One could reduce a bit the memory requirement by using non-blocking operations ...
 ABI_MALLOC_OR_DIE(values, (nradintmax, nspinor, mlang**2, natsph), ierr)
 values = czero

 my_nspinor = max(1,nspinor/mpi_enreg%nproc_spinor)
 ABI_ALLOCATE(tmppsia, (npw_k,my_nspinor))
 ABI_ALLOCATE(tmppsim, (npw_k,my_nspinor))
 ABI_ALLOCATE(dotc, (my_nspinor))
 ABI_ALLOCATE(ispinors, (my_nspinor))
 if (my_nspinor == 2) then
   ispinors(1) = 1
   ispinors(2) = 2
 else
   ispinors(1) = mpi_enreg%me_spinor+1
 end if

 sum_1lm_1atom = zero
 cplx_1lm_1atom = zero

 do ll=0,mlang-1
   do mm=-ll,ll
     ilm = (ll+1)**2-ll+mm
     ilang(ilm) = ll+1
   end do
 end do

 ! Big loop on all atoms
 do iat=1,natsph
   itypat = typat_extra(iat)
   lm_size = mlang_type(itypat) ** 2
   dr = rmax(iat) / (nradint(iat)-1)

   ! u(G) e^{i(k+G).Ra}
   ! tmppsia = Temporary array for part which depends only on iat
   do ispinor=1,my_nspinor
     do ipw=1,npw_k
       tmppsia(ipw,ispinor) = dcmplx(cg_1band(1,ipw+(ispinor-1)*npw_k),cg_1band(2,ipw+(ispinor-1)*npw_k)) &
&            * dcmplx(ph3d(1,ipw,iat), ph3d(2,ipw,iat))
     end do
   end do

   ! tmppsim = temporary arrays for part of psi which does not depend on ixint = tmppsia * ylm.
   ! could remove this intermediate array to save memory...
   ! u(G) Y_LM^*(k+G) e^{i(k+G).Ra}
   ! Take into account the fact that ylm_k are REAL spherical harmonics, see initylmg.f
   ! For time-reversal states, detailed treatment show that only the real or imaginary
   ! part of tmppsia is needed here, depending on l being even or odd: only one of the coef is 1, the other 0
   do ilm=1,lm_size
     il = ilang(ilm)
     ll = ilang(ilm) - 1
     mm = ilm - (ll+1)**2 + ll

     select case (rc_ylm)
     case (1)
       ! to get PDOS for real spherical harmonics, simply multiply here by ylm_k
       do ispinor=1,my_nspinor
         do ipw=1,npw_k
           tmppsim(ipw,ispinor) = tmppsia(ipw,ispinor) * ylm_k(ipw,ilm)
         end do
       end do

       ! Handle time-reversal
       ! TODO: check if time reversal with spinors is special and may need some treatment.
       !  normally SOC will simply impose istwfk 1 in appropriate cases (kptopt 4).
       if (istwfk /= 1) then
         if (mod(ll, 2) == 0) then
            tmppsim(:,:) = dcmplx(real(tmppsim(:,:)),zero)
         else
            tmppsim(:,:) = dcmplx(aimag(tmppsim(:,:)),zero)
         end if
! f2008 version:
!         if (mod(ll, 2) == 0) then
!            tmppsim(:,:)%im = zero
!         else
!            tmppsim(:,:)%re = tmppsim(:,:)%im
!            tmppsim(:,:)%im = zero
!         end if
       end if

     case (2)
       ! to get PDOS for complex spherical harmonics, build linear combination of real ylm_k
       jlm = (ll+1)**2-ll-mm ! index of (l, -m)
       if (mm == 0) then
         vect(:) = dcmplx(ylm_k(1:npw_k,ilm),zero)
       else if (mm > 0) then
          !vect(1,:) =  invsqrt2 * ylm_k(1:npw_k,ilm) * (-1)**mm
          !vect(2,:) = +invsqrt2 * ylm_k(1:npw_k,jlm) * (-1)**mm
          c1 = sy(ll, mm, mm)
          c2 = sy(ll,-mm, mm)
          vect(:) = dcmplx(c1(1) * ylm_k(1:npw_k,ilm) + c2(1) * ylm_k(1:npw_k,jlm), &
&                         c1(2) * ylm_k(1:npw_k,ilm) + c2(2) * ylm_k(1:npw_k,jlm))

       else if (mm < 0) then
          !vect(1,:) =  invsqrt2 * ylm_k(1:npw_k,jlm) !* (-1)**mm
          !vect(2,:) = -invsqrt2 * ylm_k(1:npw_k,ilm) !* (-1)**mm
          c1 = sy(ll, mm,  mm)
          c2 = sy(ll,-mm,  mm)
          vect(:) = dcmplx(c1(1) * ylm_k(1:npw_k,ilm) + c2(1) * ylm_k(1:npw_k,jlm),&
&                         c1(2) * ylm_k(1:npw_k,ilm) + c2(2) * ylm_k(1:npw_k,jlm))
       end if
       vect(:) = dcmplx(real(vect(:)), -aimag(vect(:)))
       !vect(:)%im = -vect(:)%im

       if (istwfk == 1) then
         do ispinor=1,my_nspinor
           do ipw=1,npw_k
             tmppsim(ipw,ispinor) = tmppsia(ipw,ispinor) * vect(ipw)
           end do
         end do
       else
         ! Handle time-reversal
         if (mod(ll, 2) == 0) then
           do ispinor=1,my_nspinor
             do ipw=1,npw_k
               tmppsim(ipw,ispinor) = real(tmppsia(ipw,ispinor)) * vect(ipw)
             end do
           end do
         else
           do ispinor=1,my_nspinor
             do ipw=1,npw_k
               tmppsim(ipw,ispinor) = aimag(tmppsia(ipw,ispinor)) * vect(ipw)
             end do
           end do
         end if
       end if

     case default
       MSG_ERROR("Wrong value for rc_ylm")
     end select

     ! Compute integral $ \int_0^{rc} dr r**2 ||\sum_G u(G) Y_LM^*(k+G) e^{i(k+G).Ra} j_L(|k+G| r)||**2 $
     ! or more general spinor case integral
     !    $ \int_0^{rc} dr r**2 dotc^*_s   \sigma^x_{ss'}   dotc_{s'}
     ! where   dotc_s = \sum_G u_s (G) Y_LM^*(k+G) e^{i(k+G).Ra} j_L(|k+G| r)
     do ixint=1,nradint(iat)
       dotc = czero
       do ispinor=1, my_nspinor
         do ipw=1,npw_k
           dotc(ispinor) = dotc(ispinor) + bess_fit(ipw, ixint, il) * tmppsim(ipw, ispinor)
         end do
       end do
       if (istwfk /= 1) then
         dotc = two * dotc
         if (istwfk == 2 .and. mpi_enreg%me_g0 == 1) then
           dotc(:) = dotc(:) - bess_fit(1, ixint, il) * tmppsim(1, :)
         end if
       end if

       ! Store results to reduce number of xmpi_sum calls if MPI
       do ispinor=1, my_nspinor
         values(ixint, ispinors(ispinor), ilm, iat) = dotc(ispinor)
       end do
     end do ! ixint

   end do ! ilm
 end do ! iat

 ABI_DEALLOCATE(tmppsia)
 ABI_DEALLOCATE(tmppsim)
 ABI_DEALLOCATE(dotc)
 ABI_DEALLOCATE(ispinors)

 ! Collect results in comm_pw (data are distributed over plane waves)
 call xmpi_sum(values, mpi_enreg%comm_bandfft, ierr)
! ! Collect results in mpi_enreg%comm_spinor (data are distributed over spinor components)
 call xmpi_sum(values, mpi_enreg%comm_spinor, ierr)

 ! Multiply by r**2 and take norm, integrate
 do iat=1,natsph
   itypat = typat_extra(iat)
   lm_size = mlang_type(itypat) ** 2
   do ilm=1,lm_size

     do ipauli=0,nspinor**2-1
       do ixint=1,nradint(iat)
         func(ixint) = zero
         do is=1,nspinor
           do isp=1,nspinor
             func(ixint) =  func(ixint) + real(conjg(values(ixint, is, ilm, iat)) * pauli_mat(is,isp,ipauli)*&
&                                                    values(ixint, isp, ilm, iat))
           end do
         end do
         func(ixint) = rint(ixint)**2 * func(ixint)
       end do
       ! Here I should treat the case in which the last point /= rcut
       ! NB: indexing is from 1 not 0 for spin matrix components
       sum_1lm_1atom   (ipauli+1, ilm, iat) = simpson(dr, func(1:nradint(iat)))
     end do ! ipauli

     do is = 1, nspinor
       func_cplx(:,1) = real(values(:, is, ilm, iat))
       func_cplx(:,2) = aimag(values(:, is, ilm, iat))
       do ixint=1,nradint(iat)
         func_cplx(ixint,:) = rint(ixint)**2 * func_cplx(ixint,:)
       end do
       cplx_1lm_1atom(1, is, ilm, iat) = simpson(dr, func_cplx(1:nradint(iat),1))
       cplx_1lm_1atom(2, is, ilm, iat) = simpson(dr, func_cplx(1:nradint(iat),2))
     end do

   end do ! ilm
 end do ! iat

 ! Normalize with unit cell volume and include 4pi term coming from Rayleigh expansion.
 fact = four_pi**2 / ucvol
 sum_1lm_1atom = fact * sum_1lm_1atom
 cplx_1lm_1atom = fact * cplx_1lm_1atom

 ! sum up the m-independent fractions
 sum_1ll_1atom = zero
 do iat=1,natsph
   itypat = typat_extra(iat)
   lm_size = mlang_type(itypat) ** 2
   do ilm=1,lm_size
     il = ilang(ilm)
     sum_1ll_1atom(:,il, iat) = sum_1ll_1atom(:,il, iat) + sum_1lm_1atom(:,ilm, iat)
   end do
 end do

 ABI_FREE(values)

 ! Output
 if (prtsphere == 1) then
   sum_1ll = zero
   sum_1lm = zero
   sum_1atom = zero
   do iat=1,natsph
     sum_1atom(iat) = sum(sum_1lm_1atom(1,:,iat))
     sum_1ll(:)=sum_1ll(:)+sum_1ll_1atom(1,:,iat)
     sum_1lm(:)=sum_1lm(:)+sum_1lm_1atom(1,:,iat)
   end do
   sum_all = sum(sum_1atom)

   if (rc_ylm == 1) msg = " Angular analysis (real spherical harmonics)"
   if (rc_ylm == 2) msg = " Angular analysis (complex spherical harmonics)"
   call wrtout(std_out, msg)
   do iat=1,natsph
     call atomdata_from_znucl(atom, znucl_sph(iat))
     call wrtout(std_out, " ")
     write(msg,'(a,i3,a,a,a,f10.6)' )' Atom # ',iat, ' is  ',  atom%symbol,', in-sphere charge =',sum_1atom(iat)
     call wrtout(std_out, msg)
     do ll=0,mlang-1
       write(msg,'(a,i1,a,f9.6,a,9f6.3)' )&
&       ' l=',ll,', charge=',sum_1ll_1atom(1,ll+1,iat),&
&       ', m=-l,l splitting:',sum_1lm_1atom(1,1+ll**2:(ll+1)**2,iat)
       call wrtout(std_out, msg)
     end do ! ll
   end do ! iat
   write(msg,'(a,a)') ch10,' Sum of angular contributions for all atomic spheres '
   call wrtout(std_out, msg)
   do ll=0,mlang-1
     write(msg,'(a,i1,a,f9.6,a,f9.6)' )&
&     ' l=',ll,', charge =',sum_1ll(ll+1),' proportion =',sum_1ll(ll+1)/sum_all
     call wrtout(std_out, msg)
   end do
   write(msg,'(a,a,f10.6)' ) ch10,' Total over all atoms and l=0 to 4 :',sum_all
   call wrtout(std_out, msg)
   call wrtout(std_out, " ")
 end if

contains

 function sy(ll, mm, mp)
   use  m_paw_sphharm, only : ys
   ! Computes the matrix element <Slm|Ylm'>
   integer,intent(in) :: ll,mm, mp

   real(dp) :: sy(2)
   complex(dpc) :: ys_val

   ! Computes the matrix element <Yl'm'|Slm>
   call ys(ll,mp,ll,mm,ys_val)
   !call ys(ll,mm,ll,mp,ys_val)
   sy(1) = real(ys_val)
   sy(2) = -aimag(ys_val)

 end function sy

end subroutine recip_ylm
!!***

!!****f* m_epjdos/dens_in_sph
!! NAME
!! dens_in_sph
!!
!! FUNCTION
!!  Calculate integrated density in sphere around each atom
!!
!! INPUTS
!!  cg      = wavefunction coefficitents in recip space
!!  gmet    = metric in recip space
!!  istwfk  = storage mode for cg coefficients
!!  kg_k    = G vector indices
!!  natom   = number of atoms
!!  mpi_enreg=information about MPI parallelization
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  npw_k   = number of plane waves for this kpoint
!!  ph1d    = phase factors for different atoms for all G vectors
!!  rmax(natom) = max radius to integrate to (in bohr)
!!
!! OUTPUT
!!  cmax = integrated density for each atom for a rmax-radius sphere
!!
!! WARNING
!!  cg should not be modified by fourwf.
!!
!! PARENTS
!!      m_cut3d
!!
!! CHILDREN
!!      cwtime,wrtout
!!
!! SOURCE

subroutine dens_in_sph(cmax,cg,gmet,istwfk,kg_k,natom,ngfft,mpi_enreg,npw_k,&
&                       ph1d,rmax,ucvol)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwfk,natom,npw_k
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: kg_k(3,npw_k),ngfft(18)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: ph1d(2,(2*ngfft(1)+1+2*ngfft(2)+1+2*ngfft(3)+1)*natom)
 real(dp),intent(in) :: rmax(natom)
 real(dp),intent(inout) :: cg(2,npw_k)
 real(dp),intent(out) :: cmax(natom)

!Local variables -------------------------
!scalars
 integer,parameter :: tim_fourwf=0
 integer :: cplex,i1,i2,i3,iatom,id1,id2,id3,ifft,mgfft,n1,n2,n3,n4,n5,n6,nfft,nfftot
 real(dp) :: cmaxr,g1,g2,g3,norm,weight
!arrays
 integer :: ngfft_here(18)
 integer,allocatable :: garr(:,:),gbound(:,:)
 real(dp),allocatable :: denpot(:,:,:),fofgout(:,:),fofr(:,:,:,:),gnorm(:)
 real(dp),allocatable :: ph3d(:,:,:),phkxred(:,:),rhog(:,:),rhor(:)
 real(dp),allocatable :: sphrhog(:,:)

! *********************************************************************

 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)
 n4=ngfft(4)
 n5=ngfft(5)
 n6=ngfft(6)
 nfftot = n1*n2*n3
 nfft=n1*n2*n3
 ngfft_here(:) = ngfft(:)
!fourwf doesnt work with other options for mode 0 (fft G -> r)
 ngfft_here(7)=111
 ngfft_here(8)=256
 mgfft=maxval(ngfft_here(1:3))

 call sqnorm_g(norm,istwfk,npw_k,cg,mpi_enreg%me_g0,mpi_enreg%comm_fft)

 if (abs(one-norm) > tol6) then
   write(std_out,'(a,f8.5)' ) ' dens_in_sph : this state is not normalized : norm=',norm
 end if

!-----------------------------------------------------------------
!inverse FFT of wavefunction to real space => density in real space
!-----------------------------------------------------------------
 ABI_ALLOCATE(gbound,(2*mgfft+8,2))
 call sphereboundary(gbound,istwfk,kg_k,mgfft,npw_k)

 weight = one
 cplex=1
 ABI_ALLOCATE(denpot,(cplex*n4,n5,n6))
 denpot(:,:,:)=zero
 ABI_ALLOCATE(fofgout,(2,npw_k))
 ABI_ALLOCATE(fofr,(2,n4,n5,n6))
 call fourwf(cplex,denpot,cg,fofgout,fofr,gbound,gbound, &
& istwfk,kg_k,kg_k,mgfft,mpi_enreg,1,ngfft_here,npw_k,&
& npw_k,n4,n5,n6,1,tim_fourwf,weight,weight)
 ABI_DEALLOCATE(fofgout)
 ABI_DEALLOCATE(fofr)
 ABI_DEALLOCATE(gbound)

 norm = sum(denpot(:,:,:))/nfftot
 if (abs(one-norm) > tol6) then
   write(std_out,'(a,f8.5)') ' dens_in_sph : this state is not normalized in real space : norm=',norm
 end if

!-----------------------------------------------------------------
!FFT of new density: we obtain n(G) in rhog(1,:)
!-----------------------------------------------------------------

!Change the packing of the reciprocal space density
 ABI_ALLOCATE(rhor,(nfft))
 call fftpac(1,mpi_enreg,1,n1,n2,n3,n4,n5,n6,ngfft,rhor,denpot,1)

 ABI_ALLOCATE(rhog,(2,nfft))
 call fourdp(1,rhog,rhor,-1,mpi_enreg,nfft,1,ngfft,0)

 ABI_DEALLOCATE(rhor)
 ABI_DEALLOCATE(denpot)

 do ifft=1,nfft
   rhog(:,ifft) = rhog(:,ifft) / ucvol
 end do

!-----------------------------------------------------------------
!calculate norms of G vectors
!-----------------------------------------------------------------

 ABI_ALLOCATE(garr,(3,nfft))
 ABI_ALLOCATE(gnorm,(nfft))
 id3=ngfft(3)/2+2 ; id2=ngfft(2)/2+2 ; id1=ngfft(1)/2+2
 do i3=1,n3
   g3=i3-(i3/(id3))*ngfft(3)-one
   do i2=1,n2
     g2=i2-(i2/(id2))*ngfft(2)-one
     do i1=1,n1
       g1=i1-(i1/(id1))*ngfft(1)-one
       ifft=i1+(i2-1)*n1+(i3-1)*n1*n2
       garr(1,ifft)=nint(g1)
       garr(2,ifft)=nint(g2)
       garr(3,ifft)=nint(g3)
       gnorm(ifft)=sqrt(gmet(1,1)*g1*g1 + &
&       two*gmet(2,1)*g2*g1 + &
&       two*gmet(3,1)*g3*g1 + &
&       gmet(2,2)*g2*g2 + &
&       gmet(3,2)*g3*g2 + &
&       gmet(3,3)*g3*g3)
     end do
   end do
 end do

!-----------------------------------------------------------------
!For each atom call sphericaldens to calculate
!n(G) * 1/|G|^3  *  int_0^2*\pi*r_{max}*|G| 4 \pi y^2 j_0 (y) dy
!for all G vectors put into array sphrhog
!scalar product of phase factors with spherically convoluted density
!-----------------------------------------------------------------

!largest mem occupation = nfft * (2(sphrog) +2*1(ph3d) +3(garr) +2(rhog) +1(gnorm)) = nfft * 10
 ABI_ALLOCATE(sphrhog,(2,nfft))
 ABI_ALLOCATE(phkxred,(2,natom))
 phkxred(1,:)=one
 phkxred(2,:)=zero
 ABI_ALLOCATE(ph3d,(2,nfft,1))

 do iatom=1,natom

   call sphericaldens(rhog,gnorm,nfft,rmax(iatom),sphrhog)
!  -----------------------------------------------------------------
!  Compute the phases for the whole set of fft vectors
!  -----------------------------------------------------------------

   call ph1d3d(iatom,iatom,garr,natom,natom,nfft,ngfft(1),ngfft(2),ngfft(3),&
&   phkxred,ph1d,ph3d)

!  For the phase factors, take the compex conjugate, before evaluating the scalar product
   do ifft=1,nfft
     ph3d(2,ifft,1)=-ph3d(2,ifft,1)
   end do
   cplex=2
   call dotprod_v(cplex,cmaxr,nfft,1,0,ph3d,sphrhog,mpi_enreg%comm_fft)
   cmax(iatom) = cmaxr
!  write(std_out,'(a,i4,a,es14.6,a,es12.6)' )' dens_in_sph : At ', iatom, ' has ',cmaxr, ' el.s in a sphere of rad ', rmax
 end do

 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(gnorm)
 ABI_DEALLOCATE(garr)
 ABI_DEALLOCATE(sphrhog)
 ABI_DEALLOCATE(ph3d)
 ABI_DEALLOCATE(phkxred)

end subroutine dens_in_sph
!!***

!!****f* m_epjdos/sphericaldens
!! NAME
!! sphericaldens
!!
!! FUNCTION
!! Compute the convolution of a function with
!! the unity constant function over a sphere of radius rmax .
!! The function is to be given in reciprocal space,
!! the resulting function is also given in reciprocal space.
!! The routine needs the norm of the reciprocal space vectors.
!!
!! The resulting function in reciprocal space can give the
!! integral of the density in any sphere of that radius, centered
!! on any point, by a simple scalar product.
!!
!! INPUTS
!!  fofg(2,nfft)=initial function, in reciprocal space
!!  gnorm(nfft)=norm of the reciprocal space vectors
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  rmax=radius of the sphere
!!
!! OUTPUT
!!  sphfofg(2,nfft)=convoluted function, in reciprocal space
!!
!! PARENTS
!!      m_epjdos
!!
!! CHILDREN
!!      cwtime,wrtout
!!
!! SOURCE

subroutine sphericaldens(fofg,gnorm,nfft,rmax,sphfofg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft
 real(dp),intent(in) :: rmax
!arrays
 real(dp),intent(in) :: fofg(2,nfft),gnorm(nfft)
 real(dp),intent(out) :: sphfofg(2,nfft)

!Local variables-------------------------------
!scalars
 integer :: ifft
 real(dp) :: factor,int0yy,rmax_2pi,yy

! *************************************************************************

 rmax_2pi=two_pi*rmax
 factor=four_pi/(two_pi)**3

 do ifft=1,nfft
   if(abs(gnorm(ifft)) < tol12)then
     sphfofg(1,ifft)=fofg(1,ifft)*four_pi*third*rmax**3
     sphfofg(2,ifft)=fofg(2,ifft)*four_pi*third*rmax**3
   else
     yy=gnorm(ifft)*rmax_2pi
     int0yy=factor*(sin(yy)-yy*cos(yy))/(gnorm(ifft)**3)
     sphfofg(1,ifft)=fofg(1,ifft)*int0yy
     sphfofg(2,ifft)=fofg(2,ifft)*int0yy
   end if
 end do

end subroutine sphericaldens
!!***

!!****f* m_epjdos/prtfatbands
!! NAME
!! prtfatbands
!!
!! FUNCTION
!! Print dos_fractions_m in order to plot easily fatbands
!! if pawfatbnd=1  1 : fatbands are resolved in L.
!! if pawfatbnd=1  2 : fatbands are resolved in L and M.
!!
!! INPUTS
!!  dos_fractions_m(nkpt,mband,nsppol,ndosfraction*mbesslang*m_dos_flag)
!!               = m-resolved projected dos inside PAW sphere.
!!  dtset        = Input variables
!!  ebands<ebands_t>=Band structure data.
!!  pawfatbnd    = keyword for fatbands
!!  mbesslang    =maximum angular momentum for Bessel function expansion
!!  m_dos_flag   =option for the m-contributions to the partial DOS
!!  ndosfraction =natsph*mbesslang
!!
!! OUTPUT
!! (only writing)
!!
!! NOTES
!!  This routine should be called by master only
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      cwtime,wrtout
!!
!! SOURCE

subroutine prtfatbands(dos,dtset,ebands,fildata,pawfatbnd,pawtab)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: pawfatbnd
 type(epjdos_t),intent(in) :: dos
 type(ebands_t),intent(in) :: ebands
 type(dataset_type),intent(in) :: dtset
 character(len=fnlen),intent(in) :: fildata
!arrays
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer :: iall,il,iat,natsph,inbfatbands,iband,mband,ixfat,isppol,nkpt,lmax,ll,mm
 integer :: ikpt,nband_k,ndosfraction,mbesslang
 real(dp) :: xfatband,cpu,wall,gflops
 character(len=1) :: tag_l,tag_1m,tag_is
 character(len=2) :: tag_2m
 character(len=10) :: tag_il,tag_at,tag_grace
 character(len=1500) :: message
 character(len=fnlen) :: tmpfil
 type(atomdata_t) :: atom
!arrays
 integer,allocatable :: unitfatbands_arr(:,:)
 real(dp),allocatable :: eigenvalues(:,:,:)
 character(len=2) :: symbol

!*************************************************************************

 DBG_ENTER("COLL")

 ndosfraction = dos%ndosfraction; mbesslang = dos%mbesslang

 if(dos%prtdosm.ne.0) then
   write(message,'(3a)')&
&   'm decomposed dos is activated',ch10, &
&   'Action: deactivate it with prtdosm=0 !'
   MSG_ERROR(message)
 end if

 if(dtset%nspinor==2) then
   MSG_WARNING("Fatbands are not yet available in the case nspinor==2!")
 end if

 ABI_CHECK(allocated(dos%fractions_m), "dos%fractions_m is not allocated!")

 natsph=dtset%natsph
 nkpt=dtset%nkpt
 mband=dtset%mband

 if(natsph>1000) then
   write(message,'(3a)')&
&   'Too big number of fat bands!',ch10, &
&   'Action: decrease natsph in input file !'
   MSG_ERROR(message)
 end if

!--------------  PRINTING IN LOG
 call cwtime(cpu, wall, gflops, "start")
 write(message,'(a,a,a,a,i5,a,a,1000i5)') ch10," ***** Print of fatbands activated ****** ",ch10,&
& "  Number of atom: natsph = ",natsph,ch10, &
& "  atoms  are             = ",(dtset%iatsph(iat),iat=1,natsph)
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 iall=0;inbfatbands=0

 if(pawfatbnd==1) then
   inbfatbands=mbesslang-1
   write(message,'(3a)')"  (fatbands are in eV and are given for each value of L)",ch10
 else if(pawfatbnd==2) then
   write(message,'(3a)')"  (fatbands are in eV and are given for each value of L and M)",ch10
   inbfatbands=(mbesslang-1)**2
 end if
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 write(message,'(a,e12.5,a,e12.5,a)') "  Fermi energy is ",ebands%fermie*Ha_eV," eV = ",ebands%fermie," Ha"
 call wrtout(std_out,message,'COLL')

!--------------  OPEN AND NAME FILES FOR FATBANDS
 ABI_ALLOCATE(unitfatbands_arr,(natsph*inbfatbands,dtset%nsppol))
 unitfatbands_arr = -3

 do iat=1,natsph
   lmax=(pawtab(dtset%typat(dtset%iatsph(iat)))%l_size-1)/2
   call int2char4(dtset%iatsph(iat),tag_at)
   ABI_CHECK((tag_at(1:1)/='#'),'Bug: string length too short!')
   call atomdata_from_znucl(atom,dtset%znucl(dtset%typat(dtset%iatsph(iat))))
   symbol = atom%symbol
   do il=1,inbfatbands
     iall=iall+1
     ll=int(sqrt(float(il-1)))  ! compute l
     if(ll.le.lmax) then  ! print only angular momentum included in the PAW data
       do isppol=1,dtset%nsppol
         write(tag_is,'(i1)')isppol
         if(pawfatbnd==1) then
           call int2char4(il-1,tag_il)
           ABI_CHECK((tag_il(1:1)/='#'),'Bug: string length too short!')
           tmpfil = trim(fildata)// &
&           '_at'//trim(tag_at)//'_'//trim(adjustl(symbol))//'_is'//tag_is//'_l'//trim(tag_il)
         else if (pawfatbnd==2) then
           write(tag_l,'(i1)') ll
           mm=il-(ll**2+ll+1)      ! compute m
           if(mm<0) write(tag_2m,'(i2)') mm
           if(mm>=0) write(tag_1m,'(i1)') mm
           if(mm<0) tmpfil = trim(fildata)// &
&           '_at'//trim(tag_at)//'_'//trim(adjustl(symbol))//'_is'//tag_is//'_l'//tag_l//'_m'//tag_2m
           if(mm>=0) tmpfil = trim(fildata)// &
&           '_at'//trim(tag_at)//'_'//trim(adjustl(symbol))//'_is'//tag_is//'_l'//tag_l//'_m+'//tag_1m
         end if
         !unitfatbands_arr(iall,isppol)=tmp_unit+100+iall-1+(natsph*inbfatbands)*(isppol-1)
         !open (unit=unitfatbands_arr(iall,isppol),file=trim(tmpfil),status='unknown',form='formatted')
         if (open_file(tmpfil, message, newunit=unitfatbands_arr(iall,isppol), status='unknown',form='formatted') /= 0) then
           MSG_ERROR(message)
         end if

         write(message,'(a,a,a,i4)') 'opened file : ', trim(tmpfil), ' unit', unitfatbands_arr(iall,isppol)
         call wrtout(std_out,message,'COLL')
         write(message,'(9a)') "# ",ch10,"# ABINIT package : FATBAND file ", ch10,&
&         "# It contains, for each band: the eigenvalues in eV (and the character of the band) as a function of the k-point",&
&         ch10,"# This file can be read with xmgrace (http://plasma-gate.weizmann.ac.il/Grace/)  ",ch10,"#  "
         write(unitfatbands_arr(iall,isppol), "(a)")trim(message)
         do iband=1,mband
           call int2char4(iband-1,tag_grace)
           ABI_CHECK((tag_grace(1:1)/='#'),'Bug: string length too short!')
           write(message,'(16a)') ch10,"@    s",trim(tag_grace)," line color 1",&
&           ch10,"@    s",trim(tag_grace)," errorbar color 2",&
&           ch10,"@    s",trim(tag_grace)," errorbar riser linewidth 5.0", &
&           ch10,"@    s",trim(tag_grace)," errorbar linestyle 0"
           write(unitfatbands_arr(iall,isppol), "(a)")trim(message)
         end do  !iband
         write(unitfatbands_arr(iall,isppol), '(a,a)') ch10,'@type xydy'
       end do   ! isppol
     end if ! ll=<lmax
   end do   ! il
 end do  ! iat

 if(iall.ne.(natsph*inbfatbands)) then
   MSG_ERROR("error1 ")
 end if



!--------------  WRITE FATBANDS IN FILES
 if (pawfatbnd>0) then
   ! Store eigenvalues with nkpt as first dimension for efficiency reasons
   ABI_ALLOCATE(eigenvalues,(nkpt,mband,dtset%nsppol))
   do isppol=1,dtset%nsppol
     do ikpt=1,nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
       do iband=1,mband
         eigenvalues(ikpt,iband,isppol)= ebands%eig(iband, ikpt, isppol) - ebands%fermie
       end do
     end do
   end do
   iall=0
   do iat=1,natsph
     lmax=(pawtab(dtset%typat(dtset%iatsph(iat)))%l_size-1)/2
     do il=1,inbfatbands
       iall=iall+1
       ll=int(sqrt(float(il-1)))
       if(ll.le.lmax) then
         do isppol=1,dtset%nsppol
           do iband=1,mband
             write(unitfatbands_arr(iall,isppol),'(a,a,i8)') ch10,"# BAND number :",iband
             do ikpt=1,nkpt
               if(pawfatbnd==1) then
                 xfatband=0.d0
                 do ixfat=(il-1)**2+1,il**2
                   xfatband=xfatband+dos%fractions_m(ikpt,iband,isppol,(iat-1)*mbesslang**2+ixfat)
                 end do ! ixfat
               else if (pawfatbnd==2) then
                 xfatband=dos%fractions_m(ikpt,iband,isppol,(iat-1)*mbesslang**2+il)
               end if
               write(unitfatbands_arr(iall,isppol),'(i5,e20.5,e20.5)')&
                 ikpt-1,eigenvalues(ikpt,iband,isppol)*Ha_eV,xfatband
             end do ! ikpt
           end do  !iband
           write(unitfatbands_arr(iall,isppol),'(a)') '&'
           !close(unitfatbands_arr(iall,isppol))
         end do  !isppol
       end if
     end do ! il
   end do ! iat
   ABI_DEALLOCATE(eigenvalues)
 end if

 do isppol=1,size(unitfatbands_arr, dim=2)
   do iat=1,size(unitfatbands_arr, dim=1)
     if (unitfatbands_arr(iat, isppol) /= -3) close (unitfatbands_arr(iat, isppol))
   end do
 end do

 ABI_DEALLOCATE(unitfatbands_arr)

 call cwtime(cpu,wall,gflops,"stop")
 write(message,'(2(a,f8.2),a)')" prtfatbands: cpu_time: ",cpu,"[s], walltime: ",wall," [s]"
 call wrtout(std_out,message,"PERS")

 DBG_EXIT("COLL")

end subroutine prtfatbands
!!***

!----------------------------------------------------------------------

!!****f* m_epjdos/fatbands_ncwrite
!! NAME
!! fatbands_ncwrite
!!
!! FUNCTION
!!  Write PJDOS contributions to netcdf file.
!!
!! INPUTS
!!  crystal<crystal_t>=Object defining the unit cell and its symmetries.
!!  ebands<ebands_t>=Band structure data.
!!  hdr<hdr_t>=Abinit header
!!  dtset<dtset_type>=Dataset type
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ncid=NC file handle.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      cwtime,wrtout
!!
!! SOURCE

subroutine fatbands_ncwrite(dos, crystal, ebands, hdr, dtset, psps, pawtab, ncid)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncid
 type(epjdos_t),intent(in) :: dos
 type(crystal_t),intent(in) :: crystal
 type(ebands_t),intent(in) :: ebands
 type(hdr_type),intent(in) :: hdr
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*dtset%usepaw)

#ifdef HAVE_NETCDF
!Local variables-------------------------------
!scalars
 integer :: itype,ncerr,fform
 real(dp) :: cpu,wall,gflops
 character(len=500) :: msg
!arrays
 integer :: lmax_type(crystal%ntypat)

!*************************************************************************

 ABI_CHECK(dtset%natsph > 0, "natsph <=0")
 call cwtime(cpu, wall, gflops, "start")

 fform = fform_from_ext("FATBANDS.nc")
 ABI_CHECK(fform /= 0, "Cannot find fform associated to FATBANDS.nc")

 ! Write header, crystal structure and band energies.
 NCF_CHECK(hdr_ncwrite(hdr, ncid, fform, nc_define=.True.))
 NCF_CHECK(crystal%ncwrite(ncid))
 NCF_CHECK(ebands_ncwrite(ebands, ncid))

 ! Add fatband-specific quantities
 ncerr = nctk_def_dims(ncid, [ &
   nctkdim_t("natsph", dtset%natsph), &
   nctkdim_t("ndosfraction", dos%ndosfraction)], defmode=.True.)
 NCF_CHECK(ncerr)

 if (dos%ndosfraction*dos%mbesslang > 0) then
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("mbesslang", dos%mbesslang), &
     nctkdim_t("dos_fractions_m_lastsize", dos%ndosfraction*dos%mbesslang)])
   NCF_CHECK(ncerr)
 end if
 if (dtset%natsph_extra /= 0) then
   NCF_CHECK(nctk_def_dims(ncid, [nctkdim_t("natsph_extra", dtset%natsph_extra)]))
 end if

 ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "prtdos", "pawprtdos", "prtdosm"])
 NCF_CHECK(ncerr)
 ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "ratsph_extra"])
 NCF_CHECK(ncerr)

 ncerr = nctk_def_arrays(ncid, [&
   nctkarr_t("lmax_type", "int", "number_of_atom_species"), &
   nctkarr_t("iatsph", "int", "natsph"), &
   nctkarr_t("ratsph", "dp", "number_of_atom_species"), &
   nctkarr_t("dos_fractions", "dp", "number_of_kpoints, max_number_of_states, number_of_spins, ndosfraction") &
 ])
 NCF_CHECK(ncerr)

 if (allocated(dos%fractions_m)) then
   ncerr = nctk_def_arrays(ncid, &
     nctkarr_t("dos_fractions_m", "dp", &
               "number_of_kpoints, max_number_of_states, number_of_spins, dos_fractions_m_lastsize"))
   NCF_CHECK(ncerr)
 end if

 if (allocated(dos%fractions_paw1)) then
   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t("dos_fractions_paw1", "dp", "number_of_kpoints, max_number_of_states, number_of_spins, ndosfraction"), &
     nctkarr_t("dos_fractions_pawt1", "dp", "number_of_kpoints, max_number_of_states, number_of_spins, ndosfraction") &
   ])
   NCF_CHECK(ncerr)
 end if

 if (dtset%natsph_extra /= 0) then
   ncerr = nctk_def_arrays(ncid, [&
     nctkarr_t("xredsph_extra", "dp", "number_of_reduced_dimensions, natsph_extra") &
   ])
   NCF_CHECK(ncerr)
 end if

 ! Write variables
 NCF_CHECK(nctk_set_datamode(ncid))

 ! scalars
 NCF_CHECK(nf90_put_var(ncid, vid("pawprtdos"), dtset%pawprtdos))
 NCF_CHECK(nf90_put_var(ncid, vid("prtdos"), dos%prtdos))
 NCF_CHECK(nf90_put_var(ncid, vid("prtdosm"), dos%prtdosm))

 ! arrays
 if (dtset%usepaw == 1) then
   lmax_type = (pawtab(:)%l_size - 1) / 2
 else
   do itype=1,crystal%ntypat
     lmax_type(itype) = maxval(psps%indlmn(1, :, itype))
   end do
 end if
 NCF_CHECK(nf90_put_var(ncid, vid("lmax_type"), lmax_type))
 NCF_CHECK(nf90_put_var(ncid, vid("dos_fractions"), dos%fractions))

 if (dos%prtdos == 3) then
   NCF_CHECK(nf90_put_var(ncid, vid("iatsph"), dtset%iatsph(1:dtset%natsph)))
   NCF_CHECK(nf90_put_var(ncid, vid("ratsph"), dtset%ratsph(1:dtset%ntypat)))
   NCF_CHECK(nf90_put_var(ncid, vid("ratsph_extra"), dtset%ratsph_extra))
   if (dtset%natsph_extra /= 0) then
     NCF_CHECK(nf90_put_var(ncid, vid("xredsph_extra"), dtset%xredsph_extra(:, 1:dtset%natsph_extra)))
   end if
 end if

 if (allocated(dos%fractions_m)) then
   NCF_CHECK(nf90_put_var(ncid, vid("dos_fractions_m"), dos%fractions_m))
 end if
 if (allocated(dos%fractions_paw1)) then
   NCF_CHECK(nf90_put_var(ncid, vid("dos_fractions_paw1"), dos%fractions_paw1))
   NCF_CHECK(nf90_put_var(ncid, vid("dos_fractions_pawt1"), dos%fractions_pawt1))
 end if

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2),a)')" fatbands_ncwrite: cpu_time: ",cpu,"[s], walltime: ",wall," [s]"
 call wrtout(std_out,msg,"PERS")
#endif

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
 end function vid

end subroutine fatbands_ncwrite
!!***

!!****f* m_epjdos/partial_dos_fractions
!! NAME
!! partial_dos_fractions
!!
!! FUNCTION
!! calculate partial DOS fractions to feed to the tetrahedron method
!!  1: project states on angular momenta
!!  2: should be able to choose certain atoms or atom types, slabs of space...
!!
!! INPUTS
!!  crystal<crystal_t>= data type gathering info on symmetries and unit cell
!!  dtset<type(dataset_type)>=all input variables for this dataset
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  cg(2,mcg)=planewave coefficients of wavefunctions
!!  mcg=size of wave-functions array (cg) =mpw*my_nspinor*mband*mkmem*nsppol
!!  collect=1 if fractions should be MPI collected at the end, 0 otherwise.
!!  mpi_enreg=information about MPI parallelization
!!
!! SIDE EFFECTS
!!  dos%fractions(ikpt,iband,isppol,ndosfraction) = percentage of s, p, d..
!!    character on each atom for the wavefunction # ikpt,iband, isppol
!!  == if prtdosm /= 0
!!  dos%fractions_m(ikpt,iband,isppol,ndosfraction*mbesslang) = percentage of s, p, d..
!!    character on each atom for the wavefunction # ikpt,iband, isppol (m-resolved)
!!
!! NOTES
!!
!!   psi(r) = (4pi/sqrt(ucvol)) \sum_{LMG} i**l u(G) e^{i(k+G).Ra} x Y_{LM}^*(k+G) Y_{LM}(r-Ra) j_L(|k+G||r-Ra|)
!!
!!   int_(ratsph) |psi(r)|**2 = \sum_LM rho(LM)
!!
!!   where
!!
!!   rho_{LM} = 4pi \int_o^{rc} dr r**2 ||\sum_G u(G) Y_LM^*(k+G) e^{i(k+G).Ra} j_L(|k+G| r)||**2
!!
!!   where S is a RSH. The final expression is:
!!
!!   When k = G0/2, we have u_{G0/2}(G) = u_{G0/2}(-G-G0)^* and P can be rewritten as
!!
!!     P = (4pi i^L}/sqrt(ucvol) \sum^'_G w(G) S_{LM}(k+G) \int_0^ratsph dr r^2 j_L(|k+G|r) x
!!                                  2 Re[u_k(G) e^{i(k+G).R_atom}]  if L = 2n
!!                                  2 Im[u_k(G) e^{i(k+G).R_atom}]  if L = 2n + 1
!!
!!  where the sum over G is done on the reduced G-sphere and w(G) = 1/2 if G=G0 else 1.
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      cg_getspin,cwtime,destroy_mpi_enreg,getkpgnorm,getph,initmpi_seq
!!      initylmg,jlspline_free,ph1d3d,recip_ylm,sort_dp,splint,wrtout,xmpi_sum
!!
!! SOURCE

subroutine partial_dos_fractions(dos,crystal,dtset,eigen,occ,npwarr,kg,cg,mcg,collect,mpi_enreg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcg,collect
 type(epjdos_t),intent(inout) :: dos
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: crystal
!arrays
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),npwarr(dtset%nkpt)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: prtsphere0=0 ! do not output all the band by band details for projections.
 integer :: shift_b,shift_sk,iat,iatom,iband,ierr,ikpt,ilang,ioffkg,is1, is2, isoff
 integer :: ipw,isppol,ixint,mbess,mcg_disk,me_kpt,shift_cg
 integer :: mgfft,my_nspinor,n1,n2,n3,natsph_tot,npw_k,nradintmax
 integer :: rc_ylm,itypat,nband_k
 integer :: abs_shift_b
 integer :: unit_procar, ipauli
 real(dp),parameter :: bessint_delta = 0.1_dp
 real(dp) :: arg,bessarg,bessargmax,kpgmax,rmax
 real(dp) :: cpu,wall,gflops
 character(len=500) :: msg
 character(len=4) :: ikproc_str
 character(len=fnlen) :: filename
 type(jlspline_t) :: jlspl
 type(MPI_type) :: mpi_enreg_seq
!arrays
 integer :: iindex(dtset%mpw),nband_tmp(1),npwarr_tmp(1)
 integer,allocatable :: iatsph(:),nradint(:),atindx(:),typat_extra(:),kg_k(:,:)
 real(dp) :: kpoint(3),spin(3),ylmgr_dum(1)
 real(dp) :: xfit(dtset%mpw),yfit(dtset%mpw)
 real(dp),allocatable :: ylm_k(:,:)
 real(dp),allocatable :: bess_fit(:,:,:)
 real(dp),allocatable :: cg_1band(:,:),cg_1kpt(:,:),kpgnorm(:),ph1d(:,:)
 real(dp),allocatable :: ph3d(:,:,:),ratsph(:),rint(:),sum_1ll_1atom(:,:,:)
 real(dp),allocatable :: sum_1lm_1atom(:,:,:)
 real(dp),allocatable :: cplx_1lm_1atom(:,:,:,:)
 real(dp),allocatable :: xred_sph(:,:),znucl_sph(:),phkxred(:,:)
 complex(dpc) :: cgcmat(2,2)

!*************************************************************************

 ! for the moment, only support projection on angular momenta
 if (dos%partial_dos_flag /= 1 .and. dos%partial_dos_flag /= 2) then
   write(std_out,*) 'Error: partial_dos_fractions only supports angular '
   write(std_out,*) ' momentum projection and spinor components for the moment. return to outscfcv'
   write(std_out,*) ' partial_dos = ', dos%partial_dos_flag
   return
 end if

 ! impose all kpoints have same number of bands
 do isppol=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
     if (dtset%nband((isppol-1)*dtset%nkpt + ikpt) /= dtset%mband) then
       write(std_out,*) 'Error: partial_dos_fractions wants same number of bands at each kpoint'
       write(std_out,*) ' isppol, ikpt = ', isppol,ikpt, dtset%nband((isppol-1)*dtset%nkpt + ikpt), dtset%mband
       write(std_out,*) ' all nband = ', dtset%nband
       return
     end if
   end do
 end do

 ! Real or complex spherical harmonics?
 rc_ylm = 2; if (dos%prtdosm == 2) rc_ylm = 1

 me_kpt = mpi_enreg%me_kpt
 my_nspinor = max(1,dtset%nspinor/mpi_enreg%nproc_spinor)

 n1 = dtset%ngfft(1); n2 = dtset%ngfft(2); n3 = dtset%ngfft(3)
 mgfft = maxval(dtset%ngfft(1:3))

 call cwtime(cpu, wall, gflops, "start")

 if (dtset%prtprocar /= 0) then
   ! open file for each proc, and print header for master node
   call int2char4(me_kpt, ikproc_str)
   filename = 'PROCAR_'//ikproc_str
   if (open_file(filename, msg, newunit=unit_procar, form="formatted", action="write") /= 0) then
      MSG_ERROR(msg)
   end if
   if(mpi_enreg%me==0) then
     write (unit_procar,'(a)') 'PROCAR lm decomposed - need to merge files yourself in parallel case!!! Or use pyprocar package'
     if (dtset%prtprocar == 2) then
       write (unit_procar,'(a)') ' Requested complex output of PROCAR file (prtprocar 2)' 
     end if
     write (unit_procar,'(a,I10,a,I10,a,I10,a)') '# of k-points: ', dtset%nkpt, &
       ' # of bands:', dtset%mband, ' # of ions:', dtset%natom, ch10
   end if
 end if

!##############################################################
!FIRST CASE: project on angular momenta to get dos parts
!##############################################################

 if (dos%partial_dos_flag == 1) then
   natsph_tot = dtset%natsph + dtset%natsph_extra

   ABI_ALLOCATE(iatsph, (natsph_tot))
   ABI_ALLOCATE(typat_extra, (natsph_tot))
   ABI_ALLOCATE(ratsph, (natsph_tot))
   ABI_ALLOCATE(znucl_sph, (natsph_tot))
   ABI_ALLOCATE(nradint, (natsph_tot))
   ABI_ALLOCATE(atindx, (natsph_tot))
   ABI_ALLOCATE(phkxred, (2,natsph_tot))

   ! initialize atindx
   do iatom=1,natsph_tot
     atindx(iatom) = iatom
   end do

   iatsph(1:dtset%natsph) = dtset%iatsph(1:dtset%natsph)
   do iatom=1,dtset%natsph
     itypat = dtset%typat(iatsph(iatom))
     typat_extra(iatom) = itypat
     ratsph(iatom) = dtset%ratsph(itypat)
     znucl_sph(iatom) = dtset%znucl(itypat)
   end do

   ! fictitious atoms are declared with
   ! %natsph_extra, %ratsph_extra and %xredsph_extra(3, dtset%natsph_extra)
   ! they have atom index (natom + ii) and itype = ntypat + 1
   do iatom=1,dtset%natsph_extra
     typat_extra(iatom+dtset%natsph) = dtset%ntypat + 1
     ratsph(iatom+dtset%natsph) = dtset%ratsph_extra
     znucl_sph(iatom+dtset%natsph) = zero
     iatsph(iatom+dtset%natsph) = dtset%natom + iatom
   end do

   ! init bessel function integral for recip_ylm max ang mom + 1
   ABI_ALLOCATE(sum_1ll_1atom,(dtset%nspinor**2,dos%mbesslang,natsph_tot))
   ABI_ALLOCATE(sum_1lm_1atom,(dtset%nspinor**2,dos%mbesslang**2,natsph_tot))
   ABI_ALLOCATE(cplx_1lm_1atom,(2,dtset%nspinor,dos%mbesslang**2,natsph_tot))

   ! Note ecuteff instead of ecut.
   kpgmax = sqrt(dtset%ecut * dtset%dilatmx**2)
   rmax = zero; bessargmax = zero; nradintmax = 0
   do iatom=1,natsph_tot
     rmax = max(rmax, ratsph(iatom))
     bessarg = ratsph(iatom)*two_pi*kpgmax
     bessargmax = max(bessargmax, bessarg)
     nradint(iatom) = int (bessarg / bessint_delta) + 1
     nradintmax = max(nradintmax,nradint(iatom))
   end do
   !write(std_out,*)' partial_dos_fractions: rmax=', rmax,' nradintmax: ", nradintmax
!  use same number of grid points to calculate Bessel function and to do the integration later on r
!  and make sure bessargmax is a multiple of bessint_delta
   mbess = nradintmax
   bessargmax = bessint_delta*mbess

   ABI_ALLOCATE(rint,(nradintmax))
   ABI_ALLOCATE(bess_fit,(dtset%mpw,nradintmax,dos%mbesslang))

   ! initialize general Bessel function array on uniform grid xx, from 0 to (2 \pi |k+G|_{max} |r_{max}|)
   jlspl = jlspline_new(mbess, bessint_delta, dos%mbesslang)

   ABI_ALLOCATE(xred_sph, (3, natsph_tot))
   do iatom=1,dtset%natsph
     xred_sph(:,iatom) = crystal%xred(:,iatsph(iatom))
   end do
   do iatom=1,dtset%natsph_extra
     xred_sph(:,dtset%natsph+iatom) = dtset%xredsph_extra(:, iatom)
   end do

   ABI_ALLOCATE(ph1d,(2,(2*n1+1 + 2*n2+1 + 2*n3+1)*natsph_tot))
   call getph(atindx,natsph_tot,n1,n2,n3,ph1d,xred_sph)

   ! Fake MPI data to be used for sequential call to initylmg.
   call initmpi_seq(mpi_enreg_seq)
   mpi_enreg_seq%my_natom = dtset%natom

   shift_sk = 0
   abs_shift_b =  0 ! offset to allow for automatic update with +1 below
   do isppol=1,dtset%nsppol
     ioffkg = 0
     do ikpt=1,dtset%nkpt
       nband_k = dtset%nband((isppol-1)*dtset%nkpt + ikpt)
       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_kpt)) then
         abs_shift_b = abs_shift_b + nband_k ! jump the whole kpt in the eig and occ arrays
         cycle
       end if
       npw_k = npwarr(ikpt)
       kpoint(:) = dtset%kpt(:,ikpt)

       if (dtset%prtprocar /= 0) then
         write (unit_procar,'(a,I7,a,3F12.6,a,F12.6,a)') &
           ' k-point ', ikpt, ' : ', kpoint(:), ' weight = ', dtset%wtk(ikpt), ch10
       end if

       ! make phkred for all atoms
       do iat=1,natsph_tot
         arg=two_pi*( kpoint(1)*xred_sph(1,iat) + kpoint(2)*xred_sph(2,iat) + kpoint(3)*xred_sph(3,iat) )
         phkxred(1,iat)=cos(arg)
         phkxred(2,iat)=sin(arg)
       end do

       ABI_MALLOC(kg_k, (3, npw_k))
       kg_k = kg(:,ioffkg+1:ioffkg+npw_k)

       ! kpgnorm contains norms only for kpoints used by this processor
       ABI_ALLOCATE(kpgnorm, (npw_k))
       call getkpgnorm(crystal%gprimd, kpoint, kg_k, kpgnorm, npw_k)

       ! Now get Ylm(k, G) factors: returns "real Ylms", which are real (+m) and
       ! imaginary (-m) parts of actual complex Ylm. Yl-m = Ylm*
       ! Call initylmg for a single k-point (mind mpi_enreg_seq).
       ABI_MALLOC(ylm_k, (npw_k, dos%mbesslang**2))
       npwarr_tmp(1) = npw_k; nband_tmp(1) = nband_k
       call initylmg(crystal%gprimd,kg_k,kpoint,1,mpi_enreg_seq,dos%mbesslang,&
       npw_k,nband_tmp,1,npwarr_tmp,1,0,crystal%rprimd,ylm_k,ylmgr_dum)

       ! get phases exp (2 pi i (k+G).x_tau) in ph3d
       ABI_ALLOCATE(ph3d,(2,npw_k,natsph_tot))
       call ph1d3d(1,natsph_tot,kg_k,natsph_tot,natsph_tot,npw_k,n1,n2,n3,phkxred,ph1d,ph3d)

       ! get Bessel function factors on array of |k+G|*r distances
       ! since we need many r distances and have a large number of different
       ! |k+G|, get j_l on uniform grid (above, in array gen_besj),
       ! and spline it for each kpt Gvector set.
       ! Note that we use the same step based on rmax, this can lead to (hopefully small)
       ! inaccuracies when we integrate from 0 up to rmax(iatom)
       do ixint=1,nradintmax
         rint(ixint) = (ixint-1)*rmax / (nradintmax-1)
         do ipw=1,npw_k
           xfit(ipw) = two_pi * kpgnorm(ipw) * rint(ixint)
           iindex(ipw) = ipw
         end do

         call sort_dp(npw_k,xfit,iindex,tol14)
         do ilang=1,dos%mbesslang
           call splint(mbess, jlspl%xx, jlspl%bess_spl(:,ilang), jlspl%bess_spl_der(:,ilang), npw_k, xfit, yfit)
           ! re-order results for different G vectors
           do ipw=1,npw_k
             bess_fit(iindex(ipw),ixint,ilang) = yfit(ipw)
           end do
         end do
       end do ! ixint

       shift_b = 0
       do iband=1,nband_k
         abs_shift_b = abs_shift_b + 1
         if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband,isppol,me_kpt)) cycle
         !write(std_out,*)"in band:",iband
         ! TODO: eventually import eig and occ down to here - a pain, but printing outside would imply saving a huge array in memory
         if (dtset%prtprocar /= 0) then
           write (unit_procar,'(a,I7,a,F12.6,a,F12.6,a)') 'band ', iband, ' # energy ', &
             eigen(abs_shift_b), ' # occ. ', occ(abs_shift_b), ch10
         end if

         ! Select wavefunction in cg array
         shift_cg = shift_sk + shift_b

         call recip_ylm(bess_fit, cg(:,shift_cg+1:shift_cg+my_nspinor*npw_k), dtset%istwfk(ikpt),&
&          mpi_enreg, nradint, nradintmax, dos%mbesslang , dtset%mpw, natsph_tot, typat_extra, dos%mlang_type,&
&          npw_k, dtset%nspinor, ph3d, prtsphere0, rint, ratsph, rc_ylm, sum_1ll_1atom, sum_1lm_1atom, cplx_1lm_1atom,&
&          crystal%ucvol, ylm_k, znucl_sph)
         ! on exit the sum_1atom_* have both spinors counted

         ! Accumulate
         do iatom=1,natsph_tot
           do ilang=1,dos%mbesslang
             dos%fractions(ikpt,iband,isppol,dos%mbesslang*(iatom-1) + ilang) &
&             = dos%fractions(ikpt,iband,isppol,dos%mbesslang*(iatom-1) + ilang) &
&             + sum_1ll_1atom(1,ilang,iatom)
           end do
         end do

         if (dos%prtdosm /= 0) then
           do iatom=1,natsph_tot
             do ilang=1,dos%mbesslang**2
               dos%fractions_m(ikpt,iband,isppol,dos%mbesslang**2*(iatom-1) + ilang) &
&               = dos%fractions_m(ikpt,iband,isppol,dos%mbesslang**2*(iatom-1) + ilang) &
&               + sum_1lm_1atom(1,ilang,iatom)
             end do
           end do
         end if

         ! Increment band, spinor shift
         !shift_b = shift_b + npw_k
         shift_b = shift_b + my_nspinor*npw_k

         ! now we have both spinor components. The first option is for real values projections, eventually decomposed by Pauli spinor components
         if (dtset%prtprocar /= 0) then
           write (unit_procar,'(a)') 'ion       s      py     pz     px    dxy    dyz    dz2    dxz    dx2    tot'
           do ipauli= 1,dtset%nspinor**2
             ! Contract with Pauli matrices to get projections for this k and band, all atoms and ilang
             do iatom = 1, natsph_tot
               write (unit_procar, '(1x,I5)', advance='no') iatom
               do ilang=1,min(dos%mbesslang**2,9)
                 write (unit_procar, '(F7.3)',advance='no') sum_1lm_1atom(ipauli,ilang,iatom)
               end do
               write (unit_procar, '(F7.3)',advance='yes') sum(sum_1lm_1atom(ipauli,:,iatom))
             end do
             ! final line with sum over atoms
             write (unit_procar, '(a)', advance='no') 'tot   '
             do ilang=1,min(dos%mbesslang**2,9)
               write (unit_procar, '(F7.3)',advance='no') sum(sum_1lm_1atom(ipauli,ilang,:))
             end do
             write (unit_procar, '(F7.3)',advance='yes') sum(sum_1lm_1atom(ipauli,:,:))
           end do

           ! second option is to also print the complex projection on the atomic like orbital: <psi_nk | Y_lm> in a sphere
           !  Two blocks are printed, first real then imaginary part
           if (dtset%prtprocar == 2) then
             write (unit_procar,'(2a)') 'ion            s              py              pz              px',&
&              '             dxy             dyz             dz2             dxz             dx2             tot'
             do is1= 1,dtset%nspinor
               ! Contracted with Pauli matrices to get projections for this k and band, all atoms and ilang
               do iatom = 1, natsph_tot
                 write (unit_procar, '(1x,I5)', advance='no') iatom
                 do ilang=1,min(dos%mbesslang**2,9)
                   write (unit_procar, '(2(F7.3,1x))',advance='no') cplx_1lm_1atom(:,is1,ilang,iatom)
                 end do
                 write (unit_procar, '(2(F7.3,1x))',advance='yes') sum(cplx_1lm_1atom(1,is1,:,iatom)), &
&                   sum(cplx_1lm_1atom(2,is1,:,iatom))
               end do
               ! final line with sum over atoms
               write (unit_procar, '(a)', advance='no') 'charge'
               do ilang=1,min(dos%mbesslang**2,9)
                 write (unit_procar, '(F7.3,9x)',advance='no') sum(cplx_1lm_1atom(:,is1,ilang,:)**2)
               end do
               write (unit_procar, '(F7.3,9x)',advance='yes') sum(cplx_1lm_1atom(:,is1,:,:)**2)
             end do
             write (unit_procar,*)
           end if
         end if

       end do ! band

       ! Increment kpt and (spin, kpt) shifts
       ioffkg = ioffkg + npw_k
       shift_sk = shift_sk + nband_k*my_nspinor*npw_k

       ABI_FREE(kg_k)
       ABI_FREE(kpgnorm)
       ABI_FREE(ylm_k)
       ABI_DEALLOCATE(ph3d)
     end do ! ikpt
   end do ! isppol

   ! collect = 1 ==> gather all contributions from different processors
   if (collect == 1) then
     call xmpi_sum(dos%fractions,mpi_enreg%comm_kpt,ierr)
     if (dos%prtdosm /= 0) call xmpi_sum(dos%fractions_m,mpi_enreg%comm_kpt,ierr)

! this is now done inside recip_ylm
!     if (mpi_enreg%paral_spinor == 1)then
!       call xmpi_sum(dos%fractions,mpi_enreg%comm_spinor,ierr)
!       if (dos%prtdosm /= 0) call xmpi_sum(dos%fractions_m,mpi_enreg%comm_spinor,ierr)
!     end if
   end if

   ABI_DEALLOCATE(atindx)
   ABI_DEALLOCATE(bess_fit)
   ABI_DEALLOCATE(iatsph)
   ABI_DEALLOCATE(typat_extra)
   ABI_DEALLOCATE(nradint)
   ABI_DEALLOCATE(ph1d)
   ABI_DEALLOCATE(phkxred)
   ABI_DEALLOCATE(ratsph)
   ABI_DEALLOCATE(rint)
   ABI_DEALLOCATE(sum_1ll_1atom)
   ABI_DEALLOCATE(sum_1lm_1atom)
   ABI_DEALLOCATE(cplx_1lm_1atom)
   ABI_DEALLOCATE(xred_sph)
   ABI_DEALLOCATE(znucl_sph)

   call jlspline_free(jlspl)
   call destroy_mpi_enreg(mpi_enreg_seq)

 !##############################################################
 !2ND CASE: project on spinors
 !##############################################################

 else if (dos%partial_dos_flag == 2) then

   if (dtset%nsppol /= 1 .or. dtset%nspinor /= 2) then
     MSG_WARNING("spinor projection is meaningless if nsppol==2 or nspinor/=2. Not calculating projections.")
     return
   end if
   if (my_nspinor /= 2) then
     MSG_WARNING("spinor projection with spinor parallelization is not coded. Not calculating projections.")
     return
   end if
   ABI_CHECK(mpi_enreg%paral_spinor == 0, "prtdos 5 does not support spinor parallelism")

   ! FIXME: We should not allocate such a large chunk of memory!
   mcg_disk = dtset%mpw*my_nspinor*dtset%mband
   ABI_ALLOCATE(cg_1kpt,(2,mcg_disk))
   shift_sk = 0
   isppol = 1

   do ikpt=1,dtset%nkpt
     nband_k = dtset%nband((isppol-1)*dtset%nkpt + ikpt)
     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_kpt)) cycle
     npw_k = npwarr(ikpt)

     cg_1kpt(:,:) = cg(:,shift_sk+1:shift_sk+mcg_disk)
     ABI_ALLOCATE(cg_1band,(2,2*npw_k))
     shift_b=0
     do iband=1,nband_k
       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband,isppol,me_kpt)) cycle

       ! Select wavefunction in cg array
       !shift_cg = shift_sk + shift_b
       cg_1band(:,:) = cg_1kpt(:,shift_b+1:shift_b+2*npw_k)
       call cg_getspin(cg_1band, npw_k, spin, cgcmat=cgcmat)

       ! MG: TODO: imag part of off-diagonal terms is missing.
       ! I will add them later on.
       do is1=1,2
         do is2=1,2
           isoff = is2 + (is1-1)*2
           dos%fractions(ikpt,iband,isppol,isoff) = dos%fractions(ikpt,iband,isppol,isoff) &
&           + real(cgcmat(is1,is2))
         end do
       end do

       dos%fractions(ikpt,iband,isppol,5) = dos%fractions(ikpt,iband,isppol,5) + spin(1)
       dos%fractions(ikpt,iband,isppol,6) = dos%fractions(ikpt,iband,isppol,6) + spin(2)
       dos%fractions(ikpt,iband,isppol,7) = dos%fractions(ikpt,iband,isppol,7) + spin(3)

       shift_b = shift_b + 2*npw_k
     end do
     ABI_DEALLOCATE(cg_1band)
     shift_sk = shift_sk + nband_k*2*npw_k
   end do
   ABI_DEALLOCATE(cg_1kpt)

   ! Gather all contributions from different processors
   if (collect == 1) then
     call xmpi_sum(dos%fractions,mpi_enreg%comm_kpt,ierr)
     call xmpi_sum(dos%fractions,mpi_enreg%comm_bandfft,ierr)
     !below for future use - spinors should not be parallelized for the moment
     !if (mpi_enreg%paral_spinor == 1)then
     !  call xmpi_sum(dos%fractions,mpi_enreg%comm_spinor,ierr)
     !end if
   end if

 else
   MSG_WARNING('only partial_dos==1 or 2 is coded')
 end if

 if (dtset%prtprocar /= 0) close(unit_procar)

 call cwtime(cpu,wall,gflops,"stop")
 write(msg,'(2(a,f8.2),a)')" partial_dos_fractions: cpu_time: ",cpu,"[s], walltime: ",wall," [s]"
 call wrtout(std_out,msg,"PERS")

end subroutine partial_dos_fractions
!!***

!!****f* m_epjdos/partial_dos_fractions_paw
!! NAME
!! partial_dos_fractions_paw
!!
!! FUNCTION
!!  Calculate PAW contributions to the partial DOS fractions (tetrahedron method)
!!
!! INPUTS
!!  cprj(natom,mcprj)= <p_lmn|Cnk> coefficients for each WF |Cnk> and each |p_lmn> non-local projector
!!  dimcprj(natom)=array of dimensions of array cprj (not ordered)
!!  dtset     structured datatype, from which one uses :
!!   iatsph(nasph)=number of atoms used to project dos
!!   kpt(3,nkpt)  =irreducible kpoints
!!   mband        =max number of bands per k-point
!!   mkmem        =number of kpoints in memory
!!   natom        =number of atoms in total
!!   natsph       =number of atoms ofor which the spherical decomposition must be done
!!   nband        =number of electronic bands for each kpoint
!!   nkpt         =number of irreducible kpoints
!!   nspinor      =number of spinor components
!!   nsppol       =1 or 2 spin polarization channels
!!  fatbands_flag =1 if pawfatbnd=1 or 2
!!  mbesslang=maximum angular momentum for Bessel function expansion
!!  mpi_enreg=informations about MPI parallelization
!!  prtdosm=option for the m-contributions to the partial DOS
!!  ndosfraction=natsph*mbesslang
!!  paw_dos_flag=option for the PAW contributions to the partial DOS
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data:
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  === If paw_dos_flag==1:
!!   dos%fractions_paw1(ikpt,iband,isppol,natom*mbesslang) = contribution to
!!       dos fractions from the PAW partial waves (phi)
!!   dos%fractions_pawt1(ikpt,iband,isppol,natom*mbesslang) = contribution to
!!       dos fractions from the PAW pseudo partial waves (phi_tild)
!!
!! SIDE EFFECTS
!!  dos%fractions(ikpt,iband,isppol,ndosfraction) = percentage of s, p, d..
!!    character on each atom for the wavefunction # ikpt,iband, isppol
!!    As input: contains only the pseudo contribution
!!    As output: contains pseudo contribution + PAW corrections
!!  == if prtdosm==1
!!  dos%fractions_m(ikpt,iband,isppol,ndosfraction*mbesslang*prtdosm) =
!!              m discretization of partial DOS fractions
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      pawcprj_alloc,pawcprj_free,simp_gen,timab,xmpi_sum
!!
!! SOURCE

subroutine partial_dos_fractions_paw(dos,cprj,dimcprj,dtset,mcprj,mkmem,mpi_enreg,pawrad,pawtab)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcprj,mkmem
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(epjdos_t),intent(inout) :: dos
!arrays
 integer,intent(in) :: dimcprj(dtset%natom)
 type(pawcprj_type),intent(in) :: cprj(dtset%natom,mcprj)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat)
 type(pawtab_type),target,intent(in) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer :: bandpp,basis_size,comm_kptband,cplex,fatbands_flag,iat,iatom,iband,ibg,ibsp
 integer :: ierr,ikpt,il,ilang,ilmn,iln,im,iorder_cprj,ispinor,isppol,itypat,j0lmn,j0ln
 integer :: jl,jlmn,jln,jm,klmn,kln,lmn_size,mbesslang,me_band,me_kpt,my_nspinor
 integer :: nband_cprj_k,nband_k,ndosfraction,nprocband,nproc_spkptband,paw_dos_flag,prtdosm
 real(dp) :: cpij,one_over_nproc
 !character(len=500) :: message
!arrays
 integer ,allocatable :: dimcprj_atsph(:)
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: int1(:,:),int2(:,:),int1m2(:,:)
 type(pawcprj_type),allocatable :: cprj_k(:,:)

!******************************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(mkmem/=0,"mkmem==0 not supported anymore!")

 fatbands_flag = dos%fatbands_flag
 mbesslang = dos%mbesslang
 prtdosm = dos%prtdosm
 ndosfraction = dos%ndosfraction
 paw_dos_flag = dos%paw_dos_flag

!m-decomposed DOS not compatible with PAW-decomposed DOS
 if(prtdosm>=1.and.paw_dos_flag==1) then
   MSG_ERROR('m-decomposed DOS not compatible with PAW-decomposed DOS!')
 end if

!Prepare some useful integrals
 basis_size=pawtab(1)%basis_size
 if (dtset%ntypat>1) then
   do itypat=1,dtset%ntypat
     basis_size=max(basis_size,pawtab(itypat)%basis_size)
   end do
 end if
 ABI_ALLOCATE(int1  ,(basis_size*(basis_size+1)/2,dtset%natsph))
 ABI_ALLOCATE(int2,(basis_size*(basis_size+1)/2,dtset%natsph))
 ABI_ALLOCATE(int1m2,(basis_size*(basis_size+1)/2,dtset%natsph))
 int1=zero;int2=zero;int1m2=zero
 do iat=1,dtset%natsph
   iatom=dtset%iatsph(iat)
   itypat= dtset%typat(iatom)
   do jln=1,pawtab(itypat)%basis_size
     j0ln=jln*(jln-1)/2
     do iln=1,jln
       kln=j0ln+iln
       call simp_gen(int1(kln,iat),pawtab(itypat)%phiphj(:,kln),pawrad(itypat))
       if (dtset%pawprtdos<2) then
         call simp_gen(int2(kln,iat),pawtab(itypat)%tphitphj(:,kln),pawrad(itypat))
         int1m2(kln,iat)=int1(kln,iat)-int2(kln,iat)
       else
         int2(kln,iat)=zero;int1m2(kln,iat)=int1(kln,iat)
       end if
     end do !iln
   end do !jln
 end do

!Antiferro case
 if (dtset%nspden==2.and.dtset%nsppol==1.and.dtset%nspinor==1) then
   int1m2(:,:)=half*int1m2(:,:)
   if (paw_dos_flag==1.or.fatbands_flag==1.or.prtdosm==2) then
     int1(:,:)=half*int1(:,:);int2(:,:)=half*int2(:,:)
   end if
 end if

!Init parallelism
 comm_kptband=mpi_enreg%comm_kptband
 nproc_spkptband=xmpi_comm_size(comm_kptband)*mpi_enreg%nproc_spinor
 me_kpt=mpi_enreg%me_kpt ; me_band=mpi_enreg%me_band
 my_nspinor=max(1,dtset%nspinor/mpi_enreg%nproc_spinor)
 bandpp=1;if (mpi_enreg%paral_kgb==1) bandpp=mpi_enreg%bandpp
!Check if cprj is distributed over bands
 nprocband=my_nspinor*dtset%mband*dtset%mkmem*dtset%nsppol/mcprj
 if (nprocband/=mpi_enreg%nproc_band) then
   MSG_BUG('wrong mcprj/nproc_band!')
 end if

!Quick hack: in case of parallelism, dos_fractions have already
!  been reduced over MPI processes; they have to be prepared before
!  the next reduction (at the end of the following loop).
 if (nproc_spkptband>1) then
   one_over_nproc=one/real(nproc_spkptband,kind=dp)
!$OMP  PARALLEL DO COLLAPSE(4) &
!$OMP& DEFAULT(SHARED) PRIVATE(ilang,isppol,iband,ikpt)
   do ilang=1,ndosfraction
     do isppol=1,dtset%nsppol
       do iband=1,dtset%mband
         do ikpt=1,dtset%nkpt
           dos%fractions(ikpt,iband,isppol,ilang)= &
&           one_over_nproc*dos%fractions(ikpt,iband,isppol,ilang)
         end do
       end do
     end do
   end do
!$OMP END PARALLEL DO
   if (fatbands_flag==1.or.prtdosm==1.or.prtdosm==2) then
!$OMP  PARALLEL DO COLLAPSE(4) &
!$OMP& DEFAULT(SHARED) PRIVATE(ilang,isppol,iband,ikpt)
     do ilang=1,ndosfraction*mbesslang
       do isppol=1,dtset%nsppol
         do iband=1,dtset%mband
           do ikpt=1,dtset%nkpt
             dos%fractions_m(ikpt,iband,isppol,ilang)= &
&             one_over_nproc*dos%fractions_m(ikpt,iband,isppol,ilang)
           end do
         end do
       end do
     end do
!$OMP END PARALLEL DO
   end if
 end if

 iorder_cprj=0

!LOOPS OVER SPINS,KPTS
 ibg=0
 do isppol =1,dtset%nsppol
   do ikpt=1,dtset%nkpt

     nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_kpt)) cycle

     cplex=2;if (dtset%istwfk(ikpt)>1) cplex=1
     nband_cprj_k=nband_k/nprocband
     ABI_DATATYPE_ALLOCATE(cprj_k,(dtset%natsph,my_nspinor*nband_cprj_k))
     ABI_ALLOCATE(dimcprj_atsph,(dtset%natsph))
     do iat=1,dtset%natsph
       dimcprj_atsph(iat)=dimcprj(dtset%iatsph(iat))
     end do
     call pawcprj_alloc(cprj_k,0,dimcprj_atsph)
     ABI_DEALLOCATE(dimcprj_atsph)

!    Extract cprj for this k-point.
     ibsp=0
     do iband=1,nband_cprj_k
       do ispinor=1,my_nspinor
         ibsp=ibsp+1
         do iat=1,dtset%natsph
           iatom=dtset%iatsph(iat)
           cprj_k(iat,ibsp)%cp(:,:)=cprj(iatom,ibsp+ibg)%cp(:,:)
         end do
       end do
     end do

!    LOOP OVER ATOMS (natsph_extra is not included on purpose)
     do iat=1,dtset%natsph
       iatom=dtset%iatsph(iat)
       itypat= dtset%typat(iatom)
       lmn_size=pawtab(itypat)%lmn_size
       indlmn => pawtab(itypat)%indlmn

!      LOOP OVER BANDS
       ibsp=0
       do iband=1,nband_k

         if (mod((iband-1)/bandpp,nprocband)/=me_band) cycle
         if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband,isppol,me_kpt)) then
           ibsp=ibsp+my_nspinor;cycle
         end if

         do ispinor=1,my_nspinor
           ibsp=ibsp+1

           do ilang=1,mbesslang

             do jlmn=1,lmn_size
               jl=indlmn(1,jlmn)
               jm=indlmn(2,jlmn)
               j0lmn=jlmn*(jlmn-1)/2
               do ilmn=1,jlmn
                 il=indlmn(1,ilmn)
                 im=indlmn(2,ilmn)
                 klmn=j0lmn+ilmn
                 kln=pawtab(itypat)%indklmn(2,klmn)

                 if (il==ilang-1.and.jl==ilang-1.and.im==jm) then

                   cpij=cprj_k(iat,ibsp)%cp(1,ilmn)*cprj_k(iat,ibsp)%cp(1,jlmn)
                   if (cplex==2) cpij=cpij+cprj_k(iat,ibsp)%cp(2,ilmn)*cprj_k(iat,ibsp)%cp(2,jlmn)
                   cpij=pawtab(itypat)%dltij(klmn)*cpij

                   dos%fractions(ikpt,iband,isppol,mbesslang*(iat-1)+ilang)=  &
&                   dos%fractions(ikpt,iband,isppol,mbesslang*(iat-1)+ilang) + &
&                   cpij*int1m2(kln,iat)
                   if (prtdosm==1) then
                     dos%fractions_m(ikpt,iband,isppol,mbesslang**2*(iat-1)+il**2+il+1+im)= &
&                     dos%fractions_m(ikpt,iband,isppol,mbesslang**2*(iat-1)+il**2+il+1+im) + &
&                     cpij*int1m2(kln,iat)
                   end if
                   if (fatbands_flag==1.or.prtdosm==2) then
                     dos%fractions_m(ikpt,iband,isppol,mbesslang**2*(iat-1)+il**2+il+1+im)= &
&                     dos%fractions_m(ikpt,iband,isppol,mbesslang**2*(iat-1)+il**2+il+1+im) + &
&                     cpij*int1(kln,iat)
                   end if
                   if (paw_dos_flag==1) then
                     dos%fractions_paw1(ikpt,iband,isppol,mbesslang*(iat-1)+ilang)=  &
&                     dos%fractions_paw1(ikpt,iband,isppol,mbesslang*(iat-1)+ilang) + &
&                     cpij*int1(kln,iat)
                     dos%fractions_pawt1(ikpt,iband,isppol,mbesslang*(iat-1)+ilang)=  &
&                     dos%fractions_pawt1(ikpt,iband,isppol,mbesslang*(iat-1)+ilang) + &
&                     cpij*int2(kln,iat)
                   end if

                 end if

               end do !ilmn
             end do   !jlmn

           end do ! ilang
         end do ! ispinor
       end do ! iband

     end do !iatom

     if (mkmem/=0) ibg = ibg + my_nspinor*nband_cprj_k
     call pawcprj_free(cprj_k)
     ABI_DATATYPE_DEALLOCATE(cprj_k)
   end do ! ikpt
 end do ! isppol

 ABI_DEALLOCATE(int1)
 ABI_DEALLOCATE(int2)
 ABI_DEALLOCATE(int1m2)

!Reduce data in case of parallelism
 call timab(48,1,tsec)
 call xmpi_sum(dos%fractions,comm_kptband,ierr)
 if (prtdosm>=1.or.fatbands_flag==1) then
   call xmpi_sum(dos%fractions_m,comm_kptband,ierr)
 end if
 if (paw_dos_flag==1) then
   call xmpi_sum(dos%fractions_paw1,comm_kptband,ierr)
   call xmpi_sum(dos%fractions_pawt1,comm_kptband,ierr)
 end if
 call timab(48,2,tsec)
 if (mpi_enreg%paral_spinor==1) then
   call xmpi_sum(dos%fractions,mpi_enreg%comm_spinor,ierr)
   if (prtdosm>=1.or.fatbands_flag==1) then
     call xmpi_sum(dos%fractions_m,mpi_enreg%comm_spinor,ierr)
   end if
   if (paw_dos_flag==1) then
     call xmpi_sum(dos%fractions_paw1, mpi_enreg%comm_spinor,ierr)
     call xmpi_sum(dos%fractions_pawt1, mpi_enreg%comm_spinor,ierr)
   end if
 end if

!Averaging: A quick hack for m-decomposed LDOS:
!BA: not valid in presence of spin-orbit coupling  !
 if (prtdosm==1.and.fatbands_flag==0) then
!  if pawfatbnd is activated, one think in the cubic harmonics basis
!  whereas prtdosm=1 is in the spherical harmonics basis.
!  the following trick is done in order to have everything
!  in the complex spherical basis (not useful for pawfatbnd if we want to
!  have e.g t2g and eg d-orbitals).
   do iat=1,dtset%natsph
     do il = 0, mbesslang-1
       do im = 1, il
         dos%fractions_m(:,:,:,mbesslang**2*(iat-1)+il**2+il+1+im) = &
&         (dos%fractions_m(:,:,:,mbesslang**2*(iat-1)+il**2+il+1+im) + &
&         dos%fractions_m(:,:,:,mbesslang**2*(iat-1)+il**2+il+1-im))/2
         dos%fractions_m(:,:,:,mbesslang**2*(iat-1)+il**2+il+1-im) = &
&         dos%fractions_m(:,:,:,mbesslang**2*(iat-1)+il**2+il+1+im)
       end do
     end do
   end do !iatom
 end if

 DBG_EXIT("COLL")

end subroutine partial_dos_fractions_paw
!!***

end module m_epjdos
!!***
