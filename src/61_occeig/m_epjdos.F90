!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_epjdos
!! NAME
!!  m_epjdos
!!
!! FUNCTION
!!  Tools for the computiation of electronic PJDOSes
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2016 ABINIT group (MVer, XG, SM, MT)
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
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_tetrahedron
 use m_splines
 use m_cgtools
 use m_atomdata

 use m_io_tools,     only : open_file
 use m_fstrings,     only : int2char4

 implicit none

 private
!!***

 public :: tetrahedron         ! Calculate DOS by tetrahedron method.
 public :: gaus_dos            ! Calculate DOS by gauss method.
 public :: dos_degeneratewfs   ! Average the contribution to a generalised DOS.
 public :: recip_ylm           ! Project input wavefunctions (real space) on to Ylm
 !public :: dens_in_sph         ! Calculate integrated density in sphere around each atom
!!***

!----------------------------------------------------------------------

contains  !============================================================
!!***

!!****f* m_pjedos/tetrahedron
!! NAME
!! tetrahedron
!!
!! FUNCTION
!! calculate DOS by tetrahedron method
!!
!! INPUTS
!!  dos_fractions= projections of wavefunctions on each angular momentum Ylm
!!     which is the weight going into the DOS for an l-decomposed dos
!!  dos_fractions_m= same as dos_fractions, but m-decomposed not just l-
!!  dos_fractions_paw1= contribution to dos fractions from the PAW partial waves (phi)
!!  dos_fractions_pawt1= contribution to dos fractions from the PAW pseudo partial waves (phi_tild)
!!  dtset     structured datatype, from which one uses :
!!   kpt(3,nkpt)  =irreducible kpoints
!!   kptrlatt(3,3)=lattice vectors for full kpoint grid
!!   nband        =number of electronic bands for each kpoint
!!   nkpt         =number of irreducible kpoints
!!   nshiftk      =number of kpoint grid shifts
!!   nsym         =number of symmetries
!!   pawprtdos    =option to output the individual contributions to the partial DOS (0, 1 or 2)
!!   shiftk(3,nshiftk)=kpoint shifts
!!   symrel(3,3,nsym)=symmetry matrices in real space
!!   typat(ntypat)=types of atoms
!!   usepaw       =option for PAW
!!  eigen(mband*nkpt*nsppol)=eigenvalues at irred kpoints
!!  fermie=Fermi energy
!!  fildata=name of the DOS output file
!!  mbesslang=maximum angular momentum for Bessel function expansion
!!  prtdosm=option for the m-contributions to the partial DOS
!!  ndosfraction= number of types of DOS we are calculating, e.g. the number
!!    of l channels. Could be much more general, for other types of partial DOS
!!  paw_dos_flag= option for partial dos in PAW
!!  rprimd(3,3)  =real space unit cell vectors
!!
!! OUTPUT
!!  (no explicit output)
!!  note: could have routine return DOS for other purposes, and separate printing
!!  in another routine (MJV 8/2010)
!!
!! NOTES
!!  This routine should be called by master only
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      destroy_tetra,dos_hdr_write,get_dos_1band,get_dos_1band_m
!!      get_full_kgrid,get_tetra_weight,init_tetra,int2char4,matr3inv,metric
!!      wrtout
!!
!! SOURCE

subroutine tetrahedron (dos_fractions,dos_fractions_m,dos_fractions_paw1,dos_fractions_pawt1,&
& dtset,fermie,eigen,fildata,mbesslang,prtdosm,ndosfraction,paw_dos_flag,rprimd)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'tetrahedron'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_56_recipspace
 use interfaces_61_occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: prtdosm,mbesslang,ndosfraction,paw_dos_flag
 real(dp),intent(in) :: fermie
 character(len=*),intent(in) :: fildata
 type(dataset_type),intent(in) :: dtset
!arrays
 real(dp),intent(in) :: dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction)
 real(dp),intent(in) :: dos_fractions_m(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*mbesslang*&
& max(min(prtdosm,1),0))
 real(dp),intent(in) :: dos_fractions_paw1(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*paw_dos_flag)
 real(dp),intent(in) :: dos_fractions_pawt1(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*paw_dos_flag)
 real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol),rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer,parameter :: bcorr0=0
 integer :: iat,iband,iene,ifract,ikpt,ioffset_eig,isppol,natsph,natsph_extra
 integer :: nene,nkpt_fullbz,prtdos,unitdos,ierr
 real(dp) :: buffer,deltaene,enemax,enemin,enex,integral_DOS,max_occ,rcvol
 real(dp) :: ucvol
 logical :: bigDOS
 character(len=10) :: tag
 character(len=500) :: frmt,frmt_extra,message
 character(len=fnlen) :: tmpfil
 character(len=80) :: errstr
 type(t_tetrahedron) :: tetrahedra
!arrays
 integer,allocatable :: indkpt(:),unitdos_arr(:)
 real(dp) :: gmet(3,3),gprimd(3,3),klatt(3,3),rlatt(3,3),rmet(3,3)
 real(dp),allocatable :: dtweightde(:,:),integ_dos(:,:),integ_dos_m(:,:)
 real(dp),allocatable :: kpt_fullbz(:,:),partial_dos(:,:)
 real(dp),allocatable :: partial_dos_m(:,:),tmp_eigen(:),total_dos(:,:)
 real(dp),allocatable :: total_dos_m(:,:),total_dos_paw1(:,:)
 real(dp),allocatable :: total_dos_pawt1(:,:),total_integ_dos(:,:)
 real(dp),allocatable :: total_integ_dos_m(:,:),tweight(:,:)
 real(dp),allocatable :: work_ndos(:,:),work_ndosmbessl(:,:)

! *********************************************************************

!m-decomposed DOS not compatible with PAW-decomposed DOS
 if(prtdosm>=1.and.paw_dos_flag==1) then
   message = ' m-decomposed DOS (prtdosm>=1) not&
&   compatible with PAW-decomposed DOS (pawprtdos=1) !'
   MSG_ERROR(message)
 end if

!Refuse only 1 kpoint: the algorithms are no longer valid. DOH !
 if (dtset%nkpt == 1) then
   MSG_WARNING('You need at least 2 kpoints to use the tetrahedron method. DOS cannot be computed')
   return
 end if

!Do not support nshiftk > 1: lattice must be decomposed into boxes
!and this is not always possible (I think) with bizzare shiftks
!normally at this point we have incorporated everything into
!kptrlatt, and only 1 shift is needed (in particular for MP grids).

 if (dtset%nshiftk > 1) then
   write(std_out,*) 'tetrahedron: skip subroutine.'
   write(std_out,*) 'Problem with a composite k-point grid.'
   write(std_out,*) 'Only simple lattices are supported. Action: use nshiftk=1.'
   write(std_out,*) 'dtset%nshiftk, dtset%shiftk = ', dtset%nshiftk, dtset%shiftk
   write(std_out,*) 'dtset%kptrlatt= ', dtset%kptrlatt
   MSG_WARNING('tetrahedron: skip subroutine. See message above')
   return
 end if

!Refuse nband different for different kpoints

 do isppol=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
     if ( dtset%nband(dtset%nkpt*(isppol-1) + ikpt) /= dtset%nband(1) ) then
       write(std_out,*) 'tetrahedron: skip subroutine.'
       write(std_out,*) 'nband must be the same for all kpoints'
       write(std_out,*) 'nband=', dtset%nband
       MSG_WARNING('tetrahedron: skip subroutine. See message above')
       return
     end if
   end do
 end do

! Calculate nkpt_fullbz
 nkpt_fullbz= dtset%kptrlatt(1,1)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,3) &
& +dtset%kptrlatt(1,2)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,1) &
& +dtset%kptrlatt(1,3)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,2) &
& -dtset%kptrlatt(1,2)*dtset%kptrlatt(2,1)*dtset%kptrlatt(3,3) &
& -dtset%kptrlatt(1,3)*dtset%kptrlatt(2,2)*dtset%kptrlatt(3,1) &
& -dtset%kptrlatt(1,1)*dtset%kptrlatt(2,3)*dtset%kptrlatt(3,2)
 nkpt_fullbz = nkpt_fullbz*dtset%nshiftk

 if (nkpt_fullbz==0) then
   write(std_out,*)'tetrahedron: skip subroutine.'
   write(std_out,*)'no homogeneous grid  of k-points is defined ...'
   write(std_out,*)'in order to obtain the DOS using the tetrahedron method,'
   write(std_out,*)'you need to re-define ngkpt or kptrlatt.'
   MSG_WARNING('tetrahedron: skip subroutine. See message above')
   return
 end if

!Make klatt
 rlatt(:,:) = dtset%kptrlatt(:,:)
 call matr3inv(rlatt,klatt)

!Get metric tensors
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 rcvol = abs (gprimd(1,1)*(gprimd(2,2)*gprimd(3,3)-gprimd(3,2)*gprimd(2,3)) &
& -gprimd(2,1)*(gprimd(1,2)*gprimd(3,3)-gprimd(3,2)*gprimd(1,3)) &
& +gprimd(3,1)*(gprimd(1,2)*gprimd(2,3)-gprimd(2,2)*gprimd(1,3)))

 ABI_ALLOCATE(indkpt,(nkpt_fullbz))
 ABI_ALLOCATE(kpt_fullbz,(3,nkpt_fullbz))

!Make full kpoint grid and get equivalence to irred kpoints
 call get_full_kgrid(indkpt,dtset%kpt,kpt_fullbz,dtset%kptrlatt,&
& dtset%nkpt,nkpt_fullbz,dtset%nshiftk,dtset%nsym,dtset%shiftk,&
& dtset%symrel)

!Get tetrahedra, ie indexes of the full kpoints at their summits
 call init_tetra (indkpt,gprimd,klatt,kpt_fullbz,nkpt_fullbz,tetrahedra, ierr, errstr)
 ABI_CHECK(ierr==0,errstr)

 natsph=dtset%natsph
 natsph_extra=dtset%natsph_extra

!Open the DOS file
 if (dtset%prtdos == 2 .or. dtset%prtdos == 5) then
   if (open_file(fildata,message,newunit=unitdos,status='unknown',form='formatted') /= 0) then
     MSG_ERROR(message)
   end if
   rewind(unitdos)

 else if (dtset%prtdos == 3) then
   ABI_ALLOCATE(unitdos_arr,(natsph+natsph_extra))
   do iat=1,natsph
     call int2char4(dtset%iatsph(iat),tag)
     ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
     tmpfil = trim(fildata)//'_AT'//trim(tag)
     if (open_file(tmpfil, message, newunit=unitdos_arr(iat), status='unknown',form='formatted') /= 0) then
       MSG_ERROR(message)
     end if 
     rewind(unitdos_arr(iat))
     write(std_out,*) 'opened file : ', trim(tmpfil), '  unit', unitdos_arr(iat)
   end do
!  do extra spheres in vacuum too
   do iat=1,natsph_extra
!    FIXME: added random offset of 1000 to atom index in file name for extra spheres.
!    Someday we will have more than 1000 atoms...
     call int2char4(1000+iat,tag)
     ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
     tmpfil = trim(fildata)//'_AT'//trim(tag)
     if (open_file(tmpfil, message, newunit=unitdos_arr(natsph+iat), status='unknown',form='formatted') /= 0) then
       MSG_ERROR(message)
     end if
     rewind(unitdos_arr(natsph+iat))
     write(std_out,*) 'opened file : ', trim(tmpfil), '  unit', unitdos_arr(natsph+iat)
   end do
 end if

!Write the header of the DOS file, and determine the energy range and spacing
 prtdos=dtset%prtdos
 buffer=0.01_dp ! Size of the buffer around the min and max ranges
 if (dtset%prtdos == 2 .or. dtset%prtdos == 5) then
   call dos_hdr_write(buffer,deltaene,dtset%dosdeltae,eigen,enemax,enemin,fermie,dtset%mband,&
&   dtset%nband,nene,dtset%nkpt,dtset%nsppol,dtset%occopt,prtdos,&
&   dtset%tphysel,dtset%tsmear,unitdos)
 else if (dtset%prtdos == 3) then
   do iat=1,natsph+natsph_extra
     call dos_hdr_write(buffer,deltaene,dtset%dosdeltae,eigen,enemax,enemin,fermie,dtset%mband,&
&     dtset%nband,nene,dtset%nkpt,dtset%nsppol,dtset%occopt,prtdos,&
&     dtset%tphysel,dtset%tsmear,unitdos_arr(iat))
   end do
 end if

 ABI_ALLOCATE(partial_dos,(nene,ndosfraction))
 ABI_ALLOCATE(integ_dos,(nene,ndosfraction))
 ABI_ALLOCATE(total_dos,(nene,ndosfraction))
 ABI_ALLOCATE(total_integ_dos,(nene,ndosfraction))
 ABI_ALLOCATE(tweight,(dtset%nkpt,nene))
 ABI_ALLOCATE(dtweightde,(dtset%nkpt,nene))

 if (paw_dos_flag==1) then
   ABI_ALLOCATE(total_dos_paw1,(nene,ndosfraction))
   ABI_ALLOCATE(total_dos_pawt1,(nene,ndosfraction))
 end if
 if (prtdosm>=1) then
   ABI_ALLOCATE(partial_dos_m,(nene,ndosfraction*mbesslang))
   ABI_ALLOCATE(integ_dos_m,(nene,ndosfraction*mbesslang))
   ABI_ALLOCATE(total_dos_m,(nene,ndosfraction*mbesslang))
   ABI_ALLOCATE(total_integ_dos_m,(nene,ndosfraction*mbesslang))
 end if

!Get maximum occupation value (2 or 1)
 if (dtset%nspinor == 1 .and. dtset%nsppol == 1) then
   max_occ = two
 else
   max_occ = one
 end if


!-------------------------------------------------------------------
!For each spin polarisation and band, interpolate band over kpoints
!calculate integration weights and DOS contib from
!-------------------------------------------------------------------

 ! Workspace arrays.
 ABI_MALLOC(work_ndos, (dtset%nkpt, ndosfraction))
 ABI_MALLOC(work_ndosmbessl, (dtset%nkpt, ndosfraction*mbesslang))

 do isppol = 1, dtset%nsppol

   total_dos (:,:) = zero
   total_integ_dos (:,:) = zero
   if (prtdosm>=1) then
     total_dos_m (:,:) = zero
   end if
   if (paw_dos_flag==1) then
     total_dos_paw1(:,:)=zero;total_dos_pawt1(:,:)=zero
   end if

   if (dtset%nsppol==2) then
     if(isppol==1) write(message,'(a,16x,a)')  '#','Spin-up DOS'
     if(isppol==2) write(message,'(2a,16x,a)')  ch10,'#','Spin-dn DOS'
! NB: dtset%prtdos == 5 should not happen for nsppol==2
     if (dtset%prtdos == 2 .or. dtset%prtdos == 5) then
       call wrtout(unitdos,message,'COLL')
     else if (dtset%prtdos == 3) then
       do iat=1,natsph+natsph_extra
         call wrtout(unitdos_arr(iat),message,'COLL')
       end do
     end if
   end if

   ABI_ALLOCATE(tmp_eigen,(dtset%nkpt))
   do iband = 1, dtset%nband(1)

     ioffset_eig = dtset%mband*dtset%nkpt*(isppol-1)

     tmp_eigen(:) = zero
!    For each band get its contribution
     do ikpt=1,dtset%nkpt
       tmp_eigen(ikpt) = eigen(ioffset_eig + dtset%mband*(ikpt-1) + iband)
     end do
!
!    calculate general integration weights at each irred kpoint as in Blochl et al
!    PRB 49 16223
!
     call get_tetra_weight(tmp_eigen,enemin,enemax,&
&     max_occ,nene,dtset%nkpt,tetrahedra,bcorr0,&
&     tweight,dtweightde,xmpi_comm_self)

!
!    calculate DOS and integrated DOS projected with the input dos_fractions
!
     if (paw_dos_flag==1) then
       work_ndos = dos_fractions_paw1(:,iband,isppol,:)

       call get_dos_1band (work_ndos,enemin,enemax,integ_dos(:,:),nene,dtset%nkpt,ndosfraction,&
&       partial_dos(:,:),tweight,dtweightde)

       do ifract=1,ndosfraction
         do iene=1,nene
           total_dos_paw1(iene,ifract)=total_dos_paw1(iene,ifract)+partial_dos(iene,ifract)
         end do
       end do

       work_ndos = dos_fractions_pawt1(:,iband,isppol,:)
       call get_dos_1band (work_ndos,enemin,enemax,integ_dos(:,:),nene,dtset%nkpt,ndosfraction,&
&       partial_dos(:,:),tweight,dtweightde)

       do ifract=1,ndosfraction
         do iene=1,nene
           total_dos_pawt1(iene,ifract)=total_dos_pawt1(iene,ifract)+partial_dos(iene,ifract)
         end do
       end do
     end if

     work_ndos = dos_fractions(:,iband,isppol,:)
     call get_dos_1band (work_ndos,enemin,enemax,integ_dos(:,:),nene,dtset%nkpt,ndosfraction,&
&     partial_dos(:,:),tweight,dtweightde)

     if (prtdosm>=1) then
       work_ndosmbessl = dos_fractions_m(:,iband,isppol,:)
       call get_dos_1band_m (work_ndosmbessl,enemin,enemax,integ_dos_m(:,:),nene,dtset%nkpt,ndosfraction*mbesslang,&
&       partial_dos_m(:,:),tweight,dtweightde)
     end if

!    Add to total dos
     do ifract=1,ndosfraction
       do iene=1,nene
         total_dos(iene,ifract) = total_dos(iene,ifract) + partial_dos(iene,ifract)
         total_integ_dos(iene,ifract) = total_integ_dos(iene,ifract) + integ_dos(iene,ifract)
       end do
     end do

     if (prtdosm>=1) then
       do ifract=1,ndosfraction*mbesslang
         do iene=1,nene
           total_dos_m(iene,ifract) = total_dos_m(iene,ifract) + partial_dos_m(iene,ifract)
           total_integ_dos_m(iene,ifract) = total_integ_dos_m(iene,ifract) + integ_dos_m(iene,ifract)
         end do
       end do
     end if

   end do ! iband
   bigDOS=(maxval(total_dos)>999._dp)

   ABI_DEALLOCATE(tmp_eigen)

   call wrtout(std_out,'about to write to the DOS file ',"COLL")
! header lines depend on the type of DOS (projected etc...) which is output
   enex=enemin
   integral_DOS=zero
   if(prtdos==2)then
     write(message, '(a)' )'# energy(Ha)     DOS  integrated DOS'
     call wrtout(unitdos,message,'COLL')

   else if(prtdos==3)then
     if (paw_dos_flag/=1.or.dtset%pawprtdos==2) then
       do iat=1,natsph
         write(message, '(3a,i5,a,i5,a,a,es16.6,3a)' ) &
&         '# Local DOS (columns 2-6) and integrated local DOS (columns 7-11),',ch10,&
&         '# for atom number iat=',iat,'  iatom=',dtset%iatsph(iat),ch10,&
&         '# inside sphere of radius ratsph=',dtset%ratsph(dtset%typat(dtset%iatsph(iat))),' Bohr.',ch10,"#"
         call wrtout(unitdos_arr(iat),message,'COLL')
         if (dtset%usepaw==1.and.dtset%pawprtdos==2) then
           write(message, '(3a)' ) &
&           '# PAW: note that only all-electron on-site part has been used to compute DOS !',ch10,"#"
           call wrtout(unitdos_arr(iat),message,'COLL')
         end if
         if (bigDOS) then
           write(message, '(a,a)' ) &
&           '# energy(Ha)   l=0       l=1       l=2       l=3       l=4',&
&           '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
         else
           write(message, '(a,a)' ) &
&           '# energy(Ha)  l=0      l=1      l=2      l=3      l=4',&
&           '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
         end if
         if (prtdosm>=1) then
           write(message, '(7a)' ) trim(message),'          ',&
&           '  lm=0 0',&
&           '  lm=1-1  lm=1 0  lm=1 1',&
&           '  lm=2-2  lm=2-1  lm=2 0  lm=2 1  lm=2 2',&
&           '  lm=3-3  lm=3-2  lm=3-1  lm=3 0  lm=3 1  lm=3 2  lm=3 3',&
&           '  lm=4-4  lm=4-3  lm=4-2  lm=4-1  lm=4 0  lm=4 1  lm=4 2  lm=4 3  lm=4 4'
         end if
         call wrtout(unitdos_arr(iat),message,'COLL')
       end do
     else
       do iat=1,natsph
         write(message, '(9a,i5,a,i5,a,a,es16.6,3a)' ) &
&         '# Local DOS (columns 2-6),',ch10,&
&         '#  plane-waves contrib. to DOS (columns 7-11),',ch10,&
&         '#  AE on-site  contrib. to DOS (columns 12-16),',ch10,&
&         '# -PS on-site  contrib. to DOS (columns 17-21),',ch10,&
&         '# for atom number iat=',iat,'  iatom=',dtset%iatsph(iat),ch10,&
&         '# inside sphere of radius ratsph=',dtset%ratsph(dtset%typat(dtset%iatsph(iat))),' Bohr.',ch10,"#"
         call wrtout(unitdos_arr(iat),message,'COLL')
         if (bigDOS) then
           write(message, '(4a)' ) &
&           '#energy(Ha)   l=0       l=1       l=2       l=3       l=4',&
&           '       (PW)  l=0       l=1       l=2       l=3       l=4',&
&           '      (Phi)  l=0       l=1       l=2       l=3       l=4',&
&           '     (tPhi)  l=0       l=1       l=2       l=3       l=4'
         else
           write(message, '(4a)' ) &
&           '#energy(Ha)  l=0      l=1      l=2      l=3      l=4',&
&           '       (PW) l=0      l=1      l=2      l=3      l=4',&
&           '      (Phi) l=0      l=1      l=2      l=3      l=4',&
&           '     (tPhi) l=0      l=1      l=2      l=3      l=4'
         end if
         call wrtout(unitdos_arr(iat),message,'COLL')
       end do
     end if
     do iat=1,natsph_extra
       write(message, '(3a,i5,2a,es16.6,3a)' ) &
&       '# Local DOS (columns 2-6) and integrated local DOS (columns 7-11),',ch10,&
&       '# for non-atomic sphere number iat=',iat,ch10,&
&       '# of radius ratsph=',dtset%ratsph_extra,' Bohr.',ch10,"#"
       call wrtout(unitdos_arr(natsph+iat),message,'COLL')
       if (bigDOS) then
         write(message, '(a,a)' ) &
&         '# energy(Ha)   l=0       l=1       l=2       l=3       l=4',&
&         '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
       else
         write(message, '(a,a)' ) &
&         '# energy(Ha)  l=0      l=1      l=2      l=3      l=4',&
&         '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
       end if
       if (prtdosm>=1) then
         write(message, '(7a)' ) trim(message),'          ',&
&         '  lm=0 0',&
&         '  lm=1-1  lm=1 0  lm=1 1',&
&         '  lm=2-2  lm=2-1  lm=2 0  lm=2 1  lm=2 2',&
&         '  lm=3-3  lm=3-2  lm=3-1  lm=3 0  lm=3 1  lm=3 2  lm=3 3',&
&         '  lm=4-4  lm=4-3  lm=4-2  lm=4-1  lm=4 0  lm=4 1  lm=4 2  lm=4 3  lm=4 4'
       end if
       call wrtout(unitdos_arr(natsph+iat),message,'COLL')
     end do
   else if(prtdos==5)then
     write(message, '(a)' )'# energy(Ha)     DOS up,up  up,dn  dn,up  up,up  sigma_x sigma_y sigma_z  and integrated DOS components'
     call wrtout(unitdos,message,'COLL')

   end if ! prtdos value
   !
   !  finished with header printing
   !


!  Write the DOS value in the DOS file
   if(prtdos==2)then
     do iene=1,nene
!      Print the data for this energy. Note the upper limit, to be
!      consistent with the format. The use of "E" format is not adequate,
!      for portability of the self-testing procedure.
!      write(message, '(i5,f9.4,f14.6)' ) iene-1,enex,total_dos(iene,:)
       write(message, '(f11.5,5f10.4,10x,5f8.2)') enex, &
&       min(total_dos(iene,:),9999.9999_dp),&
&       total_integ_dos(iene,:)
       call wrtout(unitdos,message,'COLL')
       enex=enex+deltaene
     end do
   else if(prtdos==3)then
     if (paw_dos_flag/=1.or.dtset%pawprtdos==2) then
       if (bigDOS) then
         frmt='(f11.5,5f10.4,10x,5f8.2,10x,25f8.2)'
       else
         frmt='(f11.5,5f9.4 ,10x,5f8.2,10x,25f8.2)'
       end if
!      for extra atoms in vacuum need more precision
       frmt_extra='(f11.5,5f20.16,10x,5f20.16,10x,25f20.16)'
       do iene=1,nene
         do iat=1,natsph
!          Note the upper limit, to be
!          consistent with the format. The use of "E" format is not adequate,
!          for portability of the self-testing procedure.
           if (prtdosm==0) then
             write(message,fmt=frmt) enex, &
&             min(total_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang),9999.9999_dp), &
&             total_integ_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang)
           else
             write(message,fmt=frmt) enex, &
&             min(total_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang),9999.9999_dp),&
&             total_integ_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang),&
&             min(total_dos_m(iene,(iat-1)*mbesslang**2+1:iat*mbesslang**2),9999.9999_dp)
           end if

           call wrtout(unitdos_arr(iat),message,'COLL')
         end do
         do iat=natsph+1,natsph+natsph_extra
!          Note the upper limit, to be
!          consistent with the format. The use of "E" format is not adequate,
!          for portability of the self-testing procedure.
           if (prtdosm==0) then
             write(message,fmt=frmt_extra) enex, &
&             total_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang), &
&             total_integ_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang)
           else
             write(message,fmt=frmt_extra) enex, &
&             total_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang),&
&             total_integ_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang),&
&             total_dos_m(iene,(iat-1)*mbesslang**2+1:iat*mbesslang**2)
           end if

           call wrtout(unitdos_arr(iat),message,'COLL')
         end do
         enex=enex+deltaene
       end do
     else
       if (bigDOS) then
         frmt='(f11.5,5f10.4,3(6x,5f10.4))'
       else
         frmt='(f11.5,5f9.4 ,3(6x,5f9.4 ))'
       end if
!      for extra atom spheres in vacuum need more precision
       frmt_extra='(f11.5,5f20.16,3(6x,5f20.16))'
       do iene=1,nene
         do iat=1,natsph
           write(message,fmt=frmt) enex, &
&           min(total_dos(iene,iat*5-4:iat*5),9999.9999_dp),&
&           min(total_dos(iene,iat*5-4:iat*5)-total_dos_paw1(iene,iat*5-4:iat*5)&
&           +total_dos_pawt1(iene,iat*5-4:iat*5),9999.9999_dp),&
&           min(total_dos_paw1(iene,iat*5-4:iat*5),9999.9999_dp),&
&           min(total_dos_pawt1(iene,iat*5-4:iat*5),9999.9999_dp)
           call wrtout(unitdos_arr(iat),message,'COLL')
         end do
         do iat=natsph+1,natsph+natsph_extra
           write(message,fmt=frmt_extra) enex, &
&           min(total_dos(iene,iat*5-4:iat*5),9999.9999_dp),&
&           min(total_dos(iene,iat*5-4:iat*5)-total_dos_paw1(iene,iat*5-4:iat*5)&
&           +total_dos_pawt1(iene,iat*5-4:iat*5),9999.9999_dp),&
&           min(total_dos_paw1(iene,iat*5-4:iat*5),9999.9999_dp),&
&           min(total_dos_pawt1(iene,iat*5-4:iat*5),9999.9999_dp)
           call wrtout(unitdos_arr(iat),message,'COLL')
         end do
         enex=enex+deltaene
       end do
     end if

   else if(prtdos==5)then
     if (bigDOS) then
       frmt='(f11.5,7f10.4,10x,7f8.2)'
     else
       frmt='(f11.5,7f9.4 ,10x,7f8.2)'
     end if
     do iene=1,nene
       write(message,fmt=frmt) enex, &
&       min(total_dos(iene,1:7),9999.9999_dp),&
&       total_integ_dos(iene,1:7)
       call wrtout(unitdos,message,'COLL')
       enex=enex+deltaene
     end do

   end if

!  integral_DOS=integral_DOS+deltaene*sum(total_dos(iene,:))
   integral_DOS=sum(total_integ_dos(nene,:))
   write(message, '(a,es16.8)' ) ' tetrahedron : integrate to',integral_DOS
   call wrtout(std_out,message,'COLL')

 end do ! isppol

 ABI_FREE(work_ndos)
 ABI_FREE(work_ndosmbessl)

 if(prtdos==2 .or. prtdos==5) then
   close(unitdos)
 else if(prtdos==3) then
   do iat=1,natsph+natsph_extra
     close(unitdos_arr(iat))
   end do
   ABI_DEALLOCATE(unitdos_arr)
 end if

 ABI_DEALLOCATE(indkpt)
 ABI_DEALLOCATE(kpt_fullbz)
 ABI_DEALLOCATE(partial_dos)
 ABI_DEALLOCATE(integ_dos)
 ABI_DEALLOCATE(total_dos)
 ABI_DEALLOCATE(total_integ_dos)
 ABI_DEALLOCATE(tweight)
 ABI_DEALLOCATE(dtweightde)

 if (prtdosm>=1)  then
   ABI_DEALLOCATE(partial_dos_m)
   ABI_DEALLOCATE(integ_dos_m)
   ABI_DEALLOCATE(total_dos_m)
   ABI_DEALLOCATE(total_integ_dos_m)
 end if

 if (paw_dos_flag==1)  then
   ABI_DEALLOCATE(total_dos_paw1)
   ABI_DEALLOCATE(total_dos_pawt1)
 end if

 call destroy_tetra(tetrahedra)

end subroutine tetrahedron
!!***

!----------------------------------------------------------------------

!!****f* m_pjedos/gaus_dos
!! NAME
!! gaus_dos
!!
!! FUNCTION
!! Calculate DOS by gauss method
!!
!! INPUTS
!!  dos_fractions= projections of wavefunctions on each angular momentum Ylm
!!     which is the weight going into the DOS for an l-decomposed dos
!!  dos_fractions_paw1= contribution to dos fractions from the PAW partial waves (phi)
!!  dos_fractions_pawt1= contribution to dos fractions from the PAW pseudo partial waves (phi_tild)
!!  dtset     structured datatype, from which one uses :
!!   kpt(3,nkpt)  =irreducible kpoints
!!   kptrlatt(3,3)=lattice vectors for full kpoint grid
!!   nband        =number of electronic bands for each kpoint
!!   nkpt         =number of irreducible kpoints
!!   nshiftk      =number of kpoint grid shifts
!!   nsym         =number of symmetries
!!   pawprtdos    =option to output the individual contributions to the partial DOS (0, 1 or 2)
!!   shiftk(3,nshiftk)=kpoint shifts
!!   symrel(3,3,nsym)=symmetry matrices in real space
!!   typat(ntypat)=types of atoms
!!   usepaw       =option for PAW
!!  eigen(mband*nkpt*nsppol)=eigenvalues at irred kpoints
!!  fermie=Fermi energy
!!  fildata=name of the DOS output file
!!  mbesslang=maximum angular momentum for Bessel function expansion
!!  m_dos_flag=option for the m-contributions to the partial DOS
!!  ndosfraction= number of types of DOS we are calculating, e.g. the number
!!    of l channels. Could be much more general, for other types of partial DOS
!!  paw_dos_flag= option for partial dos in PAW
!!
!! OUTPUT
!!  (no explicit output)
!!
!! NOTES
!!  This routine should be called by master only
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      dos_hdr_write,int2char4,splfit,wrtout
!!
!! SOURCE

subroutine gaus_dos(dos_fractions,&
& dos_fractions_paw1,dos_fractions_pawt1,dtset,fermie,eigen,&
& fildata,mbesslang,m_dos_flag,ndosfraction,paw_dos_flag) 


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gaus_dos'
 use interfaces_14_hidewrite
 use interfaces_61_occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mbesslang,m_dos_flag,ndosfraction,paw_dos_flag
 real(dp),intent(in) :: fermie
 character(len=fnlen),intent(in) :: fildata
 type(dataset_type),intent(in) :: dtset
!arrays
 real(dp),intent(in) :: dos_fractions(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction)
 real(dp),intent(in) :: dos_fractions_paw1(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*paw_dos_flag)
 real(dp),intent(in) :: dos_fractions_pawt1(dtset%nkpt,dtset%mband,dtset%nsppol,ndosfraction*paw_dos_flag)
 real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: nptsdiv2=6000
 integer :: bantot,ii,iat,iband,iene,ifract,ikpt,index,isppol,natsph,nkpt,nsppol,nene,prtdos
 real(dp) :: buffer,deltaene,enemax,enemin,integral_DOS,max_occ,enex,tsmearinv,tsmear,limit_occ,tratio
 real(dp) :: dblsmr,dsqrpi,increm,xx
 logical :: bigDOS
 character(len=10) :: tag
 character(len=500) :: frmt,message
 character(len=fnlen) :: tmpfil
!scalars
 integer,allocatable :: unitdos_arr(:)
 real(dp) :: xgrid(-nptsdiv2:nptsdiv2),smdfun(-nptsdiv2:nptsdiv2,2)
 real(dp),allocatable :: dos(:),arg(:),derfun(:),partial_dos(:,:,:),integ_dos(:,:,:),total_dos(:,:)
 real(dp),allocatable :: total_integ_dos(:,:),total_dos_paw1(:,:),total_dos_pawt1(:,:)
 real(dp),allocatable :: integ_dos_m(:,:,:),partial_dos_m(:,:,:)
 real(dp),allocatable :: total_dos_m(:,:),total_integ_dos_m(:,:)

! *********************************************************************

!Open the DOS file

 natsph=dtset%natsph
 nsppol=dtset%nsppol
 nkpt=dtset%nkpt
 prtdos=dtset%prtdos
 deltaene=dtset%dosdeltae
 tsmear=dtset%tsmear
 tsmearinv=one/tsmear

 ABI_ALLOCATE(unitdos_arr,(natsph))
 do iat=1,natsph
   call int2char4(dtset%iatsph(iat),tag)
   ABI_CHECK((tag(1:1)/='#'),'Bug: string length too short!')
   tmpfil = trim(fildata)//'_AT'//trim(tag)
   if (open_file(tmpfil, message, newunit=unitdos_arr(iat), status='unknown',form='formatted') /= 0) then
     MSG_ERROR(message)
   end if
   rewind(unitdos_arr(iat))
   write(std_out,*) 'opened file : ', tmpfil, 'unit', unitdos_arr(iat)
 end do


!Write the header of the DOS file, and determine the energy range and spacing

 buffer=0.01_dp ! Size of the buffer around the min and max ranges
 
 do iat=1,natsph
   call dos_hdr_write(buffer,deltaene,dtset%dosdeltae,eigen,enemax,enemin,fermie,dtset%mband,&
&   dtset%nband,nene,dtset%nkpt,dtset%nsppol,dtset%occopt,prtdos,&
&   dtset%tphysel,dtset%tsmear,unitdos_arr(iat))
 end do

 bantot=sum(dtset%nband(:))
 ABI_ALLOCATE(arg,(bantot))
 ABI_ALLOCATE(derfun,(bantot))
 ABI_ALLOCATE(dos,(bantot))
 
!DEBUG
!write(std_out,*) ' ndosfraction,dtset%mband,dtset%nkpt,nene',&
!&              ndosfraction,dtset%mband,dtset%nkpt,nene
!ENDDEBUG
 ABI_ALLOCATE(partial_dos,(nene,ndosfraction,dtset%mband))
 ABI_ALLOCATE(integ_dos,(nene,ndosfraction,dtset%mband))
 ABI_ALLOCATE(total_dos,(nene,ndosfraction))
 ABI_ALLOCATE(total_integ_dos,(nene,ndosfraction))

 if (paw_dos_flag==1) then
   ABI_ALLOCATE(total_dos_paw1,(nene,ndosfraction))
   ABI_ALLOCATE(total_dos_pawt1,(nene,ndosfraction))
 end if
 if (m_dos_flag>=1) then
   ABI_ALLOCATE(partial_dos_m,(nene,ndosfraction*mbesslang,dtset%mband))
   ABI_ALLOCATE(integ_dos_m,(nene,ndosfraction*mbesslang,dtset%mband))
   ABI_ALLOCATE(total_dos_m,(nene,ndosfraction*mbesslang))
   ABI_ALLOCATE(total_integ_dos_m,(nene,ndosfraction*mbesslang))
 end if
 
!Get maximum occupation value (2 or 1)
 if (dtset%nspinor == 1 .and. dtset%nsppol == 1) then
   max_occ = two
 else
   max_occ = one
 end if

!Define xgrid and gaussian on smdfun following getnel
 limit_occ=6.0_dp
 dblsmr=0
 if (abs(dtset%tphysel)>tol12) then
   if (abs(dtset%tsmear)>tol12) then
     dblsmr = 1
   end if
 end if

 if(dtset%occopt==3)limit_occ=30.0_dp
 if(dblsmr /= 0) then
   tratio = dtset%tsmear / dtset%tphysel
   limit_occ=30.0_dp + 6.0_dp*tratio
 end if

 increm=limit_occ/nptsdiv2
 do ii=-nptsdiv2,nptsdiv2
   xgrid(ii)=ii*increm
 end do
 dsqrpi=1.0_dp/sqrt(pi)
 do ii=0,nptsdiv2
   xx=xgrid(ii)
   smdfun( ii,1)=dsqrpi*exp(-xx**2)
   smdfun(-ii,1)=smdfun(ii,1)
 end do

!calculate DOS and integrated DOS projected with the input dos_fractions
!
 total_dos_paw1=0
 total_dos_pawt1=0
 total_dos=0
 enex=enemin
 do iene=1,nene
   arg(:)=(enex-eigen(1:bantot))*tsmearinv
   call splfit(xgrid,derfun,smdfun,0,arg,dos,(2*nptsdiv2+1),bantot)
   index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt
       do iband=1,dtset%nband(ikpt)
         index=index+1        
         do ifract=1,ndosfraction
           total_dos_paw1(iene,ifract)=total_dos_paw1(iene,ifract)+&
&           dos_fractions_paw1(ikpt,iband,isppol,ifract)*dtset%wtk(ikpt)* &
&           max_occ*dos(index)*tsmearinv
           total_dos_pawt1(iene,ifract)=total_dos_pawt1(iene,ifract)+&
&           dos_fractions_pawt1(ikpt,iband,isppol,ifract)*dtset%wtk(ikpt)* &
&           max_occ*dos(index)*tsmearinv
           total_dos(iene,ifract) = total_dos(iene,ifract) + &
&           dos_fractions(ikpt,iband,isppol,ifract)*dtset%wtk(ikpt)*max_occ*& 
&           dos(index)*tsmearinv   
!          write(99,*) iene,total_dos_paw1(iene,ifract),total_dos_pawt1(iene,ifract),total_dos(iene,ifract)
         end do
       end do ! iband
     end do ! ikpt
   end do ! isppol
!  write(99,*) iene,total_dos_paw1(iene,1),total_dos_pawt1(iene,1),total_dos(iene,1)
   enex=enex+deltaene
 end do   ! iene
 

 write(std_out,*) 'about to write to the DOS file '
!Write the DOS value in the DOS file
 do isppol=1,nsppol
   enex=enemin
   integral_DOS=zero
   
!  bigDOS=(maxval(total_dos)>999._dp)
   if (paw_dos_flag/=1.or.dtset%pawprtdos==2) then
     do iat=1,natsph
       write(message, '(3a,i5,a,i5,a,a,es16.6,3a)' ) &
&       '# Local DOS (columns 3-7) and integrated local DOS (columns 8-12),',ch10,&
&       '# for atom number iat=',iat,'  iatom=',dtset%iatsph(iat),ch10,&
&       '# inside sphere of radius ratsph=',dtset%ratsph(dtset%typat(dtset%iatsph(iat))),' Bohr.',ch10,"#"
       call wrtout(unitdos_arr(iat),message,'COLL')
       if (dtset%usepaw==1.and.dtset%pawprtdos==2) then
         write(message, '(3a)' ) &
&         '# PAW: note that only all-electron on-site part has been used to compute DOS !',ch10,"#"
         call wrtout(unitdos_arr(iat),message,'COLL')
       end if
       if (bigDOS) then
         write(message, '(a,a)' ) &
&         '# energy(Ha)   l=0       l=1       l=2       l=3       l=4',&
&         '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
       else
         write(message, '(a,a)' ) &
&         '# energy(Ha)  l=0      l=1      l=2      l=3      l=4',&
&         '    (integral=>)  l=0     l=1     l=2     l=3     l=4'
       end if
       if (m_dos_flag>=1) then
         write(message, '(7a)' ) trim(message),'          ',&
&         '  lm=0 0',&
&         '  lm=1-1  lm=1 0  lm=1 1',&
&         '  lm=2-2  lm=2-1  lm=2 0  lm=2 1  lm=2 2',&
&         '  lm=3-3  lm=3-2  lm=3-1  lm=3 0  lm=3 1  lm=3 2  lm=3 3',&
&         '  lm=4-4  lm=4-3  lm=4-2  lm=4-1  lm=4 0  lm=4 1  lm=4 2  lm=4 3  lm=4 4'
       end if
       call wrtout(unitdos_arr(iat),message,'COLL')
     end do
   else
     do iat=1,natsph
       write(message, '(9a,i5,a,i5,a,a,es16.6,3a)' ) &
&       '# Local DOS (columns 3-7),',ch10,&
&       '#  plane-waves contrib. to DOS (columns 8-12),',ch10,&
&       '#  AE on-site  contrib. to DOS (columns 13-17),',ch10,&
&       '# -PS on-site  contrib. to DOS (columns 18-22),',ch10,&
&       '# for atom number iat=',iat,'  iatom=',dtset%iatsph(iat),ch10,&
&       '# inside sphere of radius ratsph=',dtset%ratsph(dtset%typat(dtset%iatsph(iat))),' Bohr.',ch10,"#"
       call wrtout(unitdos_arr(iat),message,'COLL')
       if (bigDOS) then
         write(message, '(4a)' ) &
&         '#energy(Ha)   l=0       l=1       l=2       l=3       l=4',&
&         '       (PW)  l=0       l=1       l=2       l=3       l=4',&
&         '      (Phi)  l=0       l=1       l=2       l=3       l=4',&
&         '     (tPhi)  l=0       l=1       l=2       l=3       l=4'
       else
         write(message, '(4a)' ) &
&         '#energy(Ha)  l=0      l=1      l=2      l=3      l=4',&
&         '       (PW) l=0      l=1      l=2      l=3      l=4',&
&         '      (Phi) l=0      l=1      l=2      l=3      l=4',&
&         '     (tPhi) l=0      l=1      l=2      l=3      l=4'
       end if
       call wrtout(unitdos_arr(iat),message,'COLL')
     end do
   end if

   if (paw_dos_flag/=1.or.dtset%pawprtdos==2) then
     if (bigDOS) then
       frmt='(f11.5,5f10.4,10x,5f8.2,10x,25f8.2)'
     else
       frmt='(f11.5,5f9.4 ,10x,5f8.2,10x,25f8.2)'
     end if
     do iene=1,nene
       do iat=1,natsph
!        DEBUG
!        write(message, '(i6,5f9.4,10x,5f7.2))') iene-1, &
!        &     total_dos(iene,iat*5-2)
!        ENDDEBUG
!        Note the upper limit, to be
!        consistent with the format. The use of "E" format is not adequate,
!        for portability of the self-testing procedure.
         if (m_dos_flag==0) then
           write(message,fmt=frmt) enex, &
&           min(total_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang),9999.9999_dp), &
&           total_integ_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang)
         else
           write(message,fmt=frmt) enex, &
&           min(total_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang),9999.9999_dp),&
&           total_integ_dos(iene,(iat-1)*mbesslang+1:iat*mbesslang),&
&           min(total_dos_m(iene,(iat-1)*mbesslang**2+1:iat*mbesslang**2),9999.9999_dp)
         end if

         call wrtout(unitdos_arr(iat),message,'COLL')
       end do
       enex=enex+deltaene
     end do
   else
     if (bigDOS) then
       frmt='(f11.5,5f10.4,3(6x,5f10.4))'
     else
       frmt='(f11.5,5f9.4 ,3(6x,5f9.4 ))'
     end if
     do iene=1,nene
       do iat=1,natsph
         write(message,fmt=frmt) enex, &
&         min(total_dos(iene,iat*5-4:iat*5),9999.9999_dp),&
&         min(total_dos(iene,iat*5-4:iat*5)-total_dos_paw1(iene,iat*5-4:iat*5)&
&         +total_dos_pawt1(iene,iat*5-4:iat*5),9999.9999_dp),&
&         min(total_dos_paw1(iene,iat*5-4:iat*5),9999.9999_dp),&
&         min(total_dos_pawt1(iene,iat*5-4:iat*5),9999.9999_dp)
         call wrtout(unitdos_arr(iat),message,'COLL')
       end do
       enex=enex+deltaene
     end do
   end if
   

!  integral_DOS=integral_DOS+deltaene*sum(total_dos(iene,:))
   integral_DOS=sum(total_integ_dos(nene,:))
   write(message, '(a,es16.8)' ) ' tetrahedron : integrate to',integral_DOS
   call wrtout(std_out,message,'COLL')

 end do ! isppol

 do iat=1,natsph
   close(unitdos_arr(iat))
 end do

 ABI_DEALLOCATE(unitdos_arr)

 ABI_DEALLOCATE(partial_dos)
 ABI_DEALLOCATE(integ_dos)
 ABI_DEALLOCATE(total_dos)
 ABI_DEALLOCATE(total_integ_dos)

 if (m_dos_flag>=1)  then
   ABI_DEALLOCATE(partial_dos_m)
   ABI_DEALLOCATE(integ_dos_m)
   ABI_DEALLOCATE(total_dos_m)
   ABI_DEALLOCATE(total_integ_dos_m)
 end if
 if (paw_dos_flag==1)  then
   ABI_DEALLOCATE(total_dos_paw1)
   ABI_DEALLOCATE(total_dos_pawt1)
 end if

!DEBUG
!write(std_out,*)' gaus_dos : exit '
!ENDDEBUG

end subroutine gaus_dos
!!***

!!****f* m_pjedos/get_dos_1band
!! NAME
!! get_dos_1band
!!
!! FUNCTION
!! calculate DOS from tetrahedron method for 1 band and 1 sppol
!!
!! INPUTS
!! dos_fractions=fractional DOS at each irred kpoint
!! enemin=minimal energy for DOS
!! enemax=maximal energy for DOS
!! nene=number of energies for DOS
!! nkpt=number of irreducible kpoints
!! ndosfraction=number of different fractional DOSs
!! tweight=sum of tetrahedron weights for each irred kpoint
!! dtweightde=energy derivative of tweight
!!
!! OUTPUT
!!  partial_dos(nene,ndosfraction)=partial DOS, for each different channel
!!  integ_dos(nene,ndosfraction)=integrated DOS, for each different channel
!!
!! PARENTS
!!      tetrahedron
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_dos_1band (dos_fractions,enemin,enemax,&
&            integ_dos,nene,nkpt,ndosfraction,&
&            partial_dos,tweight,dtweightde)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_dos_1band'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndosfraction,nene,nkpt
 real(dp),intent(in) :: enemax,enemin
!arrays
 real(dp),intent(in) :: dos_fractions(nkpt,ndosfraction),dtweightde(nkpt,nene)
 real(dp),intent(in) :: tweight(nkpt,nene)
 real(dp),intent(out) :: integ_dos(nene,ndosfraction)
 real(dp),intent(out) :: partial_dos(nene,ndosfraction)

!Local variables-------------------------------
!scalars
 integer :: ieps,ifract,ikpt
 real(dp) :: deltaene,eps

! *********************************************************************

 partial_dos(:,:) = zero
 integ_dos(:,:) = zero

 deltaene = (enemax-enemin) / (nene-1)
!
!Calculate parameters of DOS at each point eps in [epsmin,epsmax]
!
 eps=enemin
 do ieps=1,nene

   do ifract=1,ndosfraction
     do ikpt=1,nkpt
       partial_dos(ieps,ifract) = partial_dos(ieps,ifract)+&
&       dtweightde(ikpt,ieps)*dos_fractions(ikpt,ifract)
       integ_dos(ieps,ifract) = integ_dos(ieps,ifract)+&
&       tweight(ikpt,ieps)*dos_fractions(ikpt,ifract)
     end do
   end do
   eps = eps + deltaene
 end do

end subroutine get_dos_1band
!!***

!!****f* m_pjedos/get_dos_1band_m
!! NAME
!! get_dos_1band_m
!!
!! FUNCTION
!! calculate DOS from tetrahedron method for 1 band and 1 sppol
!!
!! INPUTS
!! dos_fractions_m=fractional DOS at each irred kpoint
!! enemin=minimal energy for DOS
!! enemax=maximal energy for DOS
!! nene=number of energies for DOS
!! nkpt=number of irreducible kpoints
!! ndosfraction_m=number of different fractional DOSs
!! tweight=sum of tetrahedron weights for each irred kpoint
!! dtweightde=energy derivative of tweight
!!
!! OUTPUT
!!  partial_dos_m(nene,ndosfraction_m)=partial DOS, for each different channel
!!  integ_dos_m_m(nene,ndosfraction_m)=integrated DOS, for each different channel
!!
!! PARENTS
!!      tetrahedron
!!
!! CHILDREN
!!
!! SOURCE

subroutine get_dos_1band_m (dos_fractions_m,enemin,enemax,&
&            integ_dos_m,nene,nkpt,ndosfraction_m,&
&            partial_dos_m,tweight,dtweightde)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'get_dos_1band_m'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndosfraction_m,nene,nkpt
 real(dp),intent(in) :: enemax,enemin
!arrays
 real(dp),intent(in) :: dos_fractions_m(nkpt,ndosfraction_m)
 real(dp),intent(in) :: dtweightde(nkpt,nene),tweight(nkpt,nene)
 real(dp),intent(out) :: integ_dos_m(nene,ndosfraction_m)
 real(dp),intent(out) :: partial_dos_m(nene,ndosfraction_m)

!Local variables-------------------------------
!scalars
 integer :: ieps,ifract,ikpt
 real(dp) :: deltaene,eps
!arrays

! *********************************************************************

 partial_dos_m(:,:) = zero
 integ_dos_m(:,:) = zero

 deltaene = (enemax-enemin) / (nene-1)
!
!Calculate parameters of DOS at each point eps in [epsmin,epsmax]
!
 eps=enemin
 do ieps=1,nene

   do ifract=1,ndosfraction_m
     do ikpt=1,nkpt
       partial_dos_m(ieps,ifract) = partial_dos_m(ieps,ifract)+&
&       dtweightde(ikpt,ieps)*dos_fractions_m(ikpt,ifract)
       integ_dos_m(ieps,ifract) = integ_dos_m(ieps,ifract)+&
&       tweight(ikpt,ieps)*dos_fractions_m(ikpt,ifract)
     end do
   end do
   eps = eps + deltaene
 end do

end subroutine get_dos_1band_m
!!***

!!****f* m_pjedos/dos_degeneratewfs
!! NAME
!! dos_degeneratewfs
!!
!! FUNCTION
!! Average the contribution to a generalised DOS (that includes weighting factors from matrix elements)
!! from degenerate wavefunctions.
!! Warning : with a negative value of degeneracy_tol, this routine is desactivated. See KNOWN_PROBLEMS.
!!
!! INPUTS
!!  dos_fractions_in(nkpt,mband,nsppol,ndos)= contribution to generalized DOS from each wavefunction
!!  mband =maximum number of electronic bands for each kpoint
!!  nband(nkpt*nsppol)  =number of electronic bands for each kpoint
!!  nkpt         =number of irreducible kpoints
!!  eigen(mband*nkpt*nsppol)=eigenvalues at irred kpoints
!!  ndos= number of components of dos_fractions_in and dos_fractions_out
!!
!! OUTPUT
!!  dos_fractions_out= contribution to generalized DOS from each wavefunction
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dos_degeneratewfs(dos_fractions_in,dos_fractions_out,eigen,mband,nband,ndos,nkpt,nsppol)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dos_degeneratewfs'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,ndos,nkpt,nsppol
 integer,intent(in) :: nband(nkpt*nsppol)
!arrays
 real(dp),intent(in) :: dos_fractions_in(nkpt,mband,nsppol,ndos)
 real(dp),intent(out) :: dos_fractions_out(nkpt,mband,nsppol,ndos)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer :: degeneracy_index,degeneracy_min,iband,idos,ikpt,index_eig,isppol,nband_ksp
!integer :: test_done
 real(dp) :: average,degeneracy_tol=-tol6  ! With a negative value, the routine is desactivated ...
 logical :: degeneracy_flag


! *********************************************************************

!DEBUG
!write(std_out,*) 'dos_degeneratewfs : enter'
!test_done=0
!ENDDEBUG
 

!Initialization of dos_fractions_out
 dos_fractions_out(:,:,:,:)=dos_fractions_in(:,:,:,:)

!Average on the degenerate contributions
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_ksp=nband(ikpt+(isppol-1)*nkpt)
     index_eig=mband*(ikpt-1 + nkpt*(isppol-1))
     if(nband_ksp>1)then
       degeneracy_index=0 ; degeneracy_min=0
       do iband=2,nband_ksp
!        DEBUG
!        if(test_done==20)exit
!        ENDDEBUG
         degeneracy_flag=(eigen(iband+index_eig)-eigen(iband-1+index_eig) < degeneracy_tol)
         if(degeneracy_flag)degeneracy_index=degeneracy_index+1

!        A non-zero value of degeneracy_min will signal that it is time to perform the averaging
         if(degeneracy_index/=0 .and. (.not.degeneracy_flag))degeneracy_min=iband-degeneracy_index-1
         if(degeneracy_index/=0 .and. iband==nband_ksp)      degeneracy_min=iband-degeneracy_index

!        Perform average over all degenerate states at the end of a series of degenerate states
         if(degeneracy_min/=0)then
           do idos=1,ndos
             average=sum(dos_fractions_in(ikpt,degeneracy_min:degeneracy_min+degeneracy_index,isppol,idos))&
&             /real(degeneracy_index+1)
             dos_fractions_out(ikpt,degeneracy_min:degeneracy_min+degeneracy_index,isppol,idos)=average
           end do
!          DEBUG     
!          test_done=test_done+1
!          write(std_out,*)' dos_degeneratewfs '
!          write(std_out,*)' ikpt,isppol,iband,degeneracy_min,degeneracy_max=',&
!          &                      ikpt,isppol,iband,degeneracy_min,degeneracy_min+degeneracy_index
!          write(std_out,*)' eigen(1+index_eig:nband_ksp+index_eig)=',eigen(1+index_eig:nband_ksp+index_eig)
!          write(std_out,*)' dos_fractions_in(ikpt,:,isppol,idos)=',dos_fractions_in(ikpt,:,isppol,idos)
!          write(std_out,*)' dos_fractions_out(ikpt,:,isppol,idos)=',dos_fractions_out(ikpt,:,isppol,idos)
!          ENDDEBUG   

!          Reset degeneracy_index and degeneracy_min
           degeneracy_index=0 ; degeneracy_min=0
         end if
       end do ! iband
       
     end if
   end do ! ikpt
 end do ! isppol

!DEBUG
!write(std_out,*)' dos_degeneratewfs : exit '
!ENDDEBUG

end subroutine dos_degeneratewfs
!!***

!!****f* m_pjedos/recip_ylm
!! NAME
!! recip_ylm
!!
!! FUNCTION
!! Project input wavefunctions (real space) on to Ylm
!!
!! INPUTS
!!  bess_fit(mpw,nradintmax,ll) = bessel functions for ll, splined
!!   with arguments $2 \pi |k+G| \Delta r$, for all G vectors in sphere
!!   and all points on radial grid.
!!  cg_1band(2,npw_k)=wavefunction in recip space
!!  istwfk= storage mode of cg_1band
!!  nradint(natsph)=number of points on radial real-space grid for a given atom
!!  nradintmax=dimension of rint array
!!  mlang=maximum angular momentum
!!  mpi_enreg=information about MPI parallelization
!!  natsph=number of atoms around which ang mom projection has to be done 
!!  npw_k=number of plane waves for kpt
!!  ph3d(2,npw_k,natsph)=3-dim structure factors, for each atom and plane wave.
!!  prtsphere= if 1, print a complete analysis of the angular momenta in atomic spheres
!!  rint(nradintmax) = points on radial real-space grid for integration
!!  rc_ylm= 1 for real spherical harmonics. 2 for complex spherical harmonics, 
!!  rmax(natsph)=maximum radius for real space integration sphere
!!  ucvol=unit cell volume in bohr**3.
!!  ylm(mpw,mlang*mlang)=real spherical harmonics for each G and k point
!!  znucl_sph(natsph)=gives the nuclear number for each type of atom
!!
!! OUTPUT
!!  sum_1ll_1atom(mlang,natsph)= projected scalars for each atom and ang. mom.
!!  sum_1lm_1atom(mlang*mlang,natsph)= projected scalars for each atom and LM component.
!!
!! NOTES
!!  ph3d atoms are ordered with atindx -> by typat
!!  The atoms have to be reordered !
!!
!! PARENTS
!!      m_cut3d,partial_dos_fractions
!!
!! CHILDREN
!!      atomdata_from_znucl,dotprod_g
!!
!! SOURCE

subroutine recip_ylm (bess_fit,cg_1band,istwfk,nradint,nradintmax,mlang,mpi_enreg,&
&  mpw,natsph,ntypat,typat,npw_k,ph3d,jlkpgr_intr,prtsphere,&
&  rint,rmax,rc_ylm,sum_1ll_1atom,sum_1lm_1atom,ucvol,ylm,znucl_sph)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'recip_ylm'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwfk,mlang,mpw,natsph,npw_k,nradintmax
 integer,intent(in) :: prtsphere,rc_ylm,ntypat
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: nradint(natsph),typat(natsph)
 real(dp),intent(in) :: bess_fit(mpw,nradintmax,mlang),cg_1band(2,npw_k)
 real(dp),intent(in) :: ph3d(2,npw_k,natsph),rint(nradintmax)
 real(dp),intent(in) :: jlkpgr_intr(npw_k, mlang, ntypat)
 real(dp),intent(in) :: rmax(natsph),ylm(mpw,mlang*mlang)
 real(dp),intent(in) :: znucl_sph(natsph)
 real(dp),intent(out) :: sum_1ll_1atom(mlang,natsph)
 real(dp),intent(out) :: sum_1lm_1atom(mlang*mlang,natsph)

!Local variables-------------------------------
!scalars
 integer,parameter :: option2=2 ! option for dotprod_g
 integer :: illmm, iat, itypat, ipw, ixint , ll , mm , il,im !, jlm
 real(dp) :: doti, dotr, sum_all, integ, dr, llsign1, llsign2, fact
 type(atomdata_t) :: atom
!arrays
 integer, allocatable :: ilang(:), reylmind(:), imylmind(:)
 real(dp) :: sum_1atom(natsph),sum_1ll(mlang),sum_1lm(mlang*mlang)
 real(dp) :: tmppsia(2,npw_k),tmppsim(2,npw_k)
 real(dp) :: vect(2,npw_k),rint2(nradintmax)
 real(dp), allocatable :: coef1(:),coef2(:),mmsign(:) 

! *************************************************************************

 sum_1lm_1atom = zero
 sum_1ll_1atom = zero
 vect = zero

 do ixint=1,nradintmax
   rint2(ixint) = rint(ixint)**2
 end do

 !if (rc_ylm == 2) then
   ABI_ALLOCATE(ilang,    (mlang**2))
   ABI_ALLOCATE(reylmind, (mlang**2))
   ABI_ALLOCATE(imylmind, (mlang**2))
   ABI_ALLOCATE(mmsign,   (mlang**2))
   ABI_ALLOCATE(coef1,   (mlang**2))
   ABI_ALLOCATE(coef2,   (mlang**2))
   do ll=0,mlang-1
     llsign1 = one
     llsign2 = zero
     if (mod(ll,2) == 1) then
       llsign1 = zero
       llsign2 = one
     end if
     do mm = -ll, -1
       illmm = (ll+1)**2-ll+mm
       ilang(illmm) = ll+1
       coef1(illmm) = llsign1 ! 1 for even ll channels
       coef2(illmm) = llsign2 ! 1 for odd ll channels
       reylmind(illmm) = (ll+1)**2-ll-mm
       imylmind(illmm) = (ll+1)**2-ll+mm
       mmsign(illmm) = -one
     end do

     do mm = 0, ll
       illmm = (ll+1)**2-ll+mm
       ilang(illmm) = ll+1
       coef1(illmm) = llsign1 ! 1 for even ll channels
       coef2(illmm) = llsign2 ! 1 for odd ll channels
       reylmind(illmm) = (ll+1)**2-ll+mm
       imylmind(illmm) = (ll+1)**2-ll-mm
       mmsign(illmm) = one
     end do
   end do
   if (istwfk == 1) then
     coef1 = one
     coef2 = one
   end if
 !end if

!#define NEW

#ifdef NEW
 do iat=1,natsph
   itypat = typat(iat)

   ! Temporary array for part which depends only on iat
   !do ipw=1,npw_k
   !  tmppsia(1,ipw) = cg_1band(1,ipw) * ph3d(1,ipw,iat) - cg_1band(2,ipw) * ph3d(2,ipw,iat)
   !  tmppsia(2,ipw) = cg_1band(1,ipw) * ph3d(2,ipw,iat) + cg_1band(2,ipw) * ph3d(1,ipw,iat)
   !end do

   ! Loop over LM.
   do illmm=1,mlang*mlang
     il = ilang(illmm)

     if (rc_ylm == 2) then
       MSG_ERROR("Not coded")

     else if (rc_ylm == 1) then
       ! to get PDOS for real spherical harmonics, multiply here by ylm instead of linear combination
       dotr = zero; doti = zero
       do ipw=1,npw_k
         fact = ylm(ipw, illmm) * jlkpgr_intr(ipw, il, itypat) 
         tmppsim(1,ipw) = fact * (cg_1band(1, ipw) * ph3d(1, ipw, iat) - cg_1band(2, ipw) * ph3d(2, ipw, iat))
         tmppsim(2,ipw) = fact * (cg_1band(1, ipw) * ph3d(2, ipw, iat) + cg_1band(2, ipw) * ph3d(1, ipw, iat))
         dotr = dotr + tmppsim(1,ipw) 
         doti = doti + tmppsim(2,ipw) 
       end do

       !call dotprod_g(dotr, doti, istwfk, npw_k, option2, tmppsim, tmppsim, mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)
       !dotr = sum(tmppsim(1,:), dim=2); doti = sum(tmppsim(2,:), dim=2)
       !dotr = sum(tmppsim(1,:)); doti = sum(tmppsim(2,:))

     else 
       MSG_ERROR("Wrong value for rc_ylm")
     end if

     sum_1lm_1atom(illmm, iat) = dotr**2 + doti**2
     sum_1ll_1atom(il, iat) = sum_1ll_1atom(il, iat) + dotr**2 + doti**2
   end do ! illmm
 end do ! iat

 sum_1lm_1atom = four_pi**2 * sum_1lm_1atom / ucvol
 sum_1ll_1atom = four_pi**2 * sum_1ll_1atom / ucvol

#else
 ! Big loop on all atoms
 do iat=1,natsph
   dr = rmax(iat)/(nradint(iat)-1)

   ! Temporary array for part which depends only on iat
   do ipw=1,npw_k
     tmppsia(1,ipw) = cg_1band(1,ipw) * ph3d(1,ipw,iat) - cg_1band(2,ipw) * ph3d(2,ipw,iat)
     tmppsia(2,ipw) = cg_1band(1,ipw) * ph3d(2,ipw,iat) + cg_1band(2,ipw) * ph3d(1,ipw,iat)
   end do

   ! tmppsim = temporary arrays for part of psi which doesnt depend on ixint
   ! Take into account the fact that ylm are REAL spherical harmonics, see initylmg.f
   ! For time-reversal states, detailed treatment show that only the real or imaginary
   ! part of tmppsia is needed here, depending on l being even or odd: only one of the coef is 1, the other 0
   do illmm=1, mlang*mlang
     !ll = ilang(illmm) - 1
     !mm = illmm - (ll+1)**2 + ll

     if (rc_ylm == 2) then

       do ipw=1,npw_k
         ! to get PDOS for complex spherical harmonics, build linear combination of real ylm 
         ! TODO: check the prefactor part for istwfk /= 1! Could also be incorrect if we go to real spherical harmonics
         tmppsim(1, ipw) =  coef1(illmm) * tmppsia(1, ipw) * ylm(ipw, reylmind(illmm)) &
&         + mmsign(illmm) * coef2(illmm) * tmppsia(2, ipw) * ylm(ipw, imylmind(illmm))

         tmppsim(2, ipw) =  coef2(illmm) * tmppsia(2, ipw) * ylm(ipw, reylmind(illmm)) &
&         - mmsign(illmm) * coef1(illmm) * tmppsia(1, ipw) * ylm(ipw, imylmind(illmm))
       end do

     else if (rc_ylm == 1) then
       ! to get PDOS for real spherical harmonics, multiply here by ylm instead of linear combination
       do ipw=1,npw_k
         tmppsim(1,ipw) = tmppsia(1,ipw) * ylm(ipw,illmm)
         tmppsim(2,ipw) = tmppsia(2,ipw) * ylm(ipw,illmm)
       end do

       ! Handle time-reversal
       if (istwfk /= 1) then
         if (mod(ilang(illmm) - 1, 2) == 0) then
            tmppsim(2,:) = zero
         else
            tmppsim(1,:) = tmppsim(2,:)
            tmppsim(2,:) = zero
         end if
       end if

     else 
       MSG_ERROR("Wrong value for rc_ylm")
     end if

     integ = zero
     do ixint = 1, nradint(iat)
       vect(1, 1:npw_k) = bess_fit(1:npw_k, ixint, ilang(illmm))
       call dotprod_g(dotr, doti, istwfk, npw_k, option2, vect, tmppsim, mpi_enreg%me_g0, mpi_enreg%comm_spinorfft)

!      Multiply by r**2 and take norm, integrate
!      MJV 5.5.2012 removed call to simpson_int - just need last value, 
!      no need to allocate full space for primitive and integrand
       if (ixint ==1 .or. ixint == nradint(iat)) then
         integ = integ + 0.375_dp * rint2(ixint) * (dotr**2 + doti**2)
       else if (ixint ==2 .or. ixint == nradint(iat)-1) then
         integ = integ + 1.166666666666666666666666667_dp * rint2(ixint) * (dotr**2 + doti**2)
       else if (ixint ==3 .or. ixint == nradint(iat)-2) then
         integ = integ + 0.958333333333333333333333333_dp * rint2(ixint) * (dotr**2 + doti**2)
       else
         integ = integ + rint2(ixint) * (dotr**2 + doti**2)
       end if
     end do ! ixint
     integ = integ * dr

     sum_1lm_1atom(illmm, iat)        = sum_1lm_1atom(illmm, iat)        + integ
     sum_1ll_1atom(ilang(illmm), iat) = sum_1ll_1atom(ilang(illmm), iat) + integ
   end do ! illmm
 end do ! iat


!MJV 5.5.2012: removed 1/sqrt(2) above in tmppsim and (4 pi)^2 in integrand - just multiply by 8 pi^2 here
!Normalize with unit cell volume
 if (rc_ylm == 2) then
   sum_1lm_1atom(:,:) = eight * pi**2 * sum_1lm_1atom(:,:) / ucvol
   sum_1ll_1atom(:,:) = eight * pi**2 * sum_1ll_1atom(:,:) / ucvol
 else
   sum_1lm_1atom(:,:) = four_pi**2 * sum_1lm_1atom(:,:) / ucvol
   sum_1ll_1atom(:,:) = four_pi**2 * sum_1ll_1atom(:,:) / ucvol
 end if
#endif

 !if (rc_ylm == 2) then
   ABI_DEALLOCATE(coef1)
   ABI_DEALLOCATE(coef2)
   ABI_DEALLOCATE(reylmind)
   ABI_DEALLOCATE(imylmind)
   ABI_DEALLOCATE(mmsign)
   ABI_DEALLOCATE(ilang)
 !end if

 !do iat=1,natsph
 !  do il=0,mlang-1
 !    do im=1,il
 !      sum_1lm_1atom(il**2+il+1+im, iat) = half * &
 !      (sum_1lm_1atom(il**2+il+1+im, iat) + sum_1lm_1atom(il**2+il+1-im, iat))

 !      sum_1lm_1atom(il**2+il+1-im, iat) = sum_1lm_1atom(il**2+il+1+im, iat)
 !    end do
 !  end do
 !end do ! iat

!Output
 if (prtsphere == 1) then
   sum_1ll = zero
   sum_1lm = zero
   sum_1atom = zero
   do iat=1,natsph
     sum_1atom(iat) = sum(sum_1lm_1atom(:,iat))
     sum_1ll(:)=sum_1ll(:)+sum_1ll_1atom(:,iat)
     sum_1lm(:)=sum_1lm(:)+sum_1lm_1atom(:,iat)
   end do
   sum_all = sum(sum_1atom)

   write(std_out,'(a)' ) ' Angular analysis '
   do iat=1,natsph
     call atomdata_from_znucl(atom, znucl_sph(iat))
     write(std_out,'(a)' ) ' '
     write(std_out,'(a,i3,a,a,a,f10.6)' )&
&     ' Atom # ',iat, ' is  ',  atom%symbol,', in-sphere charge =',sum_1atom(iat)
     do ll=0,mlang-1
       write(std_out,'(a,i1,a,f9.6,a,9f6.3)' )&
&       ' l=',ll,', charge=',sum_1ll_1atom(ll+1,iat),&
&       ', m=-l,l splitting:',sum_1lm_1atom(1+ll**2:(ll+1)**2,iat)
     end do ! ll
   end do ! iat
   write(std_out,'(a,a)') ch10,' Sum of angular contributions for all atomic spheres '
   do ll=0,mlang-1
     write(std_out,'(a,i1,a,f9.6,a,f9.6)' )&
&     ' l=',ll,', charge =',sum_1ll(ll+1),' proportion =',sum_1ll(ll+1)/sum_all
   end do
   write(std_out,'(a,a,f10.6)' ) ch10,' Total over all atoms and l=0 to 4 :',sum_all
   write(std_out,'(a)' ) ' '
 end if

end subroutine recip_ylm
!!***

end module m_epjdos
!!***
