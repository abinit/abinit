!{\src2tex{textfont=tt}}
!!****f* ABINIT/tetrahedron
!! NAME
!! tetrahedron
!!
!! FUNCTION
!! calculate DOS by tetrahedron method
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MVer,XG,SM,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!      destroy_tetra,dos_hdr_write,get_dos_1band,get_dos_1band_m
!!      get_full_kgrid,get_tetra_weight,init_tetra,int2char4,matr3inv,metric
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine tetrahedron (dos_fractions,dos_fractions_m,dos_fractions_paw1,dos_fractions_pawt1,&
& dtset,fermie,eigen,fildata,mbesslang,prtdosm,ndosfraction,paw_dos_flag,rprimd)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_xmpi
 use m_tetrahedron
 use m_errors

 use m_io_tools,     only : open_file
 use m_fstrings,     only : int2char4

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'tetrahedron'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_56_recipspace
 use interfaces_61_occeig, except_this_one => tetrahedron
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
! integer :: nfiner
 real(dp) :: buffer,deltaene,enemax,enemin,enex,integral_DOS,max_occ,rcvol
! real(dp) :: tolfermi
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
   message = 'you need at least 2 kpoints to use the tetrahedron method'
   MSG_WARNING(message)
   return
 end if

!Do not support nshiftk > 1: lattice must be decomposed into boxes
!and this is not always possible (I think) with bizzare shiftks
!normally at this point we have incorporated everything into
!kptrlatt, and only 1 shift is needed (in particular for MP grids).

 if (dtset%nshiftk > 1) then
   write(std_out,*) 'tetrahedron: problem with a composite k-point grid.'
   write(std_out,*) '  Only simple lattices are supported Action : use nshiftk=1.'
   write(std_out,*) ' tetrahedron: skip subroutine.'
   write(std_out,*) '  dtset%nshiftk, dtset%shiftk = ', dtset%nshiftk, dtset%shiftk
   write(std_out,*) '  dtset%kptrlatt= ', dtset%kptrlatt
   return
 end if

!Refuse nband different for different kpoints

 do isppol=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
     if ( dtset%nband(dtset%nkpt*(isppol-1) + ikpt) /= dtset%nband(1) ) then
       write(std_out,*) ' tetrahedron : Error nband must be the same for all kpoints'
       write(std_out,*) '  nband = ', dtset%nband
       write(std_out,*) ' tetrahedron: skip subroutine.'
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
   write(std_out,*) ' tetrahedron : no homogeneous grid  of k-points is defined ...'
   write(std_out,*)'  in order to obtain the DOS using the tetrahedron method,'
   write(std_out,*)'  you need to re-define ngkpt or kptrlatt.'
   write(std_out,*) ' tetrahedron: skip subroutine.'
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
!    determine fermi surface contribution from this band,
!    interpolated in irred tetrahedra, on a finer subgrid of klatt
!    NOT YET FUNCTIONAL
!    if(dtset%prtfsurf==1)then
!    tolfermi = dtset%userrc
!    nfiner = dtset%userrb
!    call get_fsurf_1band(dtset,tmp_eigen,fermie,klatt,kpt_fullbz,&
!    &  mtetra,nfiner,nkpt_fullbz,ntetra,tetra_full,tetra_mult,tetra_wrap,tolfermi)
!    end if
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
