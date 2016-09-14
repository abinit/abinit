!{\src2tex{textfont=tt}}
!!****f* ABINIT/gaus_dos
!! NAME
!! gaus_dos
!!
!! FUNCTION
!! Calculate DOS by gauss method
!!
!! COPYRIGHT
!! Copyright (C) 2010-2016 ABINIT group (SM,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine gaus_dos(dos_fractions,&
& dos_fractions_paw1,dos_fractions_pawt1,dtset,fermie,eigen,&
& fildata,mbesslang,m_dos_flag,ndosfraction,paw_dos_flag) 

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_splines

 use m_fstrings, only : int2char4
 use m_io_tools, only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gaus_dos'
 use interfaces_14_hidewrite
 use interfaces_61_occeig, except_this_one => gaus_dos
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
