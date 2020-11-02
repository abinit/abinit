!!****m* ABINIT/m_paw_tools
!! NAME
!!  m_paw_tools
!!
!! FUNCTION
!!  This module contains miscelaneous routines used in the PAW context.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2020 ABINIT group (FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_paw_tools

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_dtset

 use m_paral_atom,       only : get_my_atmtab, free_my_atmtab
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype,EP_POSITRON
 use m_pawang,           only : pawang_type, mat_slm2ylm, mat_mlms2jmj
 use m_pawtab,           only : pawtab_type
 use m_paw_ij,           only : paw_ij_type, paw_ij_free, paw_ij_nullify, paw_ij_gather
 use m_pawdij,           only : pawdij_print_dij
 use m_pawrhoij,         only : pawrhoij_type, pawrhoij_free, pawrhoij_gather, pawrhoij_nullify, &
&                               pawrhoij_print_rhoij
 use m_paw_io,           only : pawio_print_ij
 use m_paw_sphharm,      only : mat_mlms2jmj, mat_slm2ylm
 use m_paw_correlations, only : setnoccmmp

 implicit none

 private

!public procedures.
 public :: chkpawovlp
 public :: pawprt

CONTAINS  !========================================================================================
!!***

!----------------------------------------------------------------------

!!****f* m_paw_tools/chkpawovlp
!! NAME
!! chkpawovlp
!!
!! FUNCTION
!! Verify that the PAW spheres are not overlapping
!!
!! INPUTS
!!  natom=number of atoms in cell.
!!  ntypat=number of types of atoms in unit cell.
!!  pawovlp=percentage of voluminal overlap ratio allowed to continue execution
!!          (if negative value, execution always continues)
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data:
!!  rmet(3,3)=real space metric ($\textrm{bohr}^{2}$).
!!  typat(natom)=type (integer) for each atom
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  (only checking)
!!
!! NOTES
!!
!! PARENTS
!!      m_bethe_salpeter,m_nonlinear,m_respfn_driver,m_scfcv_core
!!      m_screening_driver,m_sigma_driver,m_wfk_analyze
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,mat_mlms2jmj,mat_slm2ylm,paw_ij_free
!!      paw_ij_gather,paw_ij_nullify,pawdij_print_dij,pawrhoij_free
!!      pawrhoij_gather,pawrhoij_nullify,pawrhoij_print_rhoij,setnoccmmp,wrtout
!!      xmpi_comm_group,xmpi_group_free,xmpi_group_translate_ranks
!!
!! SOURCE

subroutine chkpawovlp(natom,ntypat,pawovlp,pawtab,rmet,typat,xred)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat
 real(dp) :: pawovlp
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: rmet(3,3),xred(3,natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: ia,ib,ii,t1,t2,t3
 logical :: stop_on_error
 real(dp) :: dd,dif1,dif2,dif3,ha,hb,norm2
 real(dp) :: ratio_percent,va,vb,vv
 character(len=750) :: message
!arrays
 integer :: iamax(2),ibmax(2),iovl(2)
 real(dp) :: norm2_min(2),r2cut(2),ratio_percent_max(2),rcuta(2),rcutb(2)


! *************************************************************************

 DBG_ENTER("COLL")

 iamax(:)=-1;ibmax(:)=-1
 norm2_min(:)=-1.d0;ratio_percent_max(:)=-1.d0
 iovl(:)=0

!Loop on "overlapping" atoms with the maximum overlap
 do ia=1,natom

   rcuta(1)=pawtab(typat(ia))%rpaw
   rcuta(2)=pawtab(typat(ia))%rshp

   do ib=ia,natom

     rcutb(1)=pawtab(typat(ib))%rpaw
     rcutb(2)=pawtab(typat(ib))%rshp
     r2cut(1)=(rcuta(1)+rcutb(1))**2
     r2cut(2)=(rcuta(2)+rcutb(2))**2

!    Visit the box and its first images:
     do t3=-1,1
       do t2=-1,1
         do t1=-1,1

           dif1=xred(1,ia)-(xred(1,ib)+dble(t1))
           dif2=xred(2,ia)-(xred(2,ib)+dble(t2))
           dif3=xred(3,ia)-(xred(3,ib)+dble(t3))
           norm2=sqnrm_pawovlp(dif1,dif2,dif3)

           do ii=1,2

             if(norm2>tol10.and.norm2<r2cut(ii)) then

               iovl(ii)=iovl(ii)+1

!              Compute the overlap ratio:
               dd=sqrt(norm2)
               va=4._dp/3._dp*pi*rcuta(ii)**3
               vb=4._dp/3._dp*pi*rcutb(ii)**3
               ha=(rcutb(ii)**2-(dd-rcuta(ii))**2)/(two*dd)
               hb=(rcuta(ii)**2-(dd-rcutb(ii))**2)/(two*dd)
               vv=pi/3.d0*(ha**2*(three*rcuta(ii)-ha)+hb**2*(three*rcutb(ii)-hb))
               ratio_percent=100._dp*min(vv/min(va,vb),one)
               if (ratio_percent>ratio_percent_max(ii)) then
                 ratio_percent_max(ii)=ratio_percent
                 norm2_min(ii)=norm2
                 iamax(ii)=ia;ibmax(ii)=ib
               end if

             end if
           end do
         end do
       end do
     end do
   end do
 end do

 stop_on_error=(abs(pawovlp)<=tol6.or.(pawovlp>tol6.and.ratio_percent_max(1)>pawovlp))

!Print adapted message with overlap value
 if (iovl(1)+iovl(2)>0) then

   !ii=1: PAW augmentation regions overlap
   !ii=2: compensation charges overlap
   if (iovl(2)==0) ii=1
   if (iovl(2)> 0) ii=2

   if (iovl(ii)>0) then

     if (ii==1) write(message,' (a)' ) 'PAW SPHERES ARE OVERLAPPING!'
     if (ii==2) write(message, '(2a)' )'PAW COMPENSATION DENSITIES ARE OVERLAPPING !!!!'

     if (iovl(ii)==1) then
       write(message, '(3a)' ) trim(message),ch10,&
&       '   There is one pair of overlapping atoms.'
     else
       write(message, '(3a,i5,a)' ) trim(message),ch10,&
&       '   There are ', iovl(1),' pairs of overlapping atoms.'
     end if
     write(message, '(3a,i3,a,i3,a)' ) trim(message),ch10,&
     '   The maximum overlap percentage is obtained for the atoms ',iamax(ii),' and ',ibmax(ii),'.'
     write(message, '(2a,2(a,i3),a,f9.5,a,2(a,i3,a,f9.5,a),a,f5.2,a)' ) trim(message),ch10,&
&     '    | Distance between atoms ',iamax(ii),' and ',ibmax(ii),' is  : ',sqrt(norm2_min(ii)),ch10,&
&     '    | PAW radius of the sphere around atom ',iamax(ii),' is: ',pawtab(typat(iamax(ii)))%rpaw,ch10,&
&     '    | PAW radius of the sphere around atom ',ibmax(ii),' is: ',pawtab(typat(ibmax(ii)))%rpaw,ch10,&
&     '    | This leads to a (voluminal) overlap ratio of ',ratio_percent_max(ii),' %'
     if (ii==2) then
       write(message, '(3a)' ) trim(message),ch10,&
&       'THIS IS DANGEROUS !, as PAW formalism assumes non-overlapping compensation densities.'
     end if

     if (stop_on_error) then
       MSG_ERROR_NOSTOP(message,ia) !ia is dummy
     else
       MSG_WARNING(message)
     end if

   end if

!  Print advice
   if (stop_on_error) then
     write(message, '(3a)' )&
&     '  Action: 1- decrease cutoff radius of PAW dataset',ch10,&
&     '    OR  2- ajust "pawovlp" input variable to allow overlap (risky)'
     MSG_ERROR(message)
   end if

!  Print last message if execution continues:
   if (pawovlp<=tol6) then
     write(message, '(6a)' ) &
&     '       Results might be approximate,',ch10,&
&     '       and even inaccurate (if overlap is too big) !',ch10,&
&     '       Assume experienced user. Execution will continue.',ch10
     call wrtout(std_out,message,'COLL')
   else if (ratio_percent_max(1)<=pawovlp) then
     write(message, '(8a)' ) &
&     '       Overlap ratio seems to be acceptable (less than value',ch10,&
&     '       of "pawovlp" input parameter): execution will continue.',ch10,&
&     '       But be aware that results might be approximate,',ch10,&
&     '       and even inaccurate (depending on your physical system) !',ch10
     call wrtout(std_out,message,'COLL')
   end if

 end if !iovl>0

 DBG_EXIT("COLL")

 contains

   function sqnrm_pawovlp(u1,u2,u3)
!squared norm of a vector
   real(dp) :: sqnrm_pawovlp
   real(dp),intent(in) :: u1,u2,u3

   sqnrm_pawovlp=rmet(1,1)*u1*u1+rmet(2,1)*u2*u1+rmet(3,1)*u3*u1&
&   +rmet(1,2)*u1*u2+rmet(2,2)*u2*u2+rmet(3,2)*u3*u2&
&   +rmet(1,3)*u1*u3+rmet(2,3)*u2*u3+rmet(3,3)*u3*u3

 end function sqnrm_pawovlp

end subroutine chkpawovlp
!!***

!----------------------------------------------------------------------

!!****f* m_paw_tools/pawprt
!! NAME
!! pawprt
!!
!! FUNCTION
!! Print out data concerning PAW formalism
!! (pseudopotential strength, augmentation occupancies...)
!! To be called at the end of the SCF cycle
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (FJ,MT,BA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | enunit=parameter determining units of output energies
!!   | kptopt=option for the generation of k points
!!   | natom=number of atoms in cell
!!   | ntypat = number of atom types
!!   | pawprtvol= printing volume
!!   | pawspnorb=flag: 1 if spin-orbit coupling is activated
!!   | typat(natom)=type of each atom
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  paw_ij(my_natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  (only printing)
!!
!! PARENTS
!!      m_bethe_salpeter,m_outscfcv,m_screening_driver,m_sigma_driver
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,mat_mlms2jmj,mat_slm2ylm,paw_ij_free
!!      paw_ij_gather,paw_ij_nullify,pawdij_print_dij,pawrhoij_free
!!      pawrhoij_gather,pawrhoij_nullify,pawrhoij_print_rhoij,setnoccmmp,wrtout
!!      xmpi_comm_group,xmpi_group_free,xmpi_group_translate_ranks
!!
!! SOURCE

subroutine pawprt(dtset,my_natom,paw_ij,pawrhoij,pawtab,&
&                 electronpositron,& ! optional argument
&                 mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom
 integer,optional,intent(in) :: comm_atom
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer,optional :: electronpositron
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 type(paw_ij_type),target,intent(inout) :: paw_ij(my_natom)
 type(pawrhoij_type),target,intent(inout) :: pawrhoij(my_natom)
 type(pawtab_type),target,intent(in) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer,parameter :: natmax=2
 integer :: cplex_dij,group1,group2,iat,iatom,ierr,ii,im1,im2,ipositron,ispden
 integer :: i_unitfi,itypat,ll,llp,me_atom,my_comm_atom,natprt,ndij,nspden,nsppol
 integer :: unitfi,unt
 real(dp) :: mnorm,mx,my,mz,ntot,valmx,localm
 logical :: my_atmtab_allocated,paral_atom,useexexch,usepawu
 type(pawang_type):: pawang_dum
 character(len=7),parameter :: dspin1(6)=(/"up     ","down   ","up-up  ","dwn-dwn","up-dwn ","dwn-up "/)
 character(len=8),parameter :: dspin2(6)=(/"up      ","down    ","dens (n)","magn (x)","magn (y)","magn (z)"/)
 character(len=500) :: msg
!arrays
 integer :: idum(1)
 integer :: idum1(0),idum3(0,0,0)
 integer,allocatable :: jatom(:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: rdum2(0,0),rdum4(0,0,0,0)
 complex(dpc),allocatable :: noccmmp_ylm(:,:,:),noccmmp_jmj(:,:),noccmmp_slm(:,:,:)
 type(paw_ij_type), ABI_CONTIGUOUS pointer :: paw_ij_all(:)
 type(pawrhoij_type),ABI_CONTIGUOUS pointer :: pawrhoij_all(:)

! *********************************************************************

 DBG_ENTER("COLL")

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=dtset%natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,dtset%natom,my_natom_ref=my_natom)
 me_atom=xmpi_comm_rank(my_comm_atom)

!Continue only if comm_atom contains the master of the output comm
 if (paral_atom) then
   call xmpi_comm_group(abinit_comm_output,group1,ierr)
   call xmpi_comm_group(my_comm_atom,group2,ierr)
   call xmpi_group_translate_ranks(group1,1,(/0/),group2,idum,ierr)
   call xmpi_group_free(group1)
   call xmpi_group_free(group2)
   if (idum(1)==xmpi_undefined) then
     call free_my_atmtab(my_atmtab,my_atmtab_allocated)
     return
   end if
 end if

!Initializations
 natprt=natmax;if (dtset%natom==1) natprt=1
 if (dtset%pawprtvol<0) natprt=dtset%natom
 ABI_ALLOCATE(jatom,(natprt))
 if (natprt==1) then
   jatom(1)=1
 else if (natprt==2) then
   jatom(1)=1;jatom(2)=dtset%natom
 else if (natprt==dtset%natom) then
   do iat=1,dtset%natom
     jatom(iat)=iat
   end do
 else
   MSG_BUG("invalid value of natprt!")
 end if
 usepawu=(count(pawtab(:)%usepawu/=0)>0)
 useexexch=(count(pawtab(:)%useexexch/=0)>0)
 ipositron=0
 if (present(electronpositron)) then
   if (associated(electronpositron)) ipositron=electronpositron%calctype
 end if

!Main title
 write(msg, '(2a)' ) ch10,&
& ' ==== Results concerning PAW augmentation regions ===='
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')
 msg=' '
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

!If atomic data are distributed, retrieve all Dij on master proc
 if (paral_atom) then
   if (me_atom==0) then
     ABI_DATATYPE_ALLOCATE(paw_ij_all,(dtset%natom))
     call paw_ij_nullify(paw_ij_all)
   else
     ABI_DATATYPE_ALLOCATE(paw_ij_all,(0))
   end if
   call paw_ij_gather(paw_ij,paw_ij_all,0,my_comm_atom)
 else
   paw_ij_all => paw_ij
 end if

!Print out pseudopotential strength
!----------------------------------
 if (me_atom==0) then
   do i_unitfi=1,2
     unitfi=ab_out;if (i_unitfi==2) unitfi=std_out
     do unt=1,2
       if (((unt==1).and.(dtset%enunit==0.or.dtset%enunit==2)).or.&
&          ((unt==2).and.(dtset%enunit==1.or.dtset%enunit==2))) then
         if ((unt==1).and.(dtset%enunit==0.or.dtset%enunit==2)) then
           write(msg,'(a)') ' Total pseudopotential strength Dij (hartree):'
         else if ((unt==2).and.(dtset%enunit==1.or.dtset%enunit==2)) then
           write(msg,'(a)') ' Total pseudopotential strength Dij (eV):'
         end if
         call wrtout(unitfi,msg,'COLL')
         if (ipositron>0) then
           if (electronpositron%has_pos_ham==0) then
             write(msg,'(a)') ' -Note: these are the electronic Dij'
           else
             write(msg,'(a)') ' -Note: these are the positronic Dij'
           end if
           call wrtout(unitfi,msg,'COLL')
         end if
         valmx=100._dp;if (ipositron>0) valmx=-1._dp
         do iat=1,natprt
           iatom=jatom(iat)
           call pawdij_print_dij(paw_ij_all(iatom)%dij,paw_ij_all(iatom)%cplex_dij,&
&                  paw_ij_all(iatom)%qphase,iatom,dtset%natom,paw_ij_all(iatom)%nspden,&
&                  test_value=valmx,unit=unitfi,Ha_or_eV=unt,opt_prtvol=dtset%pawprtvol)
         end do
       end if
       msg=' '
       call wrtout(unitfi,msg,'COLL')
     end do
   end do
 end if
 if (paral_atom.and.(.not.usepawu).and.(.not.useexexch)) then
   call paw_ij_free(paw_ij_all)
   ABI_DATATYPE_DEALLOCATE(paw_ij_all)
 end if

!If atomic data are distributed, retrieve all Rhoij on master proc
 if (paral_atom) then
   if (me_atom==0) then
     ABI_DATATYPE_ALLOCATE(pawrhoij_all,(dtset%natom))
   else
     ABI_DATATYPE_ALLOCATE(pawrhoij_all,(0))
   end if
   call pawrhoij_nullify(pawrhoij_all)
   call pawrhoij_gather(pawrhoij,pawrhoij_all,0,my_comm_atom,&
&   with_grhoij=.false.,with_lmnmix=.false.,&
&   with_rhoij_=.false.,with_rhoijres=.false.)
 else
   pawrhoij_all => pawrhoij
 end if

!Print out SYMMETRIZED occupancies of the partial waves
!------------------------------------------------------
 if (me_atom==0) then
   do i_unitfi=1,2
     unitfi=ab_out;if (i_unitfi==2) unitfi=std_out
     write(msg,'(a)') ' Augmentation waves occupancies Rhoij:'
     call wrtout(unitfi,msg,'COLL')
     if (ipositron>0) then
       if (electronpositron%particle==EP_POSITRON) then
         write(msg,'(a)') ' -Note: these are the electronic Rhoij'
       else
         write(msg,'(a)') ' -Note: these are the positronic Rhoij'
       end if
       call wrtout(unitfi,msg,'COLL')
     end if
     if (dtset%pawspnorb>0.and.pawrhoij_all(1)%cplex_rhoij==1.and.dtset%kptopt/=1.and.dtset%kptopt/=2) then
       write(msg,'(6a)') ' pawprt: - WARNING:',ch10,&
&       '       Spin-orbit coupling is activated but only real part of Rhoij occupancies',ch10,&
&       '       has been computed; they could have an imaginary part (not printed here).'
       call wrtout(unitfi,msg,'COLL')
     end if
     valmx=25._dp;if (ipositron>0) valmx=-1._dp
     do iat=1,natprt
       iatom=jatom(iat);nspden=pawrhoij_all(iatom)%nspden
       call pawrhoij_print_rhoij(pawrhoij_all(iatom)%rhoijp,pawrhoij_all(iatom)%cplex_rhoij,&
&                    pawrhoij_all(iatom)%qphase,iatom,dtset%natom,&
&                    rhoijselect=pawrhoij_all(iatom)%rhoijselect,&
&                    test_value=valmx,unit=unitfi,opt_prtvol=dtset%pawprtvol)
     end do
     msg=' '
     call wrtout(unitfi,msg,'COLL')
   end do
 end if

!PAW+U or local exact-exchange: print out +U components of occupancies
!---------------------------------------------------------------------
 if ((usepawu.or.useexexch).and.ipositron/=1.and.me_atom==0) then
   do i_unitfi=1,2
     unitfi=ab_out;if (i_unitfi==2) unitfi=std_out
     if(useexexch) write(msg,'(a)') &
&     ' "Local exact-exchange" part of augmentation waves occupancies Rhoij:'
     if(usepawu) write(msg,'(a)') &
&     ' "PAW+U" part of augmentation waves occupancies Rhoij:'
     call wrtout(unitfi,msg,'COLL')
     do iatom=1,dtset%natom
       itypat=pawrhoij_all(iatom)%itypat
       nspden=pawrhoij_all(iatom)%nspden
       ll=-1;if (pawtab(itypat)%usepawu/=0) ll=pawtab(itypat)%lpawu
       llp=-1;if (pawtab(itypat)%useexexch/=0) llp=pawtab(itypat)%lexexch
       if (ll/=llp.and.ll/=-1.and.llp/=-1) then
         MSG_BUG("lpawu/=lexexch forbidden!")
       end if
       ll=max(ll,llp)
       if (ll>=0) then
         call pawrhoij_print_rhoij(pawrhoij_all(iatom)%rhoijp,pawrhoij_all(iatom)%cplex_rhoij,&
&                      pawrhoij_all(iatom)%qphase,iatom,dtset%natom,&
&                      rhoijselect=pawrhoij_all(iatom)%rhoijselect,&
&                      l_only=ll,indlmn=pawtab(itypat)%indlmn,&
&                      unit=unitfi,opt_prtvol=dtset%pawprtvol)
       end if
     end do ! iatom
     msg=' '
     call wrtout(unitfi,msg,'COLL')
   end do
 end if

!PAW+U: print out occupations for correlated orbitals
!----------------------------------------------------
 if (usepawu.and.ipositron/=1.and.me_atom==0) then
   do i_unitfi=1,2
     unitfi=ab_out;if (i_unitfi==2) unitfi=std_out
     write(msg,'(3a)') &
&     ' ---------- DFT+U DATA --------------------------------------------------- ',ch10
     call wrtout(unitfi,msg,'COLL')
     do iatom=1,dtset%natom
       itypat=dtset%typat(iatom);ll=pawtab(itypat)%lpawu
       nspden=paw_ij_all(iatom)%nspden;ndij=paw_ij_all(iatom)%ndij
       cplex_dij=paw_ij_all(iatom)%cplex_dij
       if ((ll>=0).and.(pawtab(itypat)%usepawu/=0)) then
         write(msg,fmt='(a,i5,a,i4,a)') " ====== For Atom ", iatom,&
&         ", occupations for correlated orbitals. lpawu =",ll,ch10
         call wrtout(unitfi,msg,'COLL')
         if(pawtab(itypat)%usepawu>=10) then
           write(msg,fmt='(a)') "  (This is PAW atomic orbital occupations)"
           call wrtout(unitfi,msg,'COLL')
           write(msg,fmt='(a)') "  (For Wannier orbital occupations, refer to DFT+DMFT occupations above)"
           call wrtout(unitfi,msg,'COLL')
         end if
         if(nspden==2) then
           do ispden=1,nspden
             write(msg,fmt='(a,i4,a,i3,a,f10.5)') " Atom", iatom,&
&             ". Occ. for lpawu and for spin",ispden," =",paw_ij_all(iatom)%nocctot(ispden)
             call wrtout(unitfi,msg,'COLL')
           end do
           localm=paw_ij_all(iatom)%nocctot(2)-paw_ij_all(iatom)%nocctot(1)
           write(msg,fmt='(a,i4,a,2x,f12.6)') " => On atom",iatom,&
&           ",  local Mag. for lpawu is  ",localm
           call wrtout(unitfi,msg,'COLL')
         end if
         if(ndij==4) then
           ntot=paw_ij_all(iatom)%nocctot(1)
           mx=paw_ij_all(iatom)%nocctot(2)
           my=paw_ij_all(iatom)%nocctot(3)
           mz=paw_ij_all(iatom)%nocctot(4)
           mnorm=sqrt(mx*mx+my*my+mz*mz)
           write(msg,'(a,i4,a,2x,e15.8)') " => On atom",iatom,", for  lpawu, local Mag. x is  ",mx
           call wrtout(unitfi,msg,'COLL')
           write(msg,'(14x,a,2x,e15.8)') "               local Mag. y is  ",my
           call wrtout(unitfi,msg,'COLL')
           write(msg,'(14x,a,2x,e15.8)') "               local Mag. z is  ",mz
           call wrtout(unitfi,msg,'COLL')
           write(msg,'(14x,a,2x,e15.8)') "               norm of Mag. is  ",mnorm
           call wrtout(unitfi,msg,'COLL')
           write(msg,fmt='(8x,a,2x,f10.5)') " (along mag axis)    occ. for majority spin is = ",&
&           half*(ntot+mnorm)
           call wrtout(unitfi,msg,'COLL')
           write(msg,fmt='(8x,a,2x,f10.5)') " (along mag axis)    occ. for minority spin is = ",&
&           half*(ntot-mnorm)
           call wrtout(unitfi,msg,'COLL')
         end if
         write(msg,'(3a)') ch10," == Occupation matrix for correlated orbitals:",ch10
         call wrtout(unitfi,msg,'COLL')
         do ispden=1,ndij
           if (nspden==1.and.ndij/=4.and.(cplex_dij==1)) write(msg,fmt='(a)') " Up component only..."
           if (nspden==2) write(msg,fmt='(a,i3)')" Occupation matrix for spin",ispden
           if (ndij==4.or.(cplex_dij==2)) &
&           write(msg,fmt='(2a)')  " Occupation matrix for component ",trim(dspin1(ispden+2*(ndij/4)))
           call wrtout(unitfi,msg,'COLL')
           do im1=1,ll*2+1
             if(cplex_dij==1)&
&             write(msg,'(12(1x,9(1x,f10.5)))') (paw_ij_all(iatom)%noccmmp(1,im1,im2,ispden),im2=1,ll*2+1)
             if(cplex_dij==2)&
&             write(msg,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))')&
&             (paw_ij_all(iatom)%noccmmp(:,im1,im2,ispden),im2=1,ll*2+1)
             call wrtout(unitfi,msg,'COLL')
           end do
           write(msg,'(2a)') ch10,' '
           call wrtout(unitfi,msg,'COLL')
         end do
!        Transformation matrices: real->complex spherical harmonics
         if(paw_ij_all(iatom)%ndij==4) then
           ABI_ALLOCATE(noccmmp_ylm,(2*ll+1,2*ll+1,paw_ij_all(iatom)%ndij))
           noccmmp_ylm=czero
           ABI_ALLOCATE(noccmmp_slm,(2*ll+1,2*ll+1,paw_ij_all(iatom)%ndij))
           noccmmp_slm=czero
!          Go from real notation for complex noccmmp to complex notation in noccmmp_slm
           noccmmp_slm(:,:,:)=cmplx(paw_ij_all(iatom)%noccmmp(1,:,:,:),paw_ij_all(iatom)%noccmmp(2,:,:,:))
           ii=std_out;if (unitfi==ab_out) ii=-1
           call mat_slm2ylm(ll,noccmmp_slm,noccmmp_ylm,paw_ij_all(iatom)%ndij,&
&           1,1,dtset%pawprtvol,ii,'COLL') ! optspin=1 because up spin are first
           do ispden=1,paw_ij_all(iatom)%ndij
             write(msg,'(3a)') ch10,&
&             "== Occupation matrix in the complex harmonics basis for component ",&
&             trim(dspin1(ispden+2*(paw_ij_all(iatom)%ndij/4)))
             call wrtout(unitfi,msg,'COLL')
             do im1=1,ll*2+1
               write(msg,'(12(1x,9(1x,"(",f9.5,",",f9.5,")")))') &
&               (noccmmp_ylm(im1,im2,ispden),im2=1,ll*2+1)
               call wrtout(unitfi,msg,'COLL')
             end do
           end do
           write(msg,'(a)') ch10
           call wrtout(unitfi,msg,'COLL')
           if (dtset%pawspnorb>0) then
             ABI_ALLOCATE(noccmmp_jmj,(2*(2*ll+1),2*(2*ll+1)))
             noccmmp_jmj=czero
             ii=std_out;if (unitfi==ab_out) ii=-1
             call mat_mlms2jmj(ll,noccmmp_ylm,noccmmp_jmj,paw_ij_all(iatom)%ndij,&
&             1,1,dtset%pawprtvol,-1,'COLL') !  optspin=1: up spin are first
             write(msg,'(3a)') ch10,"== Occupation matrix in the J (= L-1/2, L+1/2) and M_J basis"
             call wrtout(unitfi,msg,'COLL')
             do im1=1,2*(ll*2+1)
               write(msg,'(12(1x,18(1x,"(",f7.3,",",f7.3,")")))') &
&               (noccmmp_jmj(im1,im2),im2=1,2*(ll*2+1))
               call wrtout(unitfi,msg,'COLL')
             end do
             write(msg,'(a)') ch10
             call wrtout(unitfi,msg,'COLL')
             ABI_DEALLOCATE(noccmmp_jmj)
           end if ! pawspnorb
           ABI_DEALLOCATE(noccmmp_ylm)
           ABI_DEALLOCATE(noccmmp_slm)
         end if ! ndij==4
       end if ! ((ll>=0).and.(pawtab(itypat)%usepawu/=0))
     end do
   end do
 end if

!Exact exchange: print out occupations for correlated orbitals
!-------------------------------------------------------------
 if (useexexch.and.ipositron/=1.and.me_atom==0) then
   nspden=paw_ij_all(1)%nspden;nsppol=paw_ij_all(1)%nsppol;ndij=paw_ij_all(1)%ndij
   do iatom=1,dtset%natom
     itypat=dtset%typat(iatom);ll=pawtab(itypat)%lexexch
     cplex_dij=paw_ij_all(iatom)%cplex_dij
     if (ll>=0.and.pawtab(itypat)%useexexch/=0) then
       ABI_ALLOCATE(paw_ij_all(iatom)%noccmmp,(cplex_dij,2*ll+1,2*ll+1,ndij))
       ABI_ALLOCATE(paw_ij_all(iatom)%nocctot,(nspden))
     end if
   end do
   call setnoccmmp(1,0,rdum4,0,0,idum3,dtset%natom,dtset%natom,0,1,nsppol,0,dtset%ntypat,&
&   paw_ij_all,pawang_dum,dtset%pawprtvol,pawrhoij_all,pawtab,rdum2,idum1,dtset%typat,1,0)
   do i_unitfi=1,2
     unitfi=ab_out;if (i_unitfi==2) unitfi=std_out
     write(msg, '(3a)' ) &
&     ' ---------- Exact Exchange --------------------------------------------------- ',ch10
     call wrtout(unitfi,msg,'COLL')
     do iatom=1,dtset%natom
       itypat=dtset%typat(iatom);ll=pawtab(itypat)%lexexch
       cplex_dij=paw_ij_all(iatom)%cplex_dij
       if ((ll>=0).and.(pawtab(itypat)%useexexch/=0)) then
         write(msg,fmt='(a,i5,a,i4,a)') " ====== For Atom",iatom,&
&         ", occupations for correlated orbitals. l =",ll,ch10
         call wrtout(unitfi,msg,'COLL')
         do ispden=1,ndij
           if (nspden==1.and.ndij/=4) write(msg,fmt='(a)')   " Up component only..."
           if (nspden==2) write(msg,fmt='(a,i3)')" Occupation matrix for spin",ispden
           if (ndij==4) write(msg,fmt='(2a)')  " Occupation matrix for component ",&
&           trim(dspin2(ispden+2*(ndij/4)))
           call wrtout(unitfi,msg,'COLL')
           do im1=1,ll*2+1
             if(cplex_dij==1)&
&             write(msg,'(12(1x,9(1x,f10.5)))')&
&             (paw_ij_all(iatom)%noccmmp(1,im1,im2,ispden),im2=1,ll*2+1)
             if(cplex_dij==2)&
&             write(msg,'(12(1x,9(1x,"(",f7.3,",",f7.3,")")))') &
&             (paw_ij_all(iatom)%noccmmp(:,im1,im2,ispden),im2=1,ll*2+1)
             call wrtout(unitfi,msg,'COLL')
           end do
           call wrtout(unitfi,' ','COLL')
         end do
       end if
     end do
   end do
   do iatom=1,dtset%natom
     if (allocated(paw_ij_all(iatom)%noccmmp)) then
       ABI_DEALLOCATE(paw_ij_all(iatom)%noccmmp)
     end if
     if (allocated(paw_ij_all(iatom)%nocctot)) then
       ABI_DEALLOCATE(paw_ij_all(iatom)%nocctot)
     end if
   end do
 end if

 msg=' '
 call wrtout(ab_out,msg,'COLL')
 call wrtout(std_out,msg,'COLL')

!Destroy temporary stored atomic data
 ABI_DEALLOCATE(jatom)
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)
 if (paral_atom) then
   call pawrhoij_free(pawrhoij_all)
   ABI_DATATYPE_DEALLOCATE(pawrhoij_all)
   if (usepawu.or.useexexch) then
     call paw_ij_free(paw_ij_all)
     ABI_DATATYPE_DEALLOCATE(paw_ij_all)
   end if
 end if

 DBG_EXIT("COLL")

end subroutine pawprt
!!***

!----------------------------------------------------------------------

END MODULE m_paw_tools
!!***
