!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawprt
!! NAME
!! pawprt
!!
!! FUNCTION
!! Print out data concerning PAW formalism
!! (pseudopotential strength, augmentation occupancies...)
!! To be called at the end of the SCF cycle
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (FJ,MT,BA)
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
!!      bethe_salpeter,outscfcv,screening,sigma
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,mat_mlms2jmj,mat_slm2ylm,paw_ij_free
!!      paw_ij_gather,paw_ij_nullify,pawio_print_ij,pawrhoij_free
!!      pawrhoij_gather,pawrhoij_nullify,setnoccmmp,wrtout,xmpi_comm_group
!!      xmpi_group_free,xmpi_group_translate_ranks
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine pawprt(dtset,my_natom,paw_ij,pawrhoij,pawtab,&
&                 electronpositron,& ! optional argument
&                 mpi_atmtab,comm_atom) ! optional arguments (parallelism)


 use defs_basis
 use defs_abitypes
 use defs_parameters
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_paral_atom,       only : get_my_atmtab, free_my_atmtab
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype,EP_POSITRON
 use m_pawang,           only : pawang_type, mat_slm2ylm, mat_mlms2jmj
 use m_pawtab,           only : pawtab_type
 use m_paw_ij,           only : paw_ij_type, paw_ij_free, paw_ij_nullify, paw_ij_gather
 use m_pawrhoij,         only : pawrhoij_type, pawrhoij_free, pawrhoij_gather, pawrhoij_nullify
 use m_paw_io,           only : pawio_print_ij
 use m_paw_sphharm,      only : mat_mlms2jmj, mat_slm2ylm

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pawprt'
 use interfaces_14_hidewrite
 use interfaces_65_paw, except_this_one => pawprt
!End of the abilint section

 implicit none

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
 integer :: group1,group2,iat,iatom,ierr,ii,im1,im2,ipositron,irhoij,ispden
 integer :: i_unitfi,itypat,jrhoij,ll,llp,me_atom,my_comm_atom,natprt,nspden,nsppol
 integer :: optsym,sz1,sz2,unitfi,unt
 real(dp) :: mnorm,mx,my,mz,ntot,valmx,localm
 logical,parameter :: debug=.false.
 logical :: my_atmtab_allocated,paral_atom,useexexch,usepawu
 type(pawang_type):: pawang_dum
 character(len=7),parameter :: dspin1(6)=(/"up     ","down   ","up-up  ","dwn-dwn","up-dwn ","dwn-up "/)
 character(len=8),parameter :: dspin2(6)=(/"up      ","down    ","dens (n)","magn (x)","magn (y)","magn (z)"/)
 character(len=9),parameter :: dspin3(6)=(/"up       ","down     ","up-up    ","down-down","Re[up-dn]","Im[up-dn]"/)
 character(len=500) :: msg0,msg
!arrays
 integer :: idum(1)
 integer :: idum1(0),idum3(0,0,0)
 integer,allocatable :: jatom(:),opt_l_index(:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: rdum2(0,0),rdum4(0,0,0,0)
 real(dp),allocatable :: rhoijs(:,:)
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
 usepawu=(count(pawtab(:)%usepawu>0)>0)
 useexexch=(count(pawtab(:)%useexexch>0)>0)
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
       if ((unt==1).and.(dtset%enunit==0.or.dtset%enunit==2)) then
         write(msg,'(a)') ' Total pseudopotential strength Dij (hartree):'
         call wrtout(unitfi,msg,'COLL')
       else if ((unt==2).and.(dtset%enunit==1.or.dtset%enunit==2)) then
         write(msg,'(a)') ' Total pseudopotential strength Dij (eV):'
         call wrtout(unitfi,msg,'COLL')
       end if
       if (ipositron>0) then
         if (((unt==1).and.(dtset%enunit==0.or.dtset%enunit==2)).or.&
&         ((unt==2).and.(dtset%enunit==1.or.dtset%enunit==2))) then
           if (electronpositron%has_pos_ham==0) then
             write(msg,'(a)') ' -Note: these are the electronic Dij'
           else
             write(msg,'(a)') ' -Note: these are the positronic Dij'
           end if
           call wrtout(unitfi,msg,'COLL')
         end if
       end if
       if (((unt==1).and.(dtset%enunit==0.or.dtset%enunit==2)).or.&
&       ((unt==2).and.(dtset%enunit==1.or.dtset%enunit==2))) then
         do iat=1,natprt
           iatom=jatom(iat);nspden=paw_ij_all(iatom)%ndij
           optsym=2;if (paw_ij_all(iatom)%cplex_dij==2.and.dtset%nspinor==1) optsym=1
           do ispden=1,nspden
             valmx=100._dp;if (ispden==1) valmx=-1._dp
             msg='' ; msg0=''
             if (dtset%natom>1.or.nspden>1) write(msg0, '(a,i3)' ) ' Atom #',iatom
             if (nspden==1) write(msg,'(a)')     trim(msg0)
             if (nspden==2) write(msg,'(2a,i1)') trim(msg0),' - Spin component ',ispden
             if (nspden==4) write(msg,'(3a)')    trim(msg0),' - Component ',trim(dspin1(ispden+2*(nspden/4)))
             if (dtset%natom>1.or.nspden>1) then
               call wrtout(unitfi,msg,'COLL')
             end if
             if (nspden/=4.or.ispden<=2) then
               call pawio_print_ij(unitfi,paw_ij_all(iatom)%dij(:,ispden),paw_ij_all(iatom)%lmn2_size,&
&               paw_ij_all(iatom)%cplex_dij,paw_ij_all(iatom)%lmn_size,-1,idum,0,dtset%pawprtvol,&
&               idum,valmx,unt,opt_sym=optsym)
             else
               call pawio_print_ij(unitfi,paw_ij_all(iatom)%dij(:,ispden),paw_ij_all(iatom)%lmn2_size,&
&               paw_ij_all(iatom)%cplex_dij,paw_ij_all(iatom)%lmn_size,-1,idum,0,dtset%pawprtvol,&
&               idum,valmx,unt,asym_ij=paw_ij_all(iatom)%dij(:,7-ispden),opt_sym=optsym)
             end if
           end do
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
     if (dtset%pawspnorb>0.and.pawrhoij_all(1)%cplex==1.and.dtset%kptopt/=1.and.dtset%kptopt/=2) then
       write(msg,'(6a)') ' pawprt: - WARNING:',ch10,&
&       '       Spin-orbit coupling is activated but only real part of Rhoij occupancies',ch10,&
&       '       has been computed; they could have an imaginary part (not printed here).'
       call wrtout(unitfi,msg,'COLL')
     end if
     do iat=1,natprt
       iatom=jatom(iat);nspden=pawrhoij_all(iatom)%nspden
       optsym=2;if (pawrhoij_all(iatom)%cplex==2.and.dtset%nspinor==1) optsym=1
       do ispden=1,nspden
         valmx=25._dp;if (ispden==1) valmx=-1._dp
         msg='' ; msg0=''
         if (dtset%natom>1.or.nspden>1) write(msg0, '(a,i3)' ) ' Atom #',iatom
         if (nspden==1) write(msg,'(a)')     trim(msg0)
         if (nspden==2) write(msg,'(2a,i1)') trim(msg0),' - Spin component ',ispden
         if (nspden==4) write(msg,'(3a)')    trim(msg0),' - Component ',dspin2(ispden+2*(nspden/4))
         if (dtset%natom>1.or.nspden>1) then
           call wrtout(unitfi,msg,'COLL')
         end if
         call pawio_print_ij(unitfi,pawrhoij_all(iatom)%rhoijp(:,ispden),pawrhoij_all(iatom)%nrhoijsel,&
&         pawrhoij_all(iatom)%cplex,pawrhoij_all(iatom)%lmn_size,-1,idum,1,dtset%pawprtvol,&
&         pawrhoij_all(iatom)%rhoijselect(:),valmx,1,opt_sym=optsym)
       end do
     end do
     msg=' '
     call wrtout(unitfi,msg,'COLL')
   end do
 end if

!PAW+U or local exact-exchange: print out +U components of occupancies
!-------------------------------------------------------------------------------
 if ((usepawu.or.useexexch).and.ipositron/=1.and.me_atom==0) then
   do i_unitfi=1,2
     unitfi=ab_out;if (i_unitfi==2) unitfi=std_out
     if(useexexch) write(msg,'(a)') &
&     ' "Local exact-exchange" part of augmentation waves occupancies Rhoij:'
     if(usepawu) write(msg,'(a)') &
&     ' "PAW+U" part of augmentation waves occupancies Rhoij:'
     call wrtout(unitfi,msg,'COLL')
     valmx=-1._dp
     do iatom=1,dtset%natom
       nspden=pawrhoij_all(iatom)%nspden;itypat=dtset%typat(iatom)
       ll=-1;llp=-1
       if (pawtab(itypat)%usepawu>0) ll=pawtab(itypat)%lpawu
       if (pawtab(itypat)%useexexch>0) llp=pawtab(itypat)%lexexch
       if (ll/=llp.and.ll/=-1.and.llp/=-1) then
         MSG_BUG(" lpawu/=lexexch forbidden!")
       end if
       ABI_ALLOCATE(opt_l_index, (pawtab(itypat)%lmn_size))
       ll=max(ll,llp)
       if (ll>=0) then
         optsym=2;if (pawrhoij_all(iatom)%cplex==2.and.dtset%nspinor==1) optsym=1
         do ispden=1,nspden
           msg='' ; msg0=''
           write(msg0,'(a,i3,a,i1,a)') ' Atom #',iatom,' - L=',ll,' ONLY'
           if (nspden==1) write(msg,'(a)')     trim(msg0)
           if (nspden==2) write(msg,'(2a,i1)') trim(msg0),' - Spin component ',ispden
           if (nspden==4) write(msg,'(3a)')    trim(msg0),' - Component ',dspin2(ispden+2*(nspden/4))
           call wrtout(unitfi,msg,'COLL')

           opt_l_index = pawtab(itypat)%indlmn(1,1:pawtab(itypat)%lmn_size)
           call pawio_print_ij(unitfi,pawrhoij_all(iatom)%rhoijp(:,ispden),pawrhoij_all(iatom)%nrhoijsel,&
&           pawrhoij_all(iatom)%cplex,pawrhoij_all(iatom)%lmn_size,ll,&
&           opt_l_index,1,dtset%pawprtvol,&
&           pawrhoij_all(iatom)%rhoijselect(:),valmx,1,opt_sym=optsym)

         end do
         if (debug.and.paw_ij_all(iatom)%ndij==4) then
           sz1=paw_ij_all(iatom)%lmn2_size*paw_ij_all(iatom)%cplex_dij
           sz2=paw_ij_all(iatom)%ndij
           ABI_ALLOCATE(rhoijs,(sz1,sz2))
           do irhoij=1,pawrhoij_all(iatom)%nrhoijsel
             jrhoij=paw_ij_all(iatom)%cplex_dij*(irhoij-1)+1
             rhoijs(jrhoij,1)=pawrhoij_all(iatom)%rhoijp(jrhoij,1)+pawrhoij_all(iatom)%rhoijp(jrhoij,4)
             rhoijs(jrhoij+1,1)=pawrhoij_all(iatom)%rhoijp(jrhoij+1,1)+pawrhoij_all(iatom)%rhoijp(jrhoij+1,4)
             rhoijs(jrhoij,2)=pawrhoij_all(iatom)%rhoijp(jrhoij,1)-pawrhoij_all(iatom)%rhoijp(jrhoij,4)
             rhoijs(jrhoij+1,2)=pawrhoij_all(iatom)%rhoijp(jrhoij+1,1)-pawrhoij_all(iatom)%rhoijp(jrhoij+1,4)
             rhoijs(jrhoij,3)=pawrhoij_all(iatom)%rhoijp(jrhoij,2)+pawrhoij_all(iatom)%rhoijp(jrhoij+1,3)
             rhoijs(jrhoij+1,3)=pawrhoij_all(iatom)%rhoijp(jrhoij+1,2)-pawrhoij_all(iatom)%rhoijp(jrhoij,3)
             rhoijs(jrhoij,4)=pawrhoij_all(iatom)%rhoijp(jrhoij ,2)-pawrhoij_all(iatom)%rhoijp(jrhoij+1,3)
             rhoijs(jrhoij+1,4)=pawrhoij_all(iatom)%rhoijp(jrhoij+1,2)+pawrhoij_all(iatom)%rhoijp(jrhoij,3)
           end do
           do ispden=1,nspden
             write(msg,'(3a)') trim(msg0),' - Component ',dspin1(ispden+2*(nspden/4))
             opt_l_index = pawtab(itypat)%indlmn(1,1:pawtab(itypat)%lmn_size)

             call pawio_print_ij(unitfi,rhoijs(:,ispden),pawrhoij_all(iatom)%nrhoijsel,&
&             pawrhoij_all(iatom)%cplex,pawrhoij_all(iatom)%lmn_size,ll,&
&             opt_l_index,1,dtset%pawprtvol,&
&             pawrhoij_all(iatom)%rhoijselect(:),valmx,1,opt_sym=optsym)

             call wrtout(unitfi,"WARNING: half of the array is not correct",'COLL')
           end do
           ABI_DEALLOCATE(rhoijs)
         end if
       end if
       ABI_DEALLOCATE(opt_l_index)
     end do ! iatom
     call wrtout(unitfi,' ','COLL')
   end do
 end if

!PAW+U: print out occupations for correlated orbitals
!----------------------------------------------------
 if (usepawu.and.ipositron/=1.and.me_atom==0) then
   do i_unitfi=1,2
     unitfi=ab_out;if (i_unitfi==2) unitfi=std_out
     write(msg,'(3a)') &
&     ' ---------- LDA+U DATA --------------------------------------------------- ',ch10
     call wrtout(unitfi,msg,'COLL')
     do iatom=1,dtset%natom
       itypat=dtset%typat(iatom);ll=pawtab(itypat)%lpawu
       nspden=paw_ij_all(iatom)%nspden
       if ((ll>=0).and.(pawtab(itypat)%usepawu>0)) then
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
         if(nspden==4) then
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
         do ispden=1,nspden
           if (nspden==1.and.(paw_ij_all(iatom)%cplex_dij==1)) write(msg,fmt='(a)') " Up component only..."
           if (nspden==2) write(msg,fmt='(a,i3)')" Occupation matrix for spin",ispden
           if (nspden==4.or.(paw_ij_all(iatom)%cplex_dij==2)) &
&           write(msg,fmt='(2a)')  " Occupation matrix for component ",trim(dspin1(ispden+2*(nspden/4)))
           call wrtout(unitfi,msg,'COLL')
           do im1=1,ll*2+1
             if(paw_ij_all(iatom)%cplex_dij==1)&
&             write(msg,'(12(1x,9(1x,f10.5)))') (paw_ij_all(iatom)%noccmmp(1,im1,im2,ispden),im2=1,ll*2+1)
             if(paw_ij_all(iatom)%cplex_dij==2)&
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
       end if ! ((ll>=0).and.(pawtab(itypat)%usepawu>0))
     end do
   end do
 end if

!Exact exchange: print out occupations for correlated orbitals
!-------------------------------------------------------------
 if (useexexch.and.ipositron/=1.and.me_atom==0) then
   nspden=paw_ij_all(1)%nspden;nsppol=paw_ij_all(1)%nsppol
   do iatom=1,dtset%natom
     itypat=dtset%typat(iatom);ll=pawtab(itypat)%lexexch
     if (ll>=0.and.pawtab(itypat)%useexexch>0) then
       ABI_ALLOCATE(paw_ij_all(iatom)%noccmmp,(paw_ij_all(iatom)%cplex_dij,2*ll+1,2*ll+1,nspden))
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
       if ((ll>=0).and.(pawtab(itypat)%useexexch>0)) then
         write(msg,fmt='(a,i5,a,i4,a)') " ====== For Atom",iatom,&
&         ", occupations for correlated orbitals. l =",ll,ch10
         call wrtout(unitfi,msg,'COLL')
         do ispden=1,nspden
           if (nspden==1) write(msg,fmt='(a)')   " Up component only..."
           if (nspden==2) write(msg,fmt='(a,i3)')" Occupation matrix for spin",ispden
           if (nspden==4) write(msg,fmt='(2a)')  " Occupation matrix for component ",&
&           trim(dspin2(ispden+2*(nspden/4)))
           call wrtout(unitfi,msg,'COLL')
           do im1=1,ll*2+1
             if(paw_ij_all(iatom)%cplex_dij==1)&
&             write(msg,'(12(1x,9(1x,f10.5)))')&
&             (paw_ij_all(iatom)%noccmmp(1,im1,im2,ispden),im2=1,ll*2+1)
             if(paw_ij_all(iatom)%cplex_dij==2)&
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
