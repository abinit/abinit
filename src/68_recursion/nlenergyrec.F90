!{\src2tex{textfont=tt}}
!!****f* ABINIT/nlenergyrec
!! NAME
!! nlenergyrec
!!
!! FUNCTION
!! During recursion, it computes the non-local energy
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2018 ABINIT group (the_author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  rset<recursion_type>=contains all recursion parameters
!!  exppot=exponential of -1/tsmear*vtrial (computed only once in vtorhorec)
!!  tsmear=temperature (Hartree)
!!  trotter=trotter parameter
!!  tol=tolerance criteria for stopping recursion_nl
!!  ngfft=information about FFT(dtset%ngfft a priori different from ngfftrec)
!!  mpi_enreg=information about MPI paralelisation
!!  rset<recursion_type> contains all parameter of recursion 
!!  typat(natom)=type of pseudo potential associated to any atom
!!  natom=number of atoms
!!
!! OUTPUT
!!  enl=non-local energy 
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      vtorhorec
!!
!! CHILDREN
!!      recursion_nl,reshape_pot,timab,wrtout,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine nlenergyrec(rset,enl,exppot,ngfft,natom,typat,&
 &                      tsmear,trotter,tol)

 use defs_basis
 use defs_rectypes
 use m_profiling_abi
 use m_xmpi
 use m_per_cond

 use m_time,       only : timab
 use m_rec_tools,  only : reshape_pot

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'nlenergyrec'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_68_recursion, except_this_one => nlenergyrec
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!Scalar
 integer , intent(in)  :: natom,trotter
 real(dp), intent(in)  :: tsmear,tol
 type(recursion_type),intent(in) :: rset
 real(dp), intent(out) :: enl
!Arrays
 integer , intent(in)  :: typat(natom),ngfft(18)
 real(dp), intent(in)  :: exppot(0:ngfft(1)*ngfft(2)*ngfft(3)-1)
!Local variables-------------------------------
 integer :: iatom,jatom
 integer :: ii,ipsp,dim_trott
 integer :: ierr,me_count
 integer :: ilmn,jlmn,ilm,jlm,in,jn,il
 character(len=500) :: msg                   
 logical  :: tronc
 real(dp) :: rho_nl,normali,mult
 type(mpi_type):: mpi_loc
!Arrays
 integer  :: gcart_loc(3,natom)
 integer  :: ngfftrec(3),trasl(3)
 real(dp) :: tsec(2)
 real(dp) :: un0(0:rset%nfftrec)
 real(dp),pointer :: projec(:,:,:,:,:)
 real(dp),allocatable ::  exppotloc(:)
 real(dp) :: proj_arr(0:rset%ngfftrec(1)-1,0:rset%ngfftrec(2)-1,0:rset%ngfftrec(3)-1)

! *************************************************************************
 

 call timab(612,1,tsec) !!--start time-counter: nlenergyrec

 if(rset%debug)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' nlenergyrec : enter'
   call wrtout(std_out,msg,'PERS')
 end if
 
 write(msg,'(a)')' -- nlenergyrec -----------------------------------'
 call wrtout(std_out,msg,'COLL')
 
!--Initialisation variables
 enl = zero
 mult = two !--is twice for non-spinned systems
 ngfftrec = rset%ngfftrec(:3)    
 gcart_loc = rset%inf%gcart
 mpi_loc = rset%mpi
 me_count = 0
 dim_trott = max(0,2*trotter-1)
 nullify(projec)
 ABI_ALLOCATE(projec,(0:rset%ngfftrec(1)-1,0:rset%ngfftrec(2)-1,0:rset%ngfftrec(3)-1,rset%nl%lmnmax,natom))
 projec = zero

 tronc = rset%tronc  !--True if troncation is used
 if(tronc)   then
   ABI_ALLOCATE(exppotloc,(0:rset%nfftrec-1))
 end if


!--LOOP ON ATOMS to create projectors-vector
 atomloop1: do iatom = 1, natom
   ipsp = typat(iatom)
!  --Aquisition,reshape,translation,rotation of the projectors vector
   do ilmn = 1,rset%nl%lmnmax  
     in = rset%nl%indlmn(3,ilmn,ipsp) 
!    --Projectors vector in 3-composant vector
     projec(:,:,:,ilmn,iatom) = reshape(rset%nl%projec(:,ilmn,ipsp),shape=shape(projec(:,:,:,1,1))) 
!    --Moving the projectors vector on the center of the grid
     do ii=1,3
       projec(:,:,:,ilmn,iatom) = cshift(projec(:,:,:,ilmn,iatom),shift=ngfftrec(ii)/2-gcart_loc(ii,iatom),dim=ii)
     end do
   end do

 end do atomloop1


!##################################################################
!--LOOP ON ATOMS (MAIN LOOP)
 atomloop: do iatom = 1, natom
   ipsp = typat(iatom)

!  --If troncation is present, the considered atom has to be in the
!  center of the grid so atoms, potential and projectors have to be translated
   if(tronc) then
     trasl = -rset%inf%gcart(:,iatom)+ngfftrec/2
!    --Translation of atoms
     do jatom=1,natom 
       gcart_loc(:,jatom) = rset%inf%gcart(:,jatom)+trasl
       gcart_loc(:,jatom) = modulo(gcart_loc(:,jatom),ngfft(:3))
!      --Translation of non-local projectors
       do ilmn = 1,rset%nl%lmnmax
         projec(:,:,:,ilmn,jatom) = reshape(rset%nl%projec(:,ilmn,typat(jatom)),shape=shape(projec(:,:,:,1,1)))
         do ii=1,3
           projec(:,:,:,ilmn,jatom) = eoshift(projec(:,:,:,ilmn,jatom),shift=ngfftrec(ii)/2-gcart_loc(ii,jatom),dim=ii)
         end do
       end do
     end do

!    --Translation of the potential
     call reshape_pot(trasl,ngfft(1)*ngfft(2)*ngfft(3),rset%nfftrec,ngfft(:3),ngfftrec,exppot,exppotloc)
   end if
   
!  --Loop on projectors
   projloop: do ilmn = 1,rset%nl%lmnmax
     me_count = iatom+ilmn*natom-2 !--counter of the number of iteration
!    --Only the proc me compute
     if(mpi_loc%me==mod(me_count,mpi_loc%nproc)) then
       ilm = rset%nl%indlmn(4,ilmn,ipsp)
       proj_arr = zero 
       do jlmn = 1,rset%nl%lmnmax  
         jlm = rset%nl%indlmn(4,jlmn,ipsp)
         if(ilm==jlm) then
           in = rset%nl%indlmn(3,ilmn,ipsp) 
           jn = rset%nl%indlmn(3,jlmn,ipsp) 
           il = rset%nl%indlmn(1,ilmn,ipsp)+1 
           proj_arr(:,:,:) = proj_arr(:,:,:) + rset%nl%eivec(jn,in,il,ipsp)*projec(:,:,:,jlmn,iatom)
!          write(std_out,*)'l,m,lm,n,n',il-1,rset%nl%indlmn(2,ilmn,ipsp),ilm,in,jn
!          write(std_out,*)'eigevectors',rset%nl%eivec(jn,in,il,ipsp)

         end if
       end do

       un0 = pack(proj_arr(:,:,:),mask=.true.)
       normali = sum(un0*un0)*rset%inf%ucvol
       un0 = (one/sqrt(normali))*un0

       if(tronc)then
         call recursion_nl(exppotloc,un0,rho_nl,rset,rset%ngfftrec,&
&         tsmear,trotter,dim_trott,tol,typat,&
&         natom,projec)
       else
         call recursion_nl(exppot,un0,rho_nl,rset,rset%ngfftrec,&
&         tsmear,trotter,dim_trott,tol,typat,&
&         natom,projec)
       end if

       enl = enl+mult*rho_nl*rset%nl%eival(in,il,ipsp)*normali
     end if

   end do projloop
 end do atomloop

!--Sum the contribution to the non-local energy computed by any procs
 call xmpi_sum(enl,mpi_loc%comm_bandfft,ierr)
 
 if(associated(projec))  then
   ABI_DEALLOCATE(projec)
 end if
 if(tronc)  then
   ABI_DEALLOCATE(exppotloc)
 end if

 if(rset%debug)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,' nlenergyrec : exit'
   call wrtout(std_out,msg,'PERS')
 end if

 call timab(612,2,tsec)  !--stop  time-counter: nlenergyrec
end subroutine nlenergyrec
!!***
