!{\src2tex{textfont=tt}}
!!****f* ABINIT/denfgr
!! NAME
!! denfgr
!!
!! FUNCTION
!!  Construct complete electron density on fine grid, by removing nhat
!!  and adding PAW corrections
!!
!! COPYRIGHT
!!   Copyright (C) 2005-2016 ABINIT group (JWZ)
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~ABINIT/COPYING
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!   atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!   gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!   mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!   comm_atom=--optional-- MPI communicator over atoms
!!   my_natom=number of atoms treated by current processor
!!   natom= number of atoms in cell
!!   nattyp(ntypat)= # atoms of each type.
!!   ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!   nhat(pawfgr%nfft,nspden)= compensation charge density used in PAW
!!   nspinor=Number of spinor components
!!   nsppol=Number of independent spin components.
!!   nspden= number of spin densities
!!   ntypat= number of types of atoms in the cell
!!   pawfgr <type(pawfgr_type)>= data about the fine grid
!!   pawrad(ntypat) <type(pawrad_type)>= radial mesh data for each type of atom
!!   pawrhoij(natom) <type(pawrhoij_type)>= rho_ij data for each atom
!!   pawtab(ntypat) <type(pawtab_type)>= PAW functions around each type of atom
!!   rhor(pawfgr%nfft,nspden)= input density ($\tilde{n}+\hat{n}$ in PAW case)
!!   rprimd(3,3)=dimensional primitive translations for real space (bohr)
!!   typat(natom)= list of atom types
!!   ucvol=unit cell volume (bohr**3)
!!   xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!! rhor_paw(pawfgr%nfft,nspden)= full electron density on the fine grid
!!
!! NOTES
!!   In PAW calculations, the valence density present in rhor includes the
!!   compensation charge density $\hat{n}$, and also doesn't include the on-site
!!   PAW contributions. For post-processing and proper visualization it is necessary
!!   to use the full electronic density, which is what this subroutine constructs.
!!   Specifically, it removes $\hat{n}$ from rhor, and also computes the on-site PAW
!!   terms. This is nothing other than the proper PAW treatment of the density
!!   operator $|\mathbf{r}\rangle\langle\mathbf{r}|$, and yields the formula
!!   $$\tilde{n}+\sum_{ij}\rho_ij\left[\varphi_i(\mathbf{r})\varphi_j(\mathbf{r})-
!!   \tilde{\varphi}_i(\mathbf{r})\tilde{\varphi}_j(\mathbf{r})\right]$$
!!   Notice that this formula is expressed on the fine grid, and requires
!!   interpolating the PAW radial functions onto this grid, as well as calling
!!   initylmr in order to get the angular functions on the grid points.
!!
!! PARENTS
!!      bethe_salpeter,outscfcv,sigma
!!
!! CHILDREN
!!      free_my_atmtab,get_my_atmtab,initylmr,nhatgrid,pawfgrtab_free
!!      pawfgrtab_init,pawrad_deducer0,pawtab_get_lsize,printxsf,sort_dp,spline
!!      splint,wrtout,xmpi_barrier,xmpi_sum,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine denfgr(atindx1,gmet,spaceComm_in,my_natom,natom,nattyp,ngfft,nhat,nspinor,nsppol,nspden,ntypat, &
& pawfgr,pawrad,pawrhoij,pawtab,prtvol,rhor,rhor_paw,rhor_n_one,rhor_nt_one,rprimd,typat,ucvol,xred,&
& abs_n_tilde_nt_diff,znucl,mpi_atmtab,comm_atom) ! Optional arguments

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi
 use m_splines

 use m_io_tools,    only : open_file
 use m_pptools,     only : printxsf
 use m_pawrad,      only : pawrad_type, pawrad_deducer0
 use m_pawtab,      only : pawtab_type,pawtab_get_lsize
 use m_pawfgrtab,   only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free
 use m_pawrhoij,    only : pawrhoij_type
 use m_pawfgr,      only : pawfgr_type
 use m_paw_sphharm, only : initylmr
 use m_paral_atom,  only : get_my_atmtab, free_my_atmtab

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'denfgr'
 use interfaces_14_hidewrite
 use interfaces_28_numeric_noabirule
 use interfaces_41_geometry
 use interfaces_65_paw, except_this_one => denfgr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom,natom,nspden,ntypat,prtvol,nsppol,nspinor
 integer,optional,intent(in) :: comm_atom
 real(dp),intent(in) :: ucvol
 type(pawfgr_type),intent(in) :: pawfgr
!arrays
 integer,intent(in) :: spaceComm_in
 integer,intent(in) :: atindx1(natom),nattyp(ntypat),ngfft(18),typat(natom)
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: gmet(3,3),nhat(pawfgr%nfft,nspden)
 real(dp),intent(in) :: rhor(pawfgr%nfft,nspden),rprimd(3,3)
 real(dp),intent(inout) :: xred(3,natom)
 real(dp),intent(out) :: rhor_paw(pawfgr%nfft,nspden)
 real(dp),intent(out) :: rhor_n_one(pawfgr%nfft,nspden)
 real(dp),intent(out) :: rhor_nt_one(pawfgr%nfft,nspden)
 real(dp),optional,intent(out) :: abs_n_tilde_nt_diff(nspden)
 real(dp),optional,intent(in) :: znucl(ntypat)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom)
 type(pawtab_type),target,intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: delta,iatom,ierr,ifgd,ifftsph,inl,inrm,ipsang,irhoij
 integer :: ispden,itypat,il,im,ilm,iln,ilmn
 integer :: jl,jlm,jln,jm,j0lmn,jlmn
 integer :: klmn,my_comm_atom,my_start_indx,my_end_indx
 integer :: nfgd,nnl,normchoice,nprocs,optcut,optgr0,optgr1,optgr2
 integer :: optrad,option,my_rank,remainder,tmp_unt
 real(dp) :: phj,phi,rR,tphj,tphi,ybcbeg,ybcend
 logical :: my_atmtab_allocated,paral_atom
 character(len=500) :: message
!arrays
 integer,allocatable :: l_size_atm(:),nrm_ifftsph(:)
 integer,ABI_CONTIGUOUS pointer :: indlmn(:,:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: ylmgr(3,3,0)
 real(dp) :: yvals(4),xcart(3,natom)
 real(dp),allocatable :: diag(:),nrm(:),phigrd(:,:),tphigrd(:,:),ylm(:,:),ypp(:)
 real(dp),allocatable :: phi_at_zero(:),tphi_at_zero(:)
 real(dp),allocatable :: rhor_tmp(:,:),tot_rhor(:)
 character(len=fnlen) :: xsf_fname
 type(pawfgrtab_type) :: local_pawfgrtab(my_natom)

! ************************************************************************

 DBG_ENTER("COLL")

!Set up parallelism over atoms (compatible only with band-FFT parallelism)
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,&
& my_natom_ref=my_natom)

!MG: FIXME It won't work if atom-parallelism is used
!but even the loop over atoms below should be rewritten in this case.

!use a local copy of pawfgrtab to make sure we use the correction in the paw spheres
!the usual pawfgrtab uses r_shape which may not be the same as r_paw
 if (my_natom>0) then
   if (paral_atom) then
     call pawtab_get_lsize(pawtab,l_size_atm,my_natom,typat,mpi_atmtab=my_atmtab)
     call pawfgrtab_init(local_pawfgrtab,pawrhoij(1)%cplex,l_size_atm,nspden,typat,&
&     mpi_atmtab=my_atmtab,comm_atom=my_comm_atom)
   else
     call pawtab_get_lsize(pawtab,l_size_atm,my_natom,typat)
     call pawfgrtab_init(local_pawfgrtab,pawrhoij(1)%cplex,l_size_atm,nspden,typat)
   end if
   ABI_DEALLOCATE(l_size_atm)
 end if

!Note: call to nhatgrid: comm_fft not used because FFT parallelism
!is done manually below
 optcut = 1 ! use rpaw to construct local_pawfgrtab
 optgr0 = 0; optgr1 = 0; optgr2 = 0 ! dont need gY terms locally
 optrad = 1 ! do store r-R
 if (paral_atom) then
   call nhatgrid(atindx1,gmet,my_natom,natom,nattyp,ngfft,ntypat,&
&   optcut,optgr0,optgr1,optgr2,optrad,local_pawfgrtab,pawtab,rprimd,typat,ucvol,xred,&
&   comm_atom=my_comm_atom,mpi_atmtab=my_atmtab)
 else
   call nhatgrid(atindx1,gmet,my_natom,natom,nattyp,ngfft,ntypat,&
&   optcut,optgr0,optgr1,optgr2,optrad,local_pawfgrtab,pawtab,rprimd,typat,ucvol,xred)
 end if
!now local_pawfgrtab is ready to use

!Initialise output arrays.
 rhor_paw=zero; rhor_n_one=zero; rhor_nt_one=zero

!Initialise and check parallell execution
 my_rank = xmpi_comm_rank(spaceComm_in)
 nprocs = xmpi_comm_size(spaceComm_in)


!loop over atoms in cell
 do iatom = 1, my_natom
   itypat = pawrhoij(iatom)%itypat
   indlmn => pawtab(itypat)%indlmn
   nfgd = local_pawfgrtab(iatom)%nfgd ! number of points in the fine grid for this PAW sphere
   nnl = pawtab(itypat)%basis_size ! number of nl elements in PAW basis

!  Division of fine grid points among processors
   if (nprocs==1) then ! Make sure everything runs with one proc
     write(message,'(a)') '  In denfgr - number of processors:     1'
     call wrtout(std_out,message,'COLL')
     write(message,'(a)') '  Calculation of PAW density done in serial'
     call wrtout(std_out,message,'COLL')
     write(message,'(a,I6)') '  Number of fine grid points:',nfgd
     call wrtout(std_out,message,'COLL')
     my_start_indx = 1
     my_end_indx = nfgd
   else ! Divide up the fine grid points among the processors
     write(message,'(a,I4)') '  In denfgr - number of processors: ',nprocs
     call wrtout(std_out,message,'COLL')
     write(message,'(a)') '  Calculation of PAW density done in parallel'
     call wrtout(std_out,message,'COLL')
     write(message,'(a,I6)') '  Number of fine grid points:',nfgd
     call wrtout(std_out,message,'COLL')
!    Divide the fine grid points among the processors
     delta = int(floor(real(nfgd)/real(nprocs)))
     remainder = nfgd-nprocs*delta
     my_start_indx = 1+my_rank*delta
     my_end_indx = (my_rank+1)*delta
!    Divide the remainder points among the processors
!    by shuffling indices
     if ((my_rank+1)>remainder) then
       my_start_indx = my_start_indx + remainder
       my_end_indx = my_end_indx + remainder
     else
       my_start_indx = my_start_indx + my_rank
       my_end_indx = my_end_indx + my_rank + 1
     end if
     if (prtvol>9) then
       write(message,'(a,I6)') '  My index Starts at: ',my_start_indx
       call wrtout(std_out,message,'PERS')
       write(message,'(a,I6)') '             Ends at: ',my_end_indx
       call wrtout(std_out,message,'PERS')
       write(message,'(a,I6)') '               # pts: ',my_end_indx+1-my_start_indx
       call wrtout(std_out,message,'PERS')
     end if
   end if

   write(message,'(a,I3,a,I3)') '  Entered loop for atom: ',iatom,' of:',natom
   call wrtout(std_out,message,'PERS')

!  obtain |r-R| values on fine grid
   ABI_ALLOCATE(nrm,(nfgd))
   do ifgd=1, nfgd
     nrm(ifgd) = sqrt(dot_product(local_pawfgrtab(iatom)%rfgd(:,ifgd),local_pawfgrtab(iatom)%rfgd(:,ifgd)))
   end do ! these are the |r-R| values

!  compute Ylm for each r-R vector.
!  ----
   ipsang = 1 + (pawtab(itypat)%l_size - 1)/2 ! recall l_size=2*l_max+1
   ABI_ALLOCATE(ylm,(ipsang*ipsang,nfgd))
   option = 1 ! compute Ylm(r-R) for vectors
   normchoice = 1 ! use computed norms of input vectors
   call initylmr(ipsang,normchoice,nfgd,nrm,option,local_pawfgrtab(iatom)%rfgd,ylm,ylmgr)

!  in order to do spline fits, the |r-R| data must be sorted
!  ----
   ABI_ALLOCATE(nrm_ifftsph,(nfgd))
   nrm_ifftsph(:) = local_pawfgrtab(iatom)%ifftsph(:) ! copy of indices of points, to be rearranged by sort_dp
   call sort_dp(nfgd,nrm,nrm_ifftsph,tol8) ! sort the nrm points, keeping track of which goes where

!  now make spline fits of phi and tphi  onto the fine grid around the atom
!  ----
   ABI_ALLOCATE(phigrd,(nfgd,nnl))
   ABI_ALLOCATE(tphigrd,(nfgd,nnl))
   ABI_ALLOCATE(phi_at_zero,(nnl))
   ABI_ALLOCATE(tphi_at_zero,(nnl))
   ABI_ALLOCATE(ypp,(pawtab(itypat)%mesh_size))
   ABI_ALLOCATE(diag,(pawtab(itypat)%mesh_size))

   do inl = 1, nnl

!    spline phi onto points
     ypp(:) = zero; diag(:) = zero; ybcbeg = zero; ybcend = zero;
     call spline(pawrad(itypat)%rad,pawtab(itypat)%phi(:,inl),pawtab(itypat)%mesh_size,ybcbeg,ybcend,ypp)
     call splint(pawtab(itypat)%mesh_size,pawrad(itypat)%rad,pawtab(itypat)%phi(:,inl),ypp,nfgd,nrm,phigrd(:,inl))

!    next splint tphi onto points
     ypp(:) = zero; diag(:) = zero; ybcbeg = zero; ybcend = zero;
     call spline(pawrad(itypat)%rad,pawtab(itypat)%tphi(:,inl),pawtab(itypat)%mesh_size,ybcbeg,ybcend,ypp)
     call splint(pawtab(itypat)%mesh_size,pawrad(itypat)%rad,pawtab(itypat)%tphi(:,inl),ypp,nfgd,nrm,tphigrd(:,inl))

!    Find out the value of the basis function at zero using extrapolation
     yvals = zero
!    Extrapolate only if this is an s-state (l=0)
     if (indlmn(1,inl)==0) then
       yvals(2:4) = pawtab(itypat)%phi(2:4,inl)/pawrad(itypat)%rad(2:4)
       call pawrad_deducer0(yvals,4,pawrad(itypat))
       write(std_out,*) 'phi_at_zero: ',yvals(1),' from:',yvals(2:4)
     end if
     phi_at_zero(inl) = yvals(1)

     yvals = zero
!    Extrapolate only if this is an s-state (l=0)
     if (indlmn(1,inl)==0) then
       yvals(2:4) = pawtab(itypat)%tphi(2:4,inl)/pawrad(itypat)%rad(2:4)
       call pawrad_deducer0(yvals,4,pawrad(itypat))
       write(std_out,*) 'tphi_at_zero: ',yvals(1),' from:',yvals(2:4)
     end if
     tphi_at_zero(inl) = yvals(1)

   end do ! end loop over nnl basis functions
   ABI_DEALLOCATE(ypp)
   ABI_DEALLOCATE(diag)

!  loop over basis elements for this atom
!  because we have to store things like <phi|r'><r'|phi>-<tphi|r'><r'|tphi> at each point of the
!  fine grid, there is no integration, and hence no simplifications of the Y_lm's. That's why
!  we have to loop through the basis elements in exhaustive detail, rather than just a loop over
!  lmn2_size or something comparable.
!  ----
   if (prtvol>9) then
     write(message,'(a,I3)') '  Entering j-loop over basis elements for atom:',iatom
     call wrtout(std_out,message,'PERS')
   end if

   do jlmn=1,pawtab(itypat)%lmn_size

     if (prtvol>9) then
       write(message,'(2(a,I3))') '  Element:',jlmn,' of:',pawtab(itypat)%lmn_size
       call wrtout(std_out,message,'PERS')
     end if

     jl=indlmn(1,jlmn)
     jm=indlmn(2,jlmn)
     jlm=indlmn(4,jlmn)
     jln=indlmn(5,jlmn)
     j0lmn=jlmn*(jlmn-1)/2

     if (prtvol>9) then
       write(message,'(a,I3)') '  Entering i-loop for j:',jlmn
       call wrtout(std_out,message,'PERS')
     end if

     do ilmn=1,jlmn

       if (prtvol>9) then
         write(message,'(2(a,I3))') '    Element:',ilmn,' of:',jlmn
         call wrtout(std_out,message,'PERS')
       end if

       il=indlmn(1,ilmn)
       im=indlmn(2,ilmn)
       iln=indlmn(5,ilmn)
       ilm=indlmn(4,ilmn)
       klmn=j0lmn+ilmn

       if (prtvol>9) then
         write(message,'(a)') '    Entering loop over nonzero elems of rhoij'
         call wrtout(std_out,message,'PERS')
       end if

!      Loop over non-zero elements of rhoij
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         if (klmn==pawrhoij(iatom)%rhoijselect(irhoij)) then ! rho_ij /= 0 for this klmn

           do ifgd=my_start_indx, my_end_indx ! loop over fine grid points in current PAW sphere
             ifftsph = local_pawfgrtab(iatom)%ifftsph(ifgd) ! index of the point on the grid

!            have to retrieve the spline point to use since these were sorted
             do inrm=1, nfgd
               if(nrm_ifftsph(inrm) == ifftsph) exit ! have found nrm point corresponding to nfgd point
             end do ! now inrm is the index of the sorted nrm vector to use

!            avoid division by zero
             if(nrm(inrm) > zero) then
               rR = nrm(inrm) ! value of |r-R| in the following
!              recall that <r|phi>=u(r)*Slm(r^)/r
               phj  = phigrd(inrm,jln)*ylm(jlm,ifgd)/rR
               phi  = phigrd(inrm,iln)*ylm(ilm,ifgd)/rR
               tphj = tphigrd(inrm,jln)*ylm(jlm,ifgd)/rR
               tphi = tphigrd(inrm,iln)*ylm(ilm,ifgd)/rR
             else
!              use precalculated <r|phi>=u(r)*Slm(r^)/r at r=0
               phj  = phi_at_zero(jln)*ylm(jlm,ifgd)
               phi  = phi_at_zero(iln)*ylm(ilm,ifgd)
               tphj = tphi_at_zero(jln)*ylm(jlm,ifgd)
               tphi = tphi_at_zero(iln)*ylm(ilm,ifgd)
             end if ! check if |r-R| = 0

             do ispden=1,nspden
               if (pawrhoij(iatom)%cplex == 1) then
                 rhor_paw(ifftsph,ispden) = rhor_paw(ifftsph,ispden) + &
&                 pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(irhoij,ispden)*(phj*phi - tphj*tphi)

                 rhor_n_one(ifftsph,ispden) = rhor_n_one(ifftsph,ispden) + &
&                 pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(irhoij,ispden)*phj*phi

                 rhor_nt_one(ifftsph,ispden) = rhor_nt_one(ifftsph,ispden) + &
&                 pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(irhoij,ispden)*tphj*tphi
               else
                 rhor_paw(ifftsph,ispden) = rhor_paw(ifftsph,ispden) + &
&                 pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(2*irhoij-1,ispden)*(phj*phi - tphj*tphi)

                 rhor_n_one(ifftsph,ispden) = rhor_n_one(ifftsph,ispden) + &
&                 pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(2*irhoij-1,ispden)*phj*phi

                 rhor_nt_one(ifftsph,ispden) = rhor_nt_one(ifftsph,ispden) + &
&                 pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(2*irhoij-1,ispden)*tphj*tphi
               end if ! end check on cplex rhoij

             end do ! end loop over nsdpen
           end do ! end loop over nfgd
         end if ! end selection on rhoij /= 0
       end do ! end loop over non-zero rhoij
     end do ! end loop over ilmn atomic basis states
   end do ! end loop over jlmn atomic basis states

   ABI_DEALLOCATE(nrm)
   ABI_DEALLOCATE(nrm_ifftsph)
   ABI_DEALLOCATE(phigrd)
   ABI_DEALLOCATE(tphigrd)
   ABI_DEALLOCATE(ylm)
   ABI_DEALLOCATE(phi_at_zero)
   ABI_DEALLOCATE(tphi_at_zero)
 end do     ! Loop on atoms

!MPI sum on each node the different contributions to the PAW densities.
 call xmpi_sum(rhor_paw,spaceComm_in,ierr)
 call xmpi_sum(rhor_n_one,spaceComm_in,ierr)
 call xmpi_sum(rhor_nt_one,spaceComm_in,ierr)
 if (paral_atom) then
   call xmpi_sum(rhor_paw,my_comm_atom,ierr)
   call xmpi_sum(rhor_n_one,my_comm_atom,ierr)
   call xmpi_sum(rhor_nt_one,my_comm_atom,ierr)
 end if

 call wrtout(std_out,' *** Partial contributions to PAW rhor summed ***','PERS')
 call xmpi_barrier(spaceComm_in)

!Add the plane-wave contribution \tilde{n} and remove \hat{n}
!BE careful here since the storage mode of rhoij and rhor is different.
 select case (nspinor)
 case (1)
   if (nsppol==1)  then
     rhor_paw = rhor_paw + rhor - nhat
   else  ! Spin-polarised case: rhor_paw contains rhor_paw(spin_up,spin_down) but we need rhor_paw(total,spin_up)
     ABI_ALLOCATE(tot_rhor,(pawfgr%nfft))
!      
!      AE rhor
     tot_rhor(:) = SUM(rhor_paw,DIM=2)
     rhor_paw(:,2) = rhor_paw(:,1)
     rhor_paw(:,1) = tot_rhor
     rhor_paw = rhor_paw + rhor - nhat
!      
!      onsite AE rhor
     tot_rhor(:) = SUM(rhor_n_one,DIM=2)
     rhor_n_one(:,2) = rhor_n_one(:,1)
     rhor_n_one(:,1) = tot_rhor
!      
!      onsite PS rhor
     tot_rhor(:) = SUM(rhor_nt_one,DIM=2)
     rhor_nt_one(:,2) = rhor_nt_one(:,1)
     rhor_nt_one(:,1) = tot_rhor

     ABI_DEALLOCATE(tot_rhor)
   end if

 case (2)
!    * if nspden==4, rhor contains (n^11, n^22, Re[n^12], Im[n^12].
!    Storage mode for rhoij is different, See pawaccrhoij.
   MSG_ERROR("nspinor 2 not coded")
 case default
   write(message,'(a,i0)')" Wrong value for nspinor=",nspinor
   MSG_ERROR(message)
 end select

!if (prtvol>9) then ! Check normalisation
!write(message,'(a,F8.4)') '  PAWDEN - NORM OF DENSITY: ',SUM(rhor_paw(:,1))*ucvol/PRODUCT(pawfgr%ngfft(1:3))
!call wrtout(std_out,message,'COLL')
!end if

 if (present(abs_n_tilde_nt_diff).AND.present(znucl)) then
   ABI_ALLOCATE(rhor_tmp,(pawfgr%nfft,nspden))
   do ispden=1,nspden
     rhor_tmp(:,ispden) = zero
     do iatom=1,my_natom
       do ifgd=1,local_pawfgrtab(iatom)%nfgd ! loop over fine grid points in current PAW sphere
         ifftsph = local_pawfgrtab(iatom)%ifftsph(ifgd) ! index of the point on the grid
         rhor_tmp(ifftsph,ispden) = rhor(ifftsph,ispden) - nhat(ifftsph,ispden) &
&         - rhor_nt_one(ifftsph,ispden)
       end do !ifgd
     end do ! iatom
   end do ! ispden
   if (paral_atom) then
     call xmpi_sum(rhor_tmp,my_comm_atom,ierr)
   end if

   if (my_rank==master) then
     do ispden=1,nspden
!      Write to xsf file
       call xred2xcart(natom,rprimd,xcart,xred)
       write(xsf_fname,'(a,I0,a)') 'N_tilde_onsite_diff_sp',ispden,'.xsf'
       if (open_file(xsf_fname,message, unit=tmp_unt,status='unknown',form='formatted') /= 0) then
         MSG_ERROR(message)
       end if
       call printxsf(ngfft(1),ngfft(2),ngfft(3),rhor_tmp(:,ispden),rprimd,&
&       (/zero,zero,zero/),natom,ntypat,typat,xcart,znucl,tmp_unt,0)
       close(tmp_unt)
       abs_n_tilde_nt_diff(ispden) = SUM(ABS(rhor_tmp(:,ispden)))/pawfgr%nfft
       write(message,'(4(a),F16.9,2(a,I0),a)') ch10,'  Wrote xsf file with \tilde{n}-\tilde{n}^1.',ch10,&
&       '  Value of norm |\tilde{n}-\tilde{n}^1|:',&
&       abs_n_tilde_nt_diff(ispden),' spin: ',ispden,' of ',nspden,ch10
       call wrtout(std_out,message,'COLL')
     end do
   end if
   ABI_DEALLOCATE(rhor_tmp)

 else if ((present(abs_n_tilde_nt_diff).AND.(.NOT.present(znucl))) &
&   .OR.(.NOT.present(abs_n_tilde_nt_diff).AND.(present(znucl)))) then
   write(message,'(a)') ' Both abs_n_tilde_nt_diff *and* znucl must be passed',ch10,&
&   'to denfgr for |\tilde{n}-\tilde{n}^1| norm evaluation.'
   MSG_ERROR(message)
 end if

 call pawfgrtab_free(local_pawfgrtab)

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

 call xmpi_barrier(spaceComm_in)

 DBG_EXIT("COLL")

 end subroutine denfgr
!!***
