!!{\src2tex{textfont=tt}}
!!****f* ABINIT/mlwfovlp_projpaw
!! NAME
!! mlwfovlp_projpaw
!!
!! FUNCTION
!! Calculates the functions that are given to 
!! Wannier90 as an starting guess.
!! Here we project them inside the PAW spheres
!!
!! COPYRIGHT
!!  Copyright (C) 2009-2017 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  band_in(mband)= logical array which indicates the bands to be excluded from the calculation
!!  cprj(natom,nspinor*mband*mkmem*nsppol)= <p_lmn|Cnk> coefficients for each WF |Cnk>
!!                                          and each |p_lmn> non-local projector
!!  just_augmentation= flag used to indicate that we are just going
!!                     to compute augmentation part of the matrix
!!                     and we are excluding the plane wave part.
!!  mband= maximum number of bands
!!  mkmem= number of k points which can fit in memory; set to 0 if use disk
!!  natom= number of atoms in cell.
!!  nband(nkpt*nsppol)= array cointaining number of bands at each k-point and isppol
!!  nkpt=number of k points.
!!  num_bands=number of bands actually used to construct the wannier function (NOT USED IN 6.7.1 SO WAS TEMPORARILY REMOVED)
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  nwan= number of wannier fonctions (read in wannier90.win).
!!  pawrad(ntypat)= type(pawrad_type) radial information of paw objects
!!  pawtab(ntypat)= For PAW, TABulated data initialized at start
!!  proj_l(mband)= angular part of the projection function (quantum number l)
!!  proj_m(mband)= angular part of the projection function (quantum number m)
!!  proj_radial(mband)= radial part of the projection.
!!  proj_site(3,mband)= site of the projection.
!!  proj_x(3,mband)= x axis for the projection.
!!  proj_z(3,mband)= z axis for the projection.
!!  proj_zona(mband)= extension of the radial part.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)= Direct lattice vectors, Bohr units.
!!  spin = just used for nsppol>1 ; 0 both, 1 just spin up, 2 just spin down
!!  typat(natom)= atom type
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  A_paw(max_num_bands,nwan,nkpt) = A matrix containing initial guess for MLWFs 
!!                          (augmentation part of the matrix)
!!
!! SIDE EFFECTS
!!
!! NOTES
!! This routine is still under developement
!!
!! PARENTS
!!      mlwfovlp
!!
!! CHILDREN
!!      mlwfovlp_ylmfar,simp_gen,simpson_int,wrtout,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mlwfovlp_projpaw(A_paw,band_in,cprj,just_augmentation,max_num_bands,mband,mkmem,&
&mwan,natom,nband,nkpt,&
&nspinor,nsppol,ntypat,nwan,pawrad,pawtab,&
&proj_l,proj_m,proj_radial,proj_site,proj_x,proj_z,proj_zona,psps,&
&rprimd,spin,typat,xred)
    
 use defs_basis
 use defs_datatypes
 use defs_wannier90
 use m_errors
 use m_profiling_abi

 use m_numeric_tools,   only : simpson_int
 use m_pawrad,  only : pawrad_type, simp_gen
 use m_pawtab,  only : pawtab_type
 use m_pawcprj, only : pawcprj_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mlwfovlp_projpaw'
 use interfaces_14_hidewrite
 use interfaces_41_geometry
 use interfaces_67_common, except_this_one => mlwfovlp_projpaw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: max_num_bands,mband,mkmem,mwan,natom,nkpt
 integer,intent(in) :: nspinor,nsppol,ntypat,spin
 !arrays
 integer,intent(in) :: nband(nsppol*nkpt),nwan(nsppol)
 integer,intent(in) :: proj_l(mband,nsppol),proj_m(mband,nsppol),proj_radial(mband,nsppol)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in):: proj_site(3,mband,nsppol)
 real(dp),intent(in) :: proj_x(3,mband,nsppol),proj_z(3,mband,nsppol),proj_zona(mband,nsppol)
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
 complex(dpc),intent(out) :: A_paw(max_num_bands,mwan,nkpt,nsppol)
 logical,intent(in) :: band_in(mband,nsppol)
 logical,intent(in)::just_augmentation(mwan,nsppol)
 type(pawcprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat) 
 type(pseudopotential_type),intent(in) :: psps

!Local variables-------------------------------
 !local variables
 integer :: basis_size,iatom,iband,ii
 integer :: ikpt,ir,isppol,itypat,iwan,jband
 integer :: ll,lm,ln,mm,ilmn
 integer :: lmn_size,max_lmax2, mesh_size,nn
 integer :: lmax(nsppol),lmax2(nsppol)
 real(dp):: aa,int_rad2,prod_real,prod_imag
 real(dp),parameter :: dx=0.015d0,rmax=10.d0,xmin=0.d0
 real(dp):: sum,wan_lm_fac,x
 complex(dpc)::prod
 character(len=500) :: message
 !arrays
 integer :: index(mband,nkpt,nsppol)
 real(dp):: dist,norm(mwan,nsppol)
 real(dp):: proj_cart(3,mwan,nsppol),proj_site_unit(3,mwan,nsppol)
 real(dp):: xcart_unit(3,natom),xred_unit(3,natom)
 real(dp),allocatable ::aux(:),ff(:),r(:),int_rad(:),rad_int(:)
 real(dp),allocatable::ylmr_fac(:,:,:)

!no_abirules
!Tables 3.1 & 3.2, User guide
 integer,save :: orb_l_defs(-5:3)=(/2,2,1,1,1,0,1,2,3/) 
!real(dp),allocatable :: ylm(:,:)

 
! *************************************************************************
 
!DEBUG
!write (std_out,*) ' mlwfovlp_projpaw : enter'
!ENDDEBUG

!DEBUG                                           ! to be uncommented, if needed
 
 write(message, '(a,a)' )ch10,&
& '** mlwfovlp_proj:  compute in-sphere part of A_matrix'
 call wrtout(std_out,message,'COLL')

!
!Check input variables
!
 do isppol=1,nsppol
   if(spin.ne.0 .and. spin.ne.isppol) cycle
   do iwan=1,nwan(nsppol)
     if(proj_radial(iwan,isppol)<1 .or. proj_radial(iwan,isppol)>4)then
       write(message,'(a,a,a,i6)')&
&       '  proj_radial should be between 1 and 4,',ch10,&
&       '  however, proj_radial=',proj_radial(iwan,isppol)
       MSG_BUG(message)
     end if
   end do
 end do

!
!Initialize
!
 A_paw(:,:,:,:)=cmplx(0.d0,0.d0)
!
!Get index for cprj
!
 ii=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     do iband=1,nband(ikpt)
       ii=ii+1
       index(iband,ikpt,isppol)=ii
     end do
   end do
 end do
!
!obtain lmax and lmax2
!
 lmax(:)=0
 lmax2(:)=0
 do isppol=1,nsppol
   if(spin.ne.0 .and. spin.ne.isppol) cycle
   do iwan=1,nwan(isppol)
     lmax(isppol)=max(lmax(isppol),orb_l_defs(proj_l(iwan,isppol)))
   end do !iwan
   lmax2(isppol)=(lmax(isppol)+1)**2
 end do
 max_lmax2=maxval(lmax2(:))
!
!get ylmfac, factor used for rotations and hybrid orbitals
!
 ABI_ALLOCATE(ylmr_fac,(max_lmax2,mwan,nsppol))

 
 do isppol=1,nsppol
   if(spin.ne.0 .and. spin.ne.isppol) cycle
   call mlwfovlp_ylmfar(ylmr_fac(1:lmax2(isppol),1:nwan(isppol),isppol),&
&   lmax(isppol),lmax2(isppol),mband,nwan(isppol),proj_l(:,isppol),proj_m(:,isppol),&
&   proj_x(:,:,isppol),proj_z(:,:,isppol))
!  
!  Shift projection centers and atom centers to the primitive cell
!  This will be useful after, when we check if the Wannier function
!  lies on one specific atom
!  
   proj_site_unit(:,:,:)=0.d0
   do iwan=1,nwan(isppol)
     do ii=1,3
       proj_site_unit(ii,iwan,isppol)=ABS(proj_site(ii,iwan,isppol)-AINT(proj_site(ii,iwan,isppol)) )
     end do
   end do
   do iatom=1,natom
     do ii=1,3
       xred_unit(ii,iatom)=ABS(xred(ii,iatom)-AINT(xred(ii,iatom)) )
     end do
   end do
   call xred2xcart(natom,rprimd,xcart_unit,xred_unit)
   call xred2xcart(mwan,rprimd,proj_cart(:,:,isppol),proj_site_unit(:,:,isppol))
!  
!  Normalize the Wannier functions 
!  
!  Radial part
   mesh_size= nint((rmax - xmin ) / dx + 1)
   ABI_ALLOCATE( ff,(mesh_size))
   ABI_ALLOCATE(r,(mesh_size))
   ABI_ALLOCATE(rad_int,(mesh_size))
   ABI_ALLOCATE(aux,(mesh_size))
   do ir=1, mesh_size
     x=xmin+DBLE(ir-1)*dx
     r(ir)=x
   end do   !ir
   do iwan=1,nwan(isppol)
!    write(std_out,*)'iwan',iwan
!    radial functions shown in table 3.3 of wannier90 manual
     if(proj_radial(iwan,isppol)==1) ff(:) = 2.d0 * proj_zona(iwan,isppol)**(1.5d0) * exp(-proj_zona(iwan,isppol)*r(:))
     if(proj_radial(iwan,isppol)==2) ff(:) = 1.d0/(2.d0*sqrt(2.d0))*proj_zona(iwan,isppol)**(1.5d0) *&
&     (2.d0 - proj_zona(iwan,isppol)*r(:))*exp(-proj_zona(iwan,isppol)*r(:)/2.d0)
     if(proj_radial(iwan,isppol)==3) ff(:) = sqrt(4.d0/27.d0)*proj_zona(iwan,isppol)**(1.5d0)&
&     * (1.d0 - 2.d0*proj_zona(iwan,isppol)*r(:)/3.d0 + 2.d0*proj_zona(iwan,isppol)**2*r(:)**2/27.d0)&
&     * exp(-proj_zona(iwan,isppol) * r(:)/3.d0)

     if(proj_radial(iwan,isppol)/=4) then
       aux(:)=ff(:)**2*r(:)**2
       call simpson_int(mesh_size,dx,aux,rad_int)
       sum=0.d0
       do ir=1,mesh_size
         sum=sum+rad_int(ir)
       end do
       int_rad2=sum/real(mesh_size,dp)
!      
!      do ir=1,mesh_size
!      if(iwan==1) write(400,*)r(ir),aux(ir),rad_int(ir)
!      end do
     else
!      
!      ==4: gaussian function
!      f(x)=\exp(-1/4(x/aa)**2)
!      \int f(x)f(x) dx = \int \exp(-1/2(x/aa)**2) = aa*sqrt(2pi)
!      
       int_rad2=sqrt(2.d0*pi)*proj_zona(iwan,isppol)
     end if

!    
!    Now angular part
!    
     prod_real=0.d0
     do lm=1,lmax2(isppol)
       wan_lm_fac=ylmr_fac(lm,iwan,isppol)
!      write(std_out,*)'wan_lm_fac',wan_lm_fac
!      write(std_out,*)'int_rad2',int_rad2
       prod_real= prod_real + wan_lm_fac**2 * int_rad2  
     end do
     norm(iwan,isppol)=sqrt(prod_real)
   end do !iwan
   ABI_DEALLOCATE(ff)
   ABI_DEALLOCATE(r)
   ABI_DEALLOCATE(rad_int)
   ABI_DEALLOCATE(aux)
!  
!  Now that we found our guiding functions
!  We proceed with the internal product of
!  our guiding functions and the wave function
!  Amn=<G_m|\Psi_n> inside the sphere.
!  The term <G_m|\Psi_n> inside the sphere is:
!  = \sum_i <G_n | \phi_i - \tphi_i> <p_im|\Psi_m>
!  
!  
!  G_n \phi and \tphi can be decomposed in 
!  a radial function times an angular function. 
!  
!  
!  Big loop on iwan and iatom
!  
   do iwan=1,nwan(isppol)
     do iatom=1,natom
!      
!      check if center of wannier function coincides 
!      with the center of the atom
!      
       dist=((proj_cart(1,iwan,isppol)-xcart_unit(1,iatom))**2 + &
&       (proj_cart(2,iwan,isppol)-xcart_unit(2,iatom))**2 + &
&       (proj_cart(3,iwan,isppol)-xcart_unit(3,iatom))**2)**0.5 
!      
!      if the distance between the centers is major than 0.1 angstroms skip
!      
       if( dist > 0.188972613) cycle 
!      
       write(message, '(2a,i4,a,i4,2a)')ch10, '   Wannier function center',iwan,' is on top of atom',&
&       iatom,ch10,'      Calculating in-sphere contribution'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,  message,'COLL')
!      
!      Get useful quantities
!      
       itypat=typat(iatom)
       lmn_size=pawtab(itypat)%lmn_size
       basis_size=pawtab(itypat)%basis_size
       mesh_size=pawtab(itypat)%mesh_size
       ABI_ALLOCATE(int_rad,(basis_size))
       ABI_ALLOCATE(ff,(mesh_size))
       ABI_ALLOCATE(aux,(mesh_size))

!      
!      Integrate first the radial part
!      and save it into an array
!      
!      
!      radial functions shown in table 3.3 of wannier90 manual
!      
       if(proj_radial(iwan,isppol)==1) aux(1:mesh_size) = 2.d0 * proj_zona(iwan,isppol)**(1.5d0) *&
&       exp(-proj_zona(iwan,isppol)*pawrad(itypat)%rad(1:mesh_size))
       if(proj_radial(iwan,isppol)==2) aux(1:mesh_size) = 1.d0/(2.d0*sqrt(2.d0))*proj_zona(iwan,isppol)**(1.5d0) *&
&       (2.d0 - proj_zona(iwan,isppol)*pawrad(itypat)%rad(1:mesh_size)) &
&       * exp(-proj_zona(iwan,isppol)*pawrad(itypat)%rad(1:mesh_size)/2.d0)
       if(proj_radial(iwan,isppol)==3) aux(1:mesh_size) = sqrt(4.d0/27.d0)*proj_zona(iwan,isppol)**(1.5d0)&
&       * (1.d0 - 2.d0*proj_zona(iwan,isppol)*pawrad(itypat)%rad(1:mesh_size)/3.d0 &
&       + 2.d0*proj_zona(iwan,isppol)**2 *pawrad(itypat)%rad(1:mesh_size)**2/27.d0)&
&       * exp(-proj_zona(iwan,isppol) * pawrad(itypat)%rad(1:mesh_size)/3.d0)
!      
!      ==4: gaussian function
!      f(x)=\exp(-1/4(x/aa)**2)
!      
       if(proj_radial(iwan,isppol)==4) then
         aa=1.d0/proj_zona(iwan,isppol)
         aux(1:mesh_size)= exp(-0.25d0*(pawrad(itypat)%rad(1:mesh_size)*aa)**2)
       end if
!      
!      Normalize aux
       aux(:)=aux(:)/norm(iwan,isppol)
!      
       do ln=1,basis_size
         if(just_augmentation(iwan,isppol)) then
!          
!          just augmentation region contribution
!          In this case there is no need to use \tphi
!          ff= \int R_wan(r) (R_phi(ln;r)/r ) r^2 dr
!          
           ff(1:mesh_size)= aux(1:mesh_size) * pawtab(itypat)%phi(1:mesh_size,ln) &
&                                            * pawrad(itypat)%rad(1:mesh_size)
         else
!          Inside sphere contribution = \phi - \tphi
!          ff= \int R_wan(r) (R_phi(ln;r)/r - R_tphi(ln;r)/r) r^2 dr
           ff(1:mesh_size)= aux(1:mesh_size) * (pawtab(itypat)%phi(1:mesh_size,ln)-pawtab(itypat)%tphi(1:mesh_size,ln)) &
&                                            * pawrad(itypat)%rad(1:mesh_size)
         end if
!        
!        Integration with simpson routine
!        
         call simp_gen(int_rad(ln),ff,pawrad(itypat))
!        do ii=1,mesh_size
!        unit_ln=400+ln
!        if( iwan==1 ) write(unit_ln,*)pawrad(itypat)%rad(ii),ff(ii),int_rad(ln)
!        end do
       end do !ln
       ABI_DEALLOCATE(ff)
       ABI_DEALLOCATE(aux)
!      
!      Now integrate the angular part
!      Cycle on i indices
!      
!      prod_real=0.d0
       do ilmn=1, lmn_size
         ll=Psps%indlmn(1,ilmn,itypat)
         mm=Psps%indlmn(2,ilmn,itypat)
         nn=Psps%indlmn(3,ilmn,itypat)
         lm=Psps%indlmn(4,ilmn,itypat)
         ln=Psps%indlmn(5,ilmn,itypat)
!        write(std_out,*)'ll ',ll,' mm ',mm,'nn',nn,"lm",lm,"ln",ln
!        
!        Get wannier factor for that lm component
         if(lm <=lmax2(isppol)) then
           wan_lm_fac=ylmr_fac(lm,iwan,isppol)
!          Make delta product
!          Here we integrate the angular part
!          Since the integral of the product of two spherical harmonics 
!          is a delta function
           if( abs(wan_lm_fac) > 0.0d0) then
!            write(std_out,*) 'll',ll,'mm',mm,'lm',lm,'ln',ln,'factor',wan_lm_fac !lm index for wannier function
!            
!            Calculate Amn_paw, now that the radial and angular integrations are done
!            
             prod=cmplx(0.d0,0.d0)
             do ikpt=1,nkpt 
               jband=0
               do iband=1,nband(ikpt)
                 if(band_in(iband,isppol)) then
                   jband=jband+1

                   prod_real= cprj(iatom,index(iband,ikpt,isppol))%cp(1,ilmn) * int_rad(ln) * wan_lm_fac
                   prod_imag= cprj(iatom,index(iband,ikpt,isppol))%cp(2,ilmn) * int_rad(ln) * wan_lm_fac
                   prod=cmplx(prod_real,prod_imag)

                   A_paw(jband,iwan,ikpt,isppol)=A_paw(jband,iwan,ikpt,isppol)+prod
                 end if !band_in
               end do !iband
             end do !ikpt
!            
           end if !lm<=lmax2
         end if  ! abs(wan_lm_fac) > 0.0d0
       end do !ilmn=1, lmn_size
       ABI_DEALLOCATE(int_rad)
     end do !iatom
   end do !iwan
 end do !isppol
!
!Deallocate quantities
!
 ABI_DEALLOCATE(ylmr_fac)


!DEBUG
!write (std_out,*) ' mlwfovlp_projpaw : exit'
!stop
!ENDDEBUG

end subroutine mlwfovlp_projpaw
!!***
