!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_paw_atom
!! NAME
!!  m_paw_atom
!!
!! FUNCTION
!!  atompaw related operations
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2019 ABINIT group (T. Rangel, MT, JWZ, GJ)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

module m_paw_atom

 USE_DEFS
 USE_MSG_HANDLING
 USE_MEMORY_PROFILING

 use m_paw_numeric, only : paw_jbessel, paw_solvbes, paw_spline, paw_splint
 use m_pawrad,      only : pawrad_type, simp_gen, poisson, pawrad_deducer0, bound_deriv, pawrad_ifromr
 use m_pawtab,      only : pawtab_type

 implicit none

 private

 public:: atompaw_shpfun
 public:: atompaw_shapebes
 public:: atompaw_vhnzc
 public:: atompaw_dij0
 public:: atompaw_kij
!!***

CONTAINS !===========================================================
!!***

!!****f* m_paw_atom/atompaw_shpfun
!! NAME
!! atompaw_shpfun
!!
!! FUNCTION
!! Compute shape function used in the definition
!! of compensation density (PAW)
!!
!! INPUTS
!!  ll= l quantum number
!!  mesh <type(pawrad_type)>=data containing radial grid information
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!
!! OUTPUT
!!  norm= factor for shape function normalization
!!
!! SIDE effects
!!  shapefunc(:)=shape function g(r)
!!    In case of numerical shape function (shape_type=-1), shapefunc
!!    array contains the shape function read in psp file at input.
!!
!! NOTES
!!  Types of shape functions:
!!   type -1: numerical shape function, given in psp file
!!   type  1: g(r)=k(r).r^l; k(r)=exp(-(r/sigma)^lambda)
!!   type  2: g(r)=k(r).r^l; k(r)=[sin(Pi.r/rshp)/(Pi.r/rshp)]^2
!!   type  3: g(r)=alpha1.jl(q1.r)+alpha2.jl(q2.r)
!!
!! PARENTS
!!      m_paw_atom,m_pawpsp,pawinit
!!
!! CHILDREN
!!      atompaw_shpfun,atompaw_vhnzc,bound_deriv,paw_spline,paw_splint,simp_gen
!!
!! SOURCE

subroutine atompaw_shpfun(ll,mesh,norm,pawtab,shapefunc)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll
 real(dp),intent(out) :: norm
 type(pawrad_type),intent(in) :: mesh
 type(pawtab_type),intent(in) :: pawtab
!arrays
 real(dp),intent(inout) :: shapefunc(:)

!Local variables ------------------------------
!scalars
 integer :: ir,ishp,mesh_size
 real(dp) :: arg,besp,bespp,jbes1,jbes2
!arrays
 real(dp) :: alpha(2),qq(2)
 real(dp),allocatable :: r2k(:)
!no_abirules

!***************************************************************************

 mesh_size=size(shapefunc)
 if (mesh_size>mesh%mesh_size) then
   MSG_BUG('wrong size!')
 end if

!Index for shape function cut-off radius
 ishp=pawrad_ifromr(mesh,pawtab%rshp)-1

!Computation of non-normalized shape function
 if (pawtab%shape_type==-1) then
   shapefunc(1:ishp)=pawtab%shapefunc(1:ishp,1+ll)
 else if (pawtab%shape_type==1) then
   if (ll==0) then
     shapefunc(1)=one
     do ir=2,ishp
       shapefunc(ir)=exp(-(mesh%rad(ir)/pawtab%shape_sigma)**pawtab%shape_lambda)
     end do
   else
     shapefunc(1)=zero
     do ir=2,ishp
       shapefunc(ir)=exp(-(mesh%rad(ir)/pawtab%shape_sigma)**pawtab%shape_lambda)*mesh%rad(ir)**ll
     end do
   end if
 else if (pawtab%shape_type==2) then
   if (ll==0) then
     shapefunc(1)=one
     do ir=2,ishp
       arg=pi*mesh%rad(ir)/pawtab%rshp
       shapefunc(ir)=(sin(arg)/arg)**2
     end do
   else
     shapefunc(1)=zero
     do ir=2,ishp
       arg=pi*mesh%rad(ir)/pawtab%rshp
       shapefunc(ir)=(sin(arg)/(arg))**2 *mesh%rad(ir)**ll
     end do
   end if
 else if (pawtab%shape_type==3) then
   alpha(1:2)=pawtab%shape_alpha(1:2,1+ll)
   qq(1:2)=pawtab%shape_q(1:2,1+ll)
   do ir=1,ishp
     call paw_jbessel(jbes1,besp,bespp,ll,0,qq(1)*mesh%rad(ir))
     call paw_jbessel(jbes2,besp,bespp,ll,0,qq(2)*mesh%rad(ir))
     shapefunc(ir)=alpha(1)*jbes1+alpha(2)*jbes2
   end do
 end if

 if (ishp<mesh_size) shapefunc(ishp+1:mesh_size)=zero

!Shape function normalization
 if (pawtab%shape_type==-1.or.pawtab%shape_type==1.or.pawtab%shape_type==2) then
   LIBPAW_ALLOCATE(r2k,(mesh_size))
   r2k=zero
   r2k(2:ishp)=shapefunc(2:ishp)*mesh%rad(2:ishp)**(2+ll)
   if (mesh%mesh_type==5) then
     call simp_gen(norm,r2k,mesh);norm=one/norm
   else
     call simp_gen(norm,r2k,mesh,r_for_intg=pawtab%rshp);norm=one/norm
   end if
   shapefunc(1:ishp)=shapefunc(1:ishp)*norm
   if (pawtab%shape_type==-1) norm=one
   LIBPAW_DEALLOCATE(r2k)
 else if (pawtab%shape_type==3) then
   norm=one
 end if

end subroutine atompaw_shpfun
!!***

!----------------------------------------------------------------------

!!****f* m_paw_atom/atompaw_atompaw_shapebes
!! NAME
!! atompaw_shapebes
!!
!! FUNCTION
!!    Find al and ql parameters for a "Bessel" shape function:
!!    Shape(r)=al1.jl(ql1.r)+al2.jl(ql2.r)
!!      such as Shape(r) and 2 derivatives are zero at r=rc
!!              Intg_0_rc[Shape(r).r^(l+2).dr]=1
!!
!! INPUTS
!!  ll= l quantum number
!!  rc= cut-off radius
!!
!! OUTPUT
!!  al(2)= al coefficients
!!  ql(2)= ql factors
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!      atompaw_shpfun,atompaw_vhnzc,bound_deriv,paw_spline,paw_splint,simp_gen
!!
!! SOURCE

 subroutine atompaw_shapebes(al,ql,ll,rc)

!Arguments ------------------------------------
!scalars
 integer :: ll
 real(dp) :: rc
!arrays
 real(dp) :: al(2),ql(2)

!Local variables-------------------------------
!scalars
 integer :: ii
 real(dp) :: alpha,beta,det,jbes,jbesp,jbespp,qr
!arrays
 real(dp) :: amat(2,2),bb(2)

! *************************************************************************

 alpha=1._dp;beta=0._dp
 call paw_solvbes(ql,alpha,beta,ll,2)
 ql(1:2)=ql(1:2)/rc

 do ii=1,2
   qr=ql(ii)*rc
   call paw_jbessel(jbes,jbesp,jbespp,ll,1,qr)
   amat(1,ii)=jbesp*ql(ii)
   call paw_jbessel(jbes,jbesp,jbespp,ll+1,0,qr)
   amat(2,ii)=jbes*rc**(ll+2)/ql(ii)  !  Intg_0_rc[jl(qr).r^(l+2).dr]
 end do

 bb(1)=zero;bb(2)=one

 det=amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1)
 al(1)=(amat(2,2)*bb(1)-amat(1,2)*bb(2))/det
 al(2)=(amat(1,1)*bb(2)-amat(2,1)*bb(1))/det

end subroutine atompaw_shapebes
!!***

!----------------------------------------------------------------------

!!****f* m_paw_atom/atompaw_vhnzc
!! NAME
!! atompaw_vhnzc
!!
!! FUNCTION
!! PAW: compute Hartree potential for n_{Zc}
!!
!! INPUTS
!!  ncore(:)=atomic core density
!!  radmesh_core <type(pawrad_type)>=radial mesh (and related data) for the core densities
!!  znucl= valence and total charge of the atomic species
!!
!! OUTPUT
!!  vhnzc(:)=Hartree potential due to Z_nc
!!
!! PARENTS
!!      m_paw_atom,m_pawpsp
!!
!! CHILDREN
!!      atompaw_shpfun,atompaw_vhnzc,bound_deriv,paw_spline,paw_splint,simp_gen
!!
!! SOURCE

 subroutine atompaw_vhnzc(ncore,radmesh_core,vhnzc,znucl)

!Arguments ---------------------------------------------
!scalars
 real(dp),intent(in) :: znucl
 type(pawrad_type),intent(in) :: radmesh_core
!arrays
 real(dp),intent(in) :: ncore(:)
 real(dp), intent(out) :: vhnzc(:)

!Local variables ---------------------------------------
 integer :: mesh_size
 real(dp),allocatable :: nwk(:)

! *********************************************************************

 mesh_size=size(ncore)
 if (mesh_size/=size(vhnzc).or.mesh_size>radmesh_core%mesh_size) then
   MSG_BUG('wrong sizes!')
 end if

 LIBPAW_ALLOCATE(nwk,(mesh_size))

 nwk(:)=ncore(:)*four_pi*radmesh_core%rad(:)**2
 call poisson(nwk,0,radmesh_core,vhnzc)
 vhnzc(2:mesh_size)=(vhnzc(2:mesh_size)-znucl)/radmesh_core%rad(2:mesh_size)
 call pawrad_deducer0(vhnzc,mesh_size,radmesh_core)

 LIBPAW_DEALLOCATE(nwk)

 end subroutine atompaw_vhnzc
!!***

!----------------------------------------------------------------------

!!****f* m_paw_atom/atompaw_dij0
!! NAME
!! atompaw_dij0
!!
!! FUNCTION
!!  PAW: Compute "frozen" values of pseudopotential strengths Dij = Dij0
!!
!! INPUTS
!!  indlmn(6,lmnmax)= array giving l,m,n,lm,ln,s for i=lmn
!!  kij(pawtab%lmn2_size)= kinetic part of Dij
!!  lmnmax=max number of (l,m,n) components over all type of psps
!!  ncore(:)=atomic core density
!!  opt_init=flag defining the storage of PAW atomic data
!!           0: PAW atomic data have not been initialized (in pawtab)
!!           1: PAW atomic data have been initialized (in pawtab)
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  radmesh <type(pawrad_type)>=paw radial mesh (and related data)
!!  radmesh_core <type(pawrad_type)>=radial mesh (and related data) for the core densities
!!  radmesh_vloc <type(pawrad_type)>=radial mesh (and related data) for the local potential (VH(tnZc))
!!  vhtnzc(:)= local potential VH(tnZc)
!!  znucl= valence and total charge of the atomic species
!!
!! OUTPUT
!!  pawtab%dij0(pawtab%lmn2_size)= Frozen part of the Dij term
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!      atompaw_shpfun,atompaw_vhnzc,bound_deriv,paw_spline,paw_splint,simp_gen
!!
!! SOURCE


 subroutine atompaw_dij0(indlmn,kij,lmnmax,ncore,opt_init,pawtab,&
&                        radmesh,radmesh_core,radmesh_vloc,vhtnzc,znucl)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lmnmax,opt_init
 real(dp),intent(in) :: znucl
 type(pawrad_type),intent(in) :: radmesh,radmesh_core,radmesh_vloc
 type(pawtab_type),intent(inout) :: pawtab
!arrays
 integer,intent(in) :: indlmn(6,lmnmax)
 real(dp),intent(in) :: kij(pawtab%lmn2_size)
 real(dp),intent(in) :: ncore(:),vhtnzc(:)
!real(dp),optional,intent(in) :: vminushalf(:)

!Local variables ---------------------------------------
 integer :: il,ilm,iln,ilmn,j0lmn,jl,jlm,jln,jlmn,klmn,lmn2_size,meshsz,meshsz_core
 integer :: meshsz_vhtnzc,meshsz_vmh
 real(dp) :: intg,intvh,yp1,ypn
 real(dp),allocatable :: ff(:),r2k(:),shpf(:),vhnzc(:),vhtnzc_sph(:),work1(:),work2(:)

! *********************************************************************

 lmn2_size=pawtab%lmn2_size
 meshsz_vhtnzc=size(vhtnzc)
 meshsz=min(radmesh%mesh_size,radmesh_core%mesh_size,radmesh_vloc%mesh_size,meshsz_vhtnzc)
 LIBPAW_ALLOCATE(ff,(meshsz))

!Retrieve VH(tnZc) on the correct radial mesh
 LIBPAW_ALLOCATE(vhtnzc_sph,(meshsz))
 if ((radmesh%mesh_type/=radmesh_vloc%mesh_type).or.&
&    (radmesh%rstep    /=radmesh_vloc%rstep)    .or.&
&    (radmesh%lstep    /=radmesh_vloc%lstep)) then
   call bound_deriv(vhtnzc,radmesh_vloc,meshsz_vhtnzc,yp1,ypn)
   LIBPAW_ALLOCATE(work1,(meshsz_vhtnzc))
   LIBPAW_ALLOCATE(work2,(meshsz_vhtnzc))
   call paw_spline(radmesh_vloc%rad,vhtnzc,meshsz_vhtnzc,yp1,ypn,work1)
   call paw_splint(meshsz_vhtnzc,radmesh_vloc%rad,vhtnzc,work1,meshsz,radmesh%rad(1:meshsz),vhtnzc_sph)
   LIBPAW_DEALLOCATE(work1)
   LIBPAW_DEALLOCATE(work2)
 else
   vhtnzc_sph(1:meshsz)=vhtnzc(1:meshsz)
 end if

!Kinetic part of Dij0
!====================
 pawtab%dij0(1:lmn2_size)=kij(1:lmn2_size)

!Computation of <phi_i|vh(nZc)|phi_j> on the PAW sphere
!======================================================
 meshsz_core=size(ncore)
 LIBPAW_ALLOCATE(vhnzc,(meshsz_core))
 call atompaw_vhnzc(ncore,radmesh_core,vhnzc,znucl)
 do jlmn=1,pawtab%lmn_size
   j0lmn=jlmn*(jlmn-1)/2
   jlm=indlmn(4,jlmn);jln=indlmn(5,jlmn)
   do ilmn=1,jlmn
     klmn=j0lmn+ilmn
     ilm=indlmn(4,ilmn);iln=indlmn(5,ilmn)
     if (jlm==ilm) then
       ff(1:meshsz)=pawtab%phi(1:meshsz,iln)*pawtab%phi(1:meshsz,jln)*vhnzc(1:meshsz)
       call simp_gen(intg,ff,radmesh)
       pawtab%dij0(klmn)=pawtab%dij0(klmn)+intg
     end if
   end do
 end do
 LIBPAW_DEALLOCATE(vhnzc)

!Computation of -<tphi_i|vh(tnZc)|tphi_j> on the PAW sphere
!==========================================================
 do jlmn=1,pawtab%lmn_size
   j0lmn=jlmn*(jlmn-1)/2
   jlm=indlmn(4,jlmn);jln=indlmn(5,jlmn)
   do ilmn=1,jlmn
     klmn=j0lmn+ilmn
     ilm=indlmn(4,ilmn);iln=indlmn(5,ilmn)
     if (jlm==ilm) then
       ff(1:meshsz)=pawtab%tphi(1:meshsz,iln)*pawtab%tphi(1:meshsz,jln)*vhtnzc_sph(1:meshsz)
       call simp_gen(intg,ff,radmesh)
       pawtab%dij0(klmn)=pawtab%dij0(klmn)-intg
     end if
   end do
 end do

!Computation of <phi_i|vminushalf|phi_j>  (if any)
!=================================================
 if(pawtab%has_vminushalf==1) then
   if(size(pawtab%vminushalf)>=1) then
     meshsz_vmh=min(meshsz,size(pawtab%vminushalf))
     do jlmn=1,pawtab%lmn_size
       j0lmn=jlmn*(jlmn-1)/2
       jlm=indlmn(4,jlmn);jln=indlmn(5,jlmn)
       do ilmn=1,jlmn
         klmn=j0lmn+ilmn
         ilm=indlmn(4,ilmn);iln=indlmn(5,ilmn)
         if (jlm==ilm) then
           ff(1:meshsz_vmh)=pawtab%phi(1:meshsz_vmh,iln)*pawtab%phi(1:meshsz_vmh,jln)*pawtab%vminushalf(1:meshsz_vmh)
           call simp_gen(intg,ff(1:meshsz_vmh),radmesh)
           pawtab%dij0(klmn)=pawtab%dij0(klmn)+intg
         end if
       end do
     end do
   end if
 end if

!Computation of -int[vh(tnzc)*Qijhat(r)dr]
!=========================================
 if (opt_init==0) then
   LIBPAW_ALLOCATE(shpf,(radmesh%mesh_size))
   call atompaw_shpfun(0,radmesh,intg,pawtab,shpf)
   if (pawtab%shape_type==3) then
     LIBPAW_ALLOCATE(r2k,(radmesh%int_meshsz))
     r2k=zero
     r2k(2:radmesh%int_meshsz)=shpf(2:radmesh%int_meshsz)*radmesh%rad(2:radmesh%int_meshsz)**2
     if(radmesh%mesh_type==5) then
       call simp_gen(intg,r2k,radmesh)
     else
       call simp_gen(intg,r2k,radmesh,r_for_intg=pawtab%rshp)
     end if
     shpf(1:meshsz)=shpf(1:meshsz)/intg
     LIBPAW_DEALLOCATE(r2k)
   end if
   ff(1:meshsz)=vhtnzc_sph(1:meshsz)*shpf(1:meshsz)*radmesh%rad(1:meshsz)**2
   LIBPAW_DEALLOCATE(shpf)
   call simp_gen(intvh,ff,radmesh)
   do jlmn=1,pawtab%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     jl=indlmn(1,jlmn);jln=indlmn(5,jlmn);jlm=indlmn(4,jlmn)
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       il=indlmn(1,ilmn);iln=indlmn(5,ilmn);ilm=indlmn(4,ilmn)
       if (ilm==jlm) then
         ff(1:meshsz)=(pawtab%phi (1:meshsz,iln)*pawtab%phi (1:meshsz,jln)&
&                     -pawtab%tphi(1:meshsz,iln)*pawtab%tphi(1:meshsz,jln))
         call simp_gen(intg,ff,radmesh)
         pawtab%dij0(klmn)=pawtab%dij0(klmn)-intvh*intg
       end if
     end do
   end do
 else
   ff(1:meshsz)=vhtnzc_sph(1:meshsz)*pawtab%shapefunc(1:meshsz,1)*radmesh%rad(1:meshsz)**2
   call simp_gen(intvh,ff,radmesh)
   do jlmn=1,pawtab%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     jl=indlmn(1,jlmn);jln=indlmn(5,jlmn);jlm=indlmn(4,jlmn)
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       il=indlmn(1,ilmn);iln=indlmn(5,ilmn);ilm=indlmn(4,ilmn)
       if (ilm==jlm) then
         intg=pawtab%qijl(1,klmn)*sqrt(four_pi)
         pawtab%dij0(klmn)=pawtab%dij0(klmn)-intvh*intg
       end if
     end do
   end do
 end if

 LIBPAW_DEALLOCATE(ff)
 LIBPAW_DEALLOCATE(vhtnzc_sph)

 end subroutine atompaw_dij0
!!***

!----------------------------------------------------------------------

!!****f* m_paw_atom/atompaw_kij
!! NAME
!! atompaw_kij
!!
!! FUNCTION
!! PAW: deduce kinetic part of psp strength (Dij) from the knowledge of frozen Dij (Dij0)
!!
!! INPUTS
!!  indlmn(6,lmnmax)= array giving l,m,n,lm,ln,s for i=lmn
!!  lmnmax=max number of (l,m,n) components over all type of psps
!!  ncore(:)=atomic core density
!!  opt_init=flag defining the storage of PAW atomic data
!!           0: PAW atomic data have not been initialized (in pawtab)
!!           1: PAW atomic data have been initialized (in pawtab)
!!  opt_vhnzc=flag defining the inclusion of VH(nZc) in computation
!!            0: VH(nZc) is not taken into account
!!            1: VH(nZc) is taken into account
!!  pawtab <type(pawtab_type)>=paw tabulated starting data
!!  radmesh <type(pawrad_type)>=paw radial mesh (and related data)
!!  radmesh_core <type(pawrad_type)>=radial mesh (and related data) for the core densities
!!  radmesh_vloc <type(pawrad_type)>=radial mesh (and related data) for the local potential (VH(tnZc))
!!  vhtnzc(:)= local potential VH(tnZc)
!!  znucl= valence and total charge of the atomic species
!!
!! OUTPUT
!!  kij(pawtab%lmn2_size)= kinetic part of Dij
!!
!! PARENTS
!!      m_pawpsp
!!
!! CHILDREN
!!      atompaw_shpfun,atompaw_vhnzc,bound_deriv,paw_spline,paw_splint,simp_gen
!!
!! SOURCE

 subroutine atompaw_kij(indlmn,kij,lmnmax,ncore,opt_init,opt_vhnzc,pawtab, &
&                       radmesh,radmesh_core,radmesh_vloc,vhtnzc,znucl)

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lmnmax,opt_init,opt_vhnzc
 real(dp),intent(in) :: znucl
 type(pawrad_type),intent(in) :: radmesh,radmesh_core,radmesh_vloc
 type(pawtab_type),intent(in) :: pawtab
!arrays
 integer,intent(in) :: indlmn(6,lmnmax)
 real(dp),intent(out) :: kij(pawtab%lmn2_size)
 real(dp),intent(in) :: ncore(:)
 real(dp),intent(in) :: vhtnzc(:)

!Local variables ---------------------------------------
 integer :: il,ilm,iln,ilmn,j0lmn,jl,jlm,jln,jlmn,klmn,lmn2_size
 integer :: meshsz,meshsz_core,meshsz_vhtnzc,meshsz_vmh
 real(dp) :: intg,intvh,yp1,ypn
 real(dp),allocatable :: ff(:),r2k(:),shpf(:),vhnzc(:),vhtnzc_sph(:),work1(:),work2(:)

! *********************************************************************

 lmn2_size=pawtab%lmn2_size
 meshsz_vhtnzc=size(vhtnzc)
 meshsz=min(radmesh%mesh_size,radmesh_core%mesh_size,radmesh_vloc%mesh_size,meshsz_vhtnzc)
 LIBPAW_ALLOCATE(ff,(meshsz))

!Retrieve VH(tnZc) on the correct radial mesh
 LIBPAW_ALLOCATE(vhtnzc_sph,(meshsz))
 if ((radmesh%mesh_type/=radmesh_vloc%mesh_type).or.&
&    (radmesh%rstep    /=radmesh_vloc%rstep)    .or.&
&    (radmesh%lstep    /=radmesh_vloc%lstep)) then
   call bound_deriv(vhtnzc,radmesh_vloc,meshsz_vhtnzc,yp1,ypn)
   LIBPAW_ALLOCATE(work1,(meshsz_vhtnzc))
   LIBPAW_ALLOCATE(work2,(meshsz_vhtnzc))
   call paw_spline(radmesh_vloc%rad,vhtnzc,meshsz_vhtnzc,yp1,ypn,work1)
   call paw_splint(meshsz_vhtnzc,radmesh_vloc%rad,vhtnzc,work1,meshsz,radmesh%rad(1:meshsz),vhtnzc_sph)
   LIBPAW_DEALLOCATE(work1)
   LIBPAW_DEALLOCATE(work2)
 else
   vhtnzc_sph(1:meshsz)=vhtnzc(1:meshsz)
 end if

!Initialize Kij with Dij0
!=========================
 kij(1:lmn2_size)=pawtab%dij0(1:lmn2_size)

!Substraction of -<phi_i|vh(nZc)|phi_j> on the PAW sphere
!========================================================
 if (opt_vhnzc/=0) then
   meshsz_core=size(ncore)
   LIBPAW_ALLOCATE(vhnzc,(meshsz_core))
   call atompaw_vhnzc(ncore,radmesh_core,vhnzc,znucl)
   do jlmn=1,pawtab%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     jlm=indlmn(4,jlmn);jln=indlmn(5,jlmn)
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       ilm=indlmn(4,ilmn);iln=indlmn(5,ilmn)
       if (jlm==ilm) then
         ff(1:meshsz)=pawtab%phi(1:meshsz,iln)*pawtab%phi(1:meshsz,jln)*vhnzc(1:meshsz)
         call simp_gen(intg,ff,radmesh)
         kij(klmn)=kij(klmn)-intg
       end if
     end do
   end do
   LIBPAW_DEALLOCATE(vhnzc)
 end if

!Substraction of <tphi_i|vh(tnZc)|tphi_j> on the PAW sphere
!==========================================================
 do jlmn=1,pawtab%lmn_size
   j0lmn=jlmn*(jlmn-1)/2
   jlm=indlmn(4,jlmn);jln=indlmn(5,jlmn)
   do ilmn=1,jlmn
     klmn=j0lmn+ilmn
     ilm=indlmn(4,ilmn);iln=indlmn(5,ilmn)
     if (jlm==ilm) then
       ff(1:meshsz)=pawtab%tphi(1:meshsz,iln)*pawtab%tphi(1:meshsz,jln)*vhtnzc_sph(1:meshsz)
       call simp_gen(intg,ff,radmesh)
       kij(klmn)=kij(klmn)+intg
     end if
   end do
 end do

!Computation of <phi_i|vminushalf|phi_j>  (if any)
!=================================================
 if(pawtab%has_vminushalf==1) then
   if(size(pawtab%vminushalf)>=1) then
     meshsz_vmh=min(meshsz,size(pawtab%vminushalf))
     do jlmn=1,pawtab%lmn_size
       j0lmn=jlmn*(jlmn-1)/2
       jlm=indlmn(4,jlmn);jln=indlmn(5,jlmn)
       do ilmn=1,jlmn
         klmn=j0lmn+ilmn
         ilm=indlmn(4,ilmn);iln=indlmn(5,ilmn)
         if (jlm==ilm) then
           ff(1:meshsz_vmh)=pawtab%phi(1:meshsz_vmh,iln)*pawtab%phi(1:meshsz_vmh,jln)*pawtab%vminushalf(1:meshsz_vmh)
           call simp_gen(intg,ff(1:meshsz_vmh),radmesh)
           kij(klmn)=kij(klmn)-intg
         end if
       end do
     end do
   end if
 end if

!Computation of int[vh(tnzc)*Qijhat(r)dr]
!==========================================
 if (opt_init==0) then
   LIBPAW_ALLOCATE(shpf,(radmesh%mesh_size))
   call atompaw_shpfun(0,radmesh,intg,pawtab,shpf)
   if (pawtab%shape_type==3) then
     LIBPAW_ALLOCATE(r2k,(radmesh%int_meshsz))
     r2k=zero
     r2k(2:radmesh%int_meshsz)=shpf(2:radmesh%int_meshsz)*radmesh%rad(2:radmesh%int_meshsz)**2
     if(radmesh%mesh_type==5) then
       call simp_gen(intg,r2k,radmesh)
     else
       call simp_gen(intg,r2k,radmesh,r_for_intg=pawtab%rshp)
     end if
     shpf(1:meshsz)=shpf(1:meshsz)/intg
     LIBPAW_DEALLOCATE(r2k)
   end if
   ff(1:meshsz)=vhtnzc_sph(1:meshsz)*shpf(1:meshsz)*radmesh%rad(1:meshsz)**2
   LIBPAW_DEALLOCATE(shpf)
   call simp_gen(intvh,ff,radmesh)
   do jlmn=1,pawtab%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     jl=indlmn(1,jlmn);jln=indlmn(5,jlmn);jlm=indlmn(4,jlmn)
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       il=indlmn(1,ilmn);iln=indlmn(5,ilmn);ilm=indlmn(4,ilmn)
       if (ilm==jlm) then
         ff(1:meshsz)=(pawtab%phi (1:meshsz,iln)*pawtab%phi (1:meshsz,jln)&
&                     -pawtab%tphi(1:meshsz,iln)*pawtab%tphi(1:meshsz,jln))
         call simp_gen(intg,ff,radmesh)
         kij(klmn)=kij(klmn)+intvh*intg
       end if
     end do
   end do
 else
   ff(1:meshsz)=vhtnzc_sph(1:meshsz)*pawtab%shapefunc(1:meshsz,1)*radmesh%rad(1:meshsz)**2
   call simp_gen(intvh,ff,radmesh)
   do jlmn=1,pawtab%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     jl=indlmn(1,jlmn);jln=indlmn(5,jlmn);jlm=indlmn(4,jlmn)
     do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       il=indlmn(1,ilmn);iln=indlmn(5,ilmn);ilm=indlmn(4,ilmn)
       if (ilm==jlm) then
         intg=pawtab%qijl(1,klmn)*sqrt(four_pi)
         kij(klmn)=kij(klmn)+intvh*intg
       end if
     end do
   end do
 end if

 LIBPAW_DEALLOCATE(ff)
 LIBPAW_DEALLOCATE(vhtnzc_sph)

 end subroutine atompaw_kij
!!***

end module m_paw_atom
!!***
