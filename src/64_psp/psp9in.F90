!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp9in
!! NAME
!! psp9in
!!
!! FUNCTION
!! Initialize pspcod=9 (pseudopotentials from the PSML XML format):
!! continue to read the corresponding file, then compute the
!! local and non-local potentials.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (JJ, MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  filpsp=filename of the PSML pseudopotential
!!  lloc=angular momentum choice of local pseudopotential
!!  lmax=value of lmax mentioned at the second line of the psp file
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  lnmax=max. number of (l,n) components over all type of psps
!!  mmax=maximum number of points in real space grid in the psp file
!!   angular momentum of nonlocal pseudopotential
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 2*maximum angular momentum for nonlocal pseudopotentials - 1
!!  mqgrid=dimension of q (or G) grid for arrays.
!!  mqgrid_vl=dimension of q (or G) grid for valence charge (array qgrid_vl)
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  qgrid(mqgrid)=values of q (or |G|) on grid from 0 to qmax
!!  qgrid_vl(psps%mqgrid_vl)=values of q on grid from 0 to qmax (bohr^-1) for valence charge
!!  pspso=spin-orbit characteristics, govern the content of ffspl and ekb
!!   if =0 : this input requires NO spin-orbit characteristics of the psp
!!   if =2 : this input requires HGH or psp8 characteristics of the psp
!!   if =3 : this input requires HFN characteristics of the psp
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  zion=nominal valence of atom as specified in psp file
!!  znucl=nuclear number of atom as specified in psp file
!!
!! OUTPUT
!!  ekb(lnmax)=Kleinman-Bylander energy,
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for each (l,n)
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r} dr]$ (hartree)
!!  ffspl(mqgrid,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  nproj(mpssoang)=number of projection functions for each angular momentum
!!  qchrg is not used, and could be suppressed later
!!  vlspl(mqgrid,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  xcccrc=XC core correction cutoff radius (bohr)
!!  xccc1d(n1xccc,6)=1D core charge function and five derivatives, from psp file
!!  nctab<nctab_t>=NC tables
!!    %has_tvale=True if the pseudo contains the pseudo valence charge
!!    %tvalespl(mqgrid_vl,2)=the pseudo valence density and 2nd derivative in reciprocal space on a regular grid 
!!
!! PARENTS
!!      pspatm
!!
!! CHILDREN
!!      nctab_eval_tvalespl,pawrad_free,pawrad_init,ps_destroy
!!      ps_get_projector_indexes,psml_reader,psp8lo,psp8nl,psp9cc,spline,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine psp9in(filpsp,ekb,epsatm,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&
&                  mmax,mpsang,mpssoang,mqgrid,mqgrid_vl,nproj,n1xccc,pspso,qchrg,qgrid,qgrid_vl,&
&                  useylm,vlspl,xcccrc,xccc1d,zion,znucl,nctab,maxrad)

 use defs_basis
 use m_splines
 use m_errors
 use m_profiling_abi

 use defs_datatypes,  only : nctab_t
 use m_pawrad,        only : pawrad_type, pawrad_init, pawrad_free
 use m_psps,          only : nctab_eval_tvalespl

#if defined HAVE_PSML
 use m_psml
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'psp9in'
 use interfaces_14_hidewrite
 use interfaces_64_psp, except_this_one => psp9in
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lloc,lmax,lmnmax,lnmax,mpsang,mpssoang,mqgrid,mqgrid_vl
 integer,intent(in) :: pspso,n1xccc,useylm
 integer,intent(out) :: mmax
 real(dp),intent(in) :: zion,znucl
 real(dp),intent(out) :: epsatm,qchrg,xcccrc,maxrad
 type(nctab_t),intent(inout) :: nctab
 character(len=fnlen),intent(in) :: filpsp
!arrays
 integer,intent(out) :: indlmn(6,lmnmax),nproj(mpssoang)
 real(dp),intent(in) :: qgrid(mqgrid),qgrid_vl(mqgrid_vl)
 real(dp),intent(out) :: ekb(lnmax),ffspl(mqgrid,2,lnmax),vlspl(mqgrid,2)
 real(dp),intent(inout) :: xccc1d(n1xccc,6) !vz_i

!Local variables-------------------------------
!no_abirules
!scalars
#if defined HAVE_PSML
 character(len=500) :: message
 character(len=30)  :: creator
 character(len=7), parameter  :: oncvpsp_name = "ONCVPSP"
 integer :: iln,pspindex,ipsang,irad,kk,ll
 integer :: mm,nn,nso,ii,ir,il,nshells
 integer :: iproj,irelt,nders
 integer :: np_dn, np_lj, np_nr, np_so, np_sr, np_up, val_l, val_n
 logical :: has_nlcc,has_spin,has_tvale,oncvpsp
 real(dp) :: amesh,damesh,fchrg,rchrg,yp1,ypn,zval
 real(dp) :: rmax,rmatch,z,chgvps,val_occ
 type(pawrad_type) :: mesh
!arrays
 integer, allocatable :: idx_so(:),idx_sr(:)
 real(dp),allocatable :: rad(:),vloc(:),vpspll(:,:),work_spl(:)
 type(ps_t) :: psxml
#endif

! ***************************************************************************

#if defined HAVE_PSML

 call ps_destroy(psxml)
 call psml_reader(filpsp,psxml,debug=.true.)

!Identify the atomic code that generated the pseudopotential
 call ps_Provenance_Get(psxml, 1, creator=creator)
!Check whether the pseudopotential has been created with ONCVPSP,
!Don Hamann's code
 oncvpsp = (trim(creator(1:7)) .eq. trim(oncvpsp_name))
!DEBUG
!write(std_out,*)' psp9in : creator : ', creator
!write(std_out,*)' psp9in : oncvpsp : ', oncvpsp
!ENDDEBUG

! SIESTA's ATOM uses spherical harmonics, while ONCVPSP uses Legendre
! polynomials, which means we have to check the consistency of input variables
! wrt the pseudos
!
! Note: commented because NC pseudos do not have non-diagonal terms
!
! if ( oncvpsp ) then
!   if ( useylm /= 0 ) then
!     write(message,'(3a)') "ONCVPSP pseudos use Legendre polynomials but we use spherical harmonics", &
!&      ch10, "ACTION: set useylm to 0 in your input file"
!     MSG_ERROR(message)
!   endif
! else
!   if ( useylm == 0 ) then
!     write(message,'(3a)') "ATOM pseudos use spherical harmonics but we use Legendre polynomials", &
!&      ch10, "ACTION: set useylm to 1 in your input file"
!     MSG_ERROR(message)
!   endif
! endif

! The atomic number is a real number instead of a simple integer
! z (in Abinit), atomic-number in the header of the PSML file.
! z      = ps_AtomicNumber(psxml)
!
! The difference between the number of protons in the nucleus and the
! sum of the populations of the core shells is the effective atomic number 
! of the pseudo-atom, Zval (in Abinit), z-pseudo in the header of the
! PSML file.
! zval   = ps_Zpseudo(psxml)

 call ps_PseudoAtomSpec_Get(psxml, &
&  atomic_number=z, z_pseudo=zval, &
&  spin_dft=has_spin, core_corrections=has_nlcc)

!---

!Feb 2015: shifted to Hamann grid for convenience - libpsml interpolates anyway
!
! The following lines are taken from the oncvpsp.f90 subroutine of the oncvpsp
! code implemented by D. Hamann
! The atomic number of the element is read from the header of the XML file
! Logarithmic grid defined by Hamann in oncvpsp code
! z    = psxml%header%z
! amesh = 1.012d0
! al    = dlog(amesh)
! rr1   = .0005d0/z
! mmax  = dlog(45.0d0 /rr1)/al
!
! ABI_ALLOCATE( rad,(mmax) )
!
! do ir = 1, mmax
!   rad(ir) = rr1 * dexp(al*(ir-1))
! end do

!Determine the maximum number of points in the grid ---
 rmax  = 6.0_dp
 amesh = 0.01_dp
 mmax  = int(rmax/amesh)
! if(mod(mmax,2) .eq. 0) mmax = mmax + 1

!Print core charge info, for compatibility with psp8
 rchrg  = zero
 fchrg  = zero
 if (has_nlcc) then
   rchrg = amesh * (mmax - 2)
!  PSML does not store fchrg for now but we know we have core corrections,
!  then let's set it arbitrarily to 1.0
   fchrg = one
 else
   write(message, '(a)' ) '- psp9in: No XC core correction.'
   call wrtout(std_out,message,'COLL')
 end if
 write(message, '(3f20.14,t64,a)' ) rchrg,fchrg,zero,'rchrg,fchrg,qchrg'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

!Do we have a valence charge?
 call ps_ValenceConfiguration_Get(psxml, nshells=nshells)
 has_tvale = (nshells > 0)

! Compute the valence charge of the reference configuration used to 
! generate the pseudopotential
 chgvps = 0.0_dp
! Loop on all the shells included in the valence
 do il = 1, nshells
!  Sum the corresponding occupation of each shell
!  FIXME: What if there is spin?
   call ps_ValenceShell_Get(psxml, il, n=val_n, l=val_l, occupation=val_occ)
   chgvps = chgvps + val_occ
   write(std_out,*)' psp9in : n, l, occupation = ',   &
&   val_n, val_l, val_occ
 end do

!DEBUG
!write(std_out,*)' psp9in : atomic number'
!write(std_out,*)' psp9in :   z = ', z
!write(std_out,*)' psp9in : valence charge of the reference configuration'
!write(std_out,*)' psp9in :   chgvps = ', chgvps
!write(std_out,*)' psp9in : nominal valence charge'
!write(std_out,*)' psp9in :   zval = ', zval
!write(std_out,*)' psp9in :   mqgrid_vl = ', mqgrid_vl
!write(std_out,*)' psp9in : parameters to define the points of the grid'
!write(std_out,*)' psp9in :   amesh = ', amesh
!write(std_out,*)' psp9in :   rmax = ', rmax
!write(std_out,*)' psp9in :   mmax = ', mmax
!ENDDEBUG

! TODO: should be simple to average these and get difference for SREL+SOC,
! but also the Ekb etc...
 call ps_NonlocalProjectors_Filter(psxml, set=SET_DOWN, number=np_dn)
 call ps_NonlocalProjectors_Filter(psxml, set=SET_LJ, number=np_lj)
 call ps_NonlocalProjectors_Filter(psxml, set=SET_NONREL, number=np_nr)
 call ps_NonlocalProjectors_Filter(psxml, set=SET_SO, number=np_so)
 call ps_NonlocalProjectors_Filter(psxml, set=SET_SREL, number=np_sr)
 call ps_NonlocalProjectors_Filter(psxml, set=SET_UP, number=np_up)
 if (np_lj > 0) then
   message = 'For the moment LJ format projectors are not supported; SREL + SO is the internal abinit format'
   MSG_BUG(message)
 end if

 if (np_up > 0 .or. np_dn > 0) then
   write (message,'(3a)') 'For the moment separate spin up and down format projectors are not supported;',ch10,&
&   ' spin average is the internal abinit format'
   MSG_BUG(message)
 end if

!--------------------------------------------------------------------

!Initialize array indlmn giving l,m,n,lm,ln,s for i=lmn
 if(pspso==2) then
   nso=2
 else
   nso=1
 end if

!Find the number of projectors per angular momentum shell
 nproj(:)=0
 if (np_nr > 0) then
   call ps_NonlocalProjectors_Filter(psxml, set=SET_NONREL, indexes=idx_sr)
   do iproj = 1, np_nr
     call ps_Projector_Get(psxml, idx_sr(iproj), l=il)
     nproj(il+1) = nproj(il+1) + 1
   end do
 else
   if (np_sr > 0) then
     call ps_NonlocalProjectors_Filter(psxml, set=SET_SREL, indexes=idx_sr)
     do iproj = 1, np_sr
       call ps_Projector_Get(psxml, idx_sr(iproj), l=il)
       nproj(il+1) = nproj(il+1) + 1
     end do
   else ! this should not happen
     MSG_BUG('Your psml potential should have either scalar- or non- relativistic projectors')
   end if
 end if

 write(message, '(a,5i6)' ) '     nproj',nproj(1:lmax+1)
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,  message,'COLL')

 irelt = 0
 if (nso == 2) then
   call ps_NonlocalProjectors_Filter(psxml, set=SET_SO, indexes=idx_so)
   do iproj = 1, np_so
     call ps_Projector_Get(psxml, idx_so(iproj), l=il)
     nproj(il+lmax+2) = nproj(il+lmax+2) + 1
     irelt = 1
   end do
 end if

 pspindex=0;iln=0;indlmn(:,:)=0
 do nn=1,nso
   do ipsang=1+(nn-1)*(lmax+1),nn*lmax+1
     ll=ipsang-(nn-1)*lmax-1
     if (nproj(ipsang)>0) then
       do kk=1,nproj(ipsang)
         iln=iln+1
         do mm=1,2*ll*useylm+1
           pspindex=pspindex+1
           indlmn(1,pspindex)=ll                      ! l angular momentum channel
           indlmn(2,pspindex)=mm-ll*useylm-1          ! hash of position in m
           indlmn(3,pspindex)=kk                      ! index of projector
           indlmn(4,pspindex)=ll*ll+(1-useylm)*ll+mm  ! hash of position in l(l+1) array
           indlmn(5,pspindex)=iln                     ! absolute index of l, n disregarding m values
           indlmn(6,pspindex)=nn                      ! spin orbit index!!! NOT the n shell index
         end do
       end do
     end if
   end do
 end do

! Determine whether the atomic calculation to generate the pseudopotential
! is relativistic or not

!DEBUG
!write(std_out,*)' psp9in : pseudopotential generation relativity ', ps_Relativity(psxml) 
!write(std_out,*)' psp9in : SOC pseudopotential? (1=yes, 0 =no) '
!write(std_out,*)' psp9in : irelt = ', irelt
!write(ab_out,*)' psp9in : irelt = ', irelt
!ENDDEBUG

!Can now allocate grids, potentials and projectors
 ABI_ALLOCATE(rad,(mmax))
 ABI_ALLOCATE(vloc,(mmax))
 ABI_ALLOCATE(vpspll,(mmax,lnmax))

!Feb 2015: shifted to Hamann grid for convenience - libpsml interpolates anyway
 do ir=1,mmax
   rad(ir) = amesh * (ir - 1)
 end do
!! DEBUG
! do ir = 2, mmax
!   write(std_out,'(i5,f20.12)')ir, rad(ir)
! end do
!! ENDDEBUG
!---
 write(ab_out, '(a,i5,es16.6,es16.6)')'  psp9in : mmax, amesh, rad(mmax) = ', mmax, amesh, rad(mmax)
 
!Check that rad grid is linear starting at zero
 amesh=rad(2)-rad(1)
 damesh=zero
 do irad=2,mmax-1
   damesh=max(damesh,abs(rad(irad)+amesh-rad(irad+1)))
 end do
 if(damesh>tol8 .or. rad(1)/=zero) then
   write(message, '(5a)' )&
&   'Pseudopotential input file requires linear radial mesh',ch10,&
&   'starting at zero.',ch10,&
&   'Action: check your pseudopotential input file.'
   MSG_ERROR(message)
 end if

!Take care of the non-linear core corrections
!----------------------------------------------------------------------------
! xcccrc           : XC core correction cutoff radius (bohr)
!                    It is defined as the radius where the pseudo-core
!                    charge density becomes zero
!                    (here we have set up a tolerance of 1.d-12).

 rmatch = zero
 nders  = 0
 if (has_nlcc) then

!    In Abinit, at least for the Troullier-Martins pseudopotential,
!    the pseudocore charge density and its derivatives (xccc1d)
!    are introduced in a linear grid.
!    This grid is normalized, so the radial coordinates run between
!    from 0 and 1 (from 0 to xcccrc, where xcccrc is the radius
!    where the pseudo-core becomes zero).

   call ps_CoreCharge_get(psxml, rc=rmatch, nderivs=nders)
   write (message,'(1X,A,A,5X,A,1X,F8.3,A,5X,A,I8,A)') &
&   "Reading pseudocore charge",ch10, &
&   "- matching radius:",rmatch,ch10, &
&   "- number of continuous derivatives",nders,ch10
   call wrtout(std_out,message,'COLL')

!Get core charge function and derivatives, if needed
   if(fchrg>1.0d-15)then
     call psp9cc(psxml,mmax,n1xccc,rad,rchrg,xccc1d)
!  The core charge function for pspcod=9
!  becomes zero beyond rchrg. Thus xcccrc must be set
!  equal to rchrg.
     xcccrc=rchrg
   else
     xccc1d(:,:) = zero
     xcccrc = zero
     fchrg = zero
     qchrg = zero
   end if

   maxrad = rad(mmax)

 end if ! has_nlcc

!!   DEBUG
!    write(std_out,*)' xcccrc = ', xcccrc, rchrg
!    write(std_out,*)
!    write(std_out,*) '# psp8in NLCC data ', n1xccc, xcccrc
!    do ii = 1, n1xccc
!    write(std_out,'(7e20.8)')xcccrc*(ii-1.d0)/(n1xccc-1.d0),xccc1d(ii,1),&
! &         xccc1d(ii,2),xccc1d(ii,3),xccc1d(ii,4),xccc1d(ii,5),xccc1d(ii,6)
!    enddo
!    write(std_out,*)
!    stop
!!   ENDDEBUG


!--------------------------------------------------------------------
!Carry out calculations for local (lloc) pseudopotential.
!Obtain Fourier transform (1-d sine transform)
!to get q^2 V(q).

!Read and process vlocal:
!The local potential is given by a <radfunc> element under the <local-potential>
!element.
!After reading, this is a copy of the treatment to the 
!local part carry out in psp8 
!i.e. (as in Hamann pseudopotential)
!
!Read the local component of the pseudopotential
 vloc = zero
 do ir = 1, mmax
   vloc(ir) = ps_LocalPotential_Value(psxml, rad(ir))
 end do 

 call psp8lo(amesh,epsatm,mmax,mqgrid,qgrid,&
& vlspl(:,1),rad,vloc,yp1,ypn,zion)

!Fit spline to q^2 V(q) (Numerical Recipes subroutine)
 ABI_ALLOCATE(work_spl,(mqgrid))
 call spline (qgrid,vlspl(:,1),mqgrid,yp1,ypn,work_spl)
 vlspl(:,2)=work_spl(:)
 ABI_DEALLOCATE(work_spl)

!!  DEBUG         
! write(std_out,*)'# Vlocal = '
! write(std_out,*)' amesh  = ', amesh
! write(std_out,*)' epsatm = ', epsatm
! write(std_out,*)' mmax   = ', mmax  
! write(std_out,*)' mqgrid = ', mqgrid
! do ir = 1, mqgrid
!   write(std_out,*)'   qgrid = ', ir, qgrid(ir)
! enddo
! do ir = 1, mqgrid
!   write(std_out,'(a,i5,2f20.12)')'   iq, vlspl = ', ir, vlspl(ir,1), vlspl(ir,2)
! enddo
! write(std_out,*)
! do ir = 1, mmax
!   write(std_out,*)'   rad   = ', rad(ir), vloc(ir)
! enddo
! write(std_out,*)
! write(std_out,*)' yp1    = ', yp1
! write(std_out,*)' ypn    = ', ypn
! write(std_out,*)' zion   = ', zion
! stop
!!  ENDDEBUG      


!--------------------------------------------------------------------
!Take care of non-local part

!Zero out all Kleinman-Bylander energies to initialize
 do ii = 1, lmnmax ! loop over all possible projectors
   if (indlmn(6,ii) == 1) then
     call ps_Projector_Get(psxml, idx_sr(indlmn(5,ii)), ekb=ekb(indlmn(5,ii)))
   else if (indlmn(6,ii) == 2) then
     call ps_Projector_Get(psxml, idx_so(indlmn(5,ii)), ekb=ekb(indlmn(5,ii)))
   end if
 end do

!Read the KB projectors from the PSML file
!Note than in the PSML file the radial part of the projector is stored,
!while Abinit expects the radial part times the radii.
!We have to multiply by r after reading it.
!Note than in Hamann's format (psp8), Abinit directly reads r * radial_part_KB
 vpspll = zero
 do ii = 1, lmnmax
   if (indlmn(6,ii) == 1) then
     do ir = 1, mmax
       vpspll(ir, indlmn(5,ii)) = ps_Projector_Value(psxml, idx_sr(indlmn(5,ii)), rad(ir))
       vpspll(ir, indlmn(5,ii)) = rad(ir) * vpspll(ir, indlmn(5,ii))
     end do 
   else if (indlmn(6,ii) == 2) then
     do ir = 1, mmax
       vpspll(ir, indlmn(5,ii)) = ps_Projector_Value(psxml, idx_so(indlmn(5,ii)), rad(ir))
       vpspll(ir, indlmn(5,ii)) = rad(ir) * vpspll(ir, indlmn(5,ii))
     end do 
   end if
 end do 

!Allow for option of no nonlocal corrections (lloc=lmax=0)
 if (lloc==0.and.lmax==0) then
   write(message, '(a,f5.1)' ) ' Note: local psp for atom with Z=',znucl
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 else

!  ----------------------------------------------------------------------
!  Compute Vanderbilt-KB form factors and fit splines

   call psp8nl(amesh,ffspl,indlmn,lmax,lmnmax,lnmax,mmax,&
&   mqgrid,qgrid,rad,vpspll)

 end if

!!  DEBUG         
! write(std_out,*)'# KB Projectors = '
! write(std_out,*)' amesh  = ', amesh
! do ir = 1, mqgrid
!   do il = 1, lnmax
!     write(std_out,*)' iq, il, ffspl = ', ir, il, ffspl(ir,1,il), ffspl(ir,2,il)
!   enddo
! enddo
! do il = 1, lmnmax
!   write(std_out,*)' indlmn = ', il, indlmn(:,il)
! enddo
! write(std_out,*)' lmax   = ', lmax
! write(std_out,*)' lmnmax = ', lmnmax
! write(std_out,*)' lnmax  = ', lnmax
! write(std_out,*)' mmax   = ', mmax
! write(std_out,*)' mqgrid = ', mqgrid
! do ir = 1, mqgrid
!   write(std_out,*)'   qgrid = ', ir, qgrid(ir)
! enddo
! do il = 1, lnmax
!   write(std_out,*)
!   write(std_out,*)'# il = ', il
!   do ir = 1, mmax
!     write(std_out,*)'   rad   = ', rad(ir), vpspll(ir,il)
!   enddo
! enddo
! stop
!!  ENDDEBUG      

! Read pseudo valence charge in real space on the linear mesh
! and transform it to reciprocal space on a regular grid. Use vloc as workspace.
 vloc(:) = zero
 if (has_tvale) then
   do irad=1,mmax
     vloc(irad) = ps_ValenceCharge_Value(psxml,rad(irad))
     vloc(irad) = vloc(irad) / four_pi
   end do

!! DEBUG
!  do irad = 1, mmax
!    write(std_out,*)' Valence Charge  = ', rad(irad), vloc(irad)
!  enddo
!  stop
!! ENDDEBUG


   ! Check that rad grid is linear starting at zero
   amesh=rad(2)-rad(1)
   damesh=zero
   do irad=2,mmax-1
     damesh=max(damesh,abs(rad(irad)+amesh-rad(irad+1)))
   end do
   if(damesh>tol8 .or. rad(1)/=zero) then
     write(message, '(5a)' )&
&     'Pseudopotential input file requires linear radial mesh',ch10,&
&     'starting at zero.',ch10,&
&     'Action: check your pseudopotential input file.'
     MSG_ERROR(message)
   end if

   !  Evaluate spline-fit of the atomic pseudo valence charge in reciprocal space.
   call pawrad_init(mesh,mesh_size=mmax,mesh_type=1,rstep=amesh)
   call nctab_eval_tvalespl(nctab, zion, mesh, vloc, mqgrid_vl, qgrid_vl)
   call pawrad_free(mesh)
 end if

 ABI_DEALLOCATE(vpspll)
 ABI_DEALLOCATE(vloc)
 ABI_DEALLOCATE(rad)
 if (allocated(idx_sr)) then
   ABI_FREE_NOCOUNT(idx_sr)
 end if
 if (allocated(idx_so)) then
   ABI_FREE_NOCOUNT(idx_so)
 end if

 call ps_destroy(psxml)

!--------------------------------------------------------------------

#else
!Initialize some arguments, for portability at compile time
 indlmn=0 ; mmax=0 ; nproj=0
 ekb=zero ; epsatm=zero ; ffspl=zero ; qchrg=zero ; vlspl=zero ; xcccrc=zero ; xccc1d=zero

 if(.false.)write(std_out,*)filpsp ! Just to keep filpsp when HAVE_PSML is false
 if(.false.)write(std_out,*)lloc   ! Just to keep lloc when HAVE_PSML is false
 if(.false.)write(std_out,*)level  ! Just to keep level when HAVE_PSML is false
 if(.false.)write(std_out,*)lmax   ! Just to keep lmax when HAVE_PSML is false
 if(.false.)write(std_out,*)mpsang ! Just to keep mpsang when HAVE_PSML is false
 if(.false.)write(std_out,*)nctab  ! Just to keep nctab when HAVE_PSML is false
 if(.false.)write(std_out,*)pspso  ! Just to keep pspso when HAVE_PSML is false
 if(.false.)write(std_out,*)qgrid  ! Just to keep qgrid when HAVE_PSML is false
 if(.false.)write(std_out,*)qgrid_vl ! Just to keep qgrid_vl when HAVE_PSML is false
 if(.false.)write(std_out,*)useylm ! Just to keep useylm when HAVE_PSML is false
 if(.false.)write(std_out,*)zion   ! Just to keep zion when HAVE_PSML is false
 if(.false.)write(std_out,*)znucl  ! Just to keep znucl when HAVE_PSML is false
#endif

end subroutine psp9in
!!***
