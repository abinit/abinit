!!****m* ABINIT/m_evdw_wannier
!! NAME
!! m_evdw_wannier
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2010-2020 ABINIT group (CE, TR, AR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_evdw_wannier

 use defs_basis
 use m_abicore
 use m_errors

 use m_special_funcs,   only : abi_derf
 use m_numeric_tools,   only : simpson_int
 use m_geometry,        only : xcart2xred, xred2xcart

 implicit none

 private
!!***

 public :: evdw_wannier
!!***

contains
!!***

!!****f* ABINIT/evdw_wannier
!! NAME
!! evdw_wannier
!!
!! FUNCTION
!!  FIXME: Evaluates the van der Waals correlation energy using maximally
!!         localized Wannier functions (MLWF) as proposed by:
!!         P. L. Silvestrelli in PRL 100:053002 (2008) [[cite:Sivestrelli2008]] vdw_xc=10 and
!!         A. Ambrosetti and P. L. Silvestrelli in PRB 85:073101 (2012) [[cite:Ambrosetti2012]] vdw_xc=11.
!!         P. L. Silvestrelli in J.Chem.Phys. 139:054106 (2013) [[cite:Silvestrelli2013]] vdw_xc=14.
!!
!! INPUTS
!!   nsppol          = Spin polarization.
!!   nwan(nsppol)    = Total number of MLWF in the system per spin component.
!!   origmwan        = max[nwan(nsppol)] from mlwfovlp.F90.
!!   tdocc_wan       = MLWFs occupation matrix diagonal terms
!!   vdw_nfrag       = Number of vdW interating fragments in the unit cell.
!!   vdw_supercell(3)     = Distance along each rprimd components for
!!                          which vdW interactions between MLWF will be taken into account.
!!   vdw_typfrag(natom)   = Fragment to which each atom belongs to.
!!   vdw_xc               = vdW-WF version.
!!   rprimd               = Real space primitive translations.
!!   wann_centres(3,origmwan,nsppol) = The centers of MLWFs  in a.u.
!!   wann_spreads(origmwan,nsppol)   = Spread of the MLWFs, in Ang**2. (from wannier90).
!!   xcart           = Coordinates of unit cell atoms in atomic units.
!!
!! OUTPUT
!!   csix(origmwan,origmwan,nsppol,nsppol) = dispersion coefficient between each pair of MLWF.
!!   corrvdw           = van der Waals correction to the energy.
!!
!! PARENTS
!!      m_mlwfovlp
!!
!! CHILDREN
!!
!! SOURCE

 subroutine evdw_wannier(csix,corrvdw,origmwan,natom,nsppol,orignwan,tdocc_wan,vdw_nfrag,&
& vdw_supercell,vdw_typfrag,vdw_xc,rprimd,wann_centres,wann_spreads,xcart)

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: origmwan,nsppol,natom,orignwan(nsppol)
 integer , intent(in)  :: vdw_nfrag,vdw_supercell(3),vdw_typfrag(natom),vdw_xc
 real(dp), intent(in)  :: rprimd(3,3),wann_centres(3,origmwan,nsppol),wann_spreads(origmwan,nsppol)
 real(dp), intent(in)  :: xcart(3,natom)
 real(dp), intent(out) :: corrvdw
 real(dp), intent(out) :: csix(origmwan,origmwan,nsppol,nsppol)
 real(dpc), intent(in) :: tdocc_wan(origmwan,nsppol)

!Local variables-------------------------------
 integer  :: ier,igr,icx,icy,icz,ii,inx,iny,inz,isppol,iwan
 integer  :: jwan,jj,ll,mm,mwan,nc,ngr,tmp_mwan,mwan_half
 integer, allocatable:: amagr(:,:,:),inwan(:,:),nw(:,:),nwan(:),npwf(:),ord(:,:)
 integer, allocatable:: tmp_nwan(:)
 real(dp) :: dnrm2,fij,rij,rij_c(3),fu,shift,erfValue
 real(dp), parameter :: a = 20.d0 !Parameter related to the damping function.
 real(dp), parameter :: gama = 4.5d0/(sqrt3**3) !alpha=gama*S**3.
 real(dp), parameter :: gama1 = 0.88d0 !alpha=gama*S**3.
 real(dp), parameter :: zeta = 1.30d0 !polar=zeta*(Z/omega**2).
 real(dp), parameter :: beta = 1.39d0 !    .
 real(dp), allocatable:: amawf(:,:),amaspr(:),amaocc(:),dcenters(:,:,:),rc(:,:)
 real(dp), allocatable:: tmp_cent(:,:,:),tmp_spr(:,:),tmp_occ(:,:)
 real(dp), allocatable:: rv(:,:),wanncent(:,:,:),wannspr(:,:),wc_rec(:,:,:),xi(:,:)
 real(dp), allocatable:: c_QHO(:,:),Tij_dip(:,:),polar(:),omega(:),eigv(:),zhpev2(:)
 real(dpc), allocatable :: newocc_wan(:,:)
 complex(dpc), allocatable :: eigvec(:,:),matrx(:),zhpev1(:)
 character(len=500) :: message                   ! to be uncommented, if needed
! *************************************************************************

!Determining presence p-like MLWFs see J.Chem.Phys.135:154105 (2011) [[cite:Andrinopoulos2011]]
 ABI_MALLOC(npwf,(nsppol))
 ABI_MALLOC(inwan,(origmwan,nsppol))
 ABI_MALLOC(nwan,(nsppol))

 ll = 0
 npwf(:) = 0
 inwan(:,:) = 0
 do jj=1,nsppol
   do iwan=1,orignwan(jj)
     if(tdocc_wan(iwan,jj)*nsppol<=1.50d0) then
       npwf(jj) = npwf(jj) + 1
       ll = ll+1
       inwan(ll,jj) = iwan
     end if
   end do
 end do

 write(std_out,*) ch10,'Number of p-like MLWFs per spin pol:',ch10
 write(std_out,*) (npwf(ii),ii=1,nsppol), ch10

 mwan=origmwan+(sum(npwf(:))) !two new MLWFs per p-like MLWF
 nwan(:)=orignwan(:)+npwf(:)


 ABI_MALLOC(wanncent,(3,mwan,nsppol))
 ABI_MALLOC(wannspr,(mwan,nsppol))
 ABI_MALLOC(wc_rec,(3,mwan,nsppol))
 ABI_MALLOC(ord,(mwan,nsppol))
 ABI_MALLOC(newocc_wan,(mwan,nsppol))

 wanncent(:,:,:) = zero
 wannspr(:,:) = zero
 wc_rec(:,:,:) = zero
 newocc_wan(:,:) = zero
 ord(:,:) = zero

!The vdW correction is calculated in atomic units:
 do ii=1,nsppol
   do iwan=1,orignwan(ii)
!    converting to bohr**2 and then squared
     wanncent(:,iwan,ii)=wann_centres(:,iwan,ii)/Bohr_Ang
!    write(std_out,*) "spread of WF",i, "=", wann_spreads(i)
     wannspr(iwan,ii)=sqrt(wann_spreads(iwan,ii)/Bohr_Ang**2)
     newocc_wan(iwan,ii)=tdocc_wan(iwan,ii)
   end do
 end do

!write(std_out,*) 'Number of MLWFs:',ch10
!do ii=1,nsppol
!write(std_out,*) 'nsppol=',ii, 'nwan(nsppol)=',nwan(nsppol),ch10
!end do

 write(std_out,*) 'Original Wannier centres and spreads:',ch10
 do ii=1,nsppol
   write(std_out,*) 'nsppol=',ii,ch10
   do iwan=1,orignwan(ii)
     write(std_out,*) (wanncent(jj,iwan,ii),jj=1,3), wannspr(iwan,ii),ch10
   end do
 end do

!Translate MLWFs to the original unit cell if vdw_nfrag > 0 :

 if(vdw_nfrag>0)then
   do jj=1,nsppol
     call xcart2xred(orignwan(jj),rprimd,wanncent(:,1:orignwan(jj),jj), &
&     wc_rec(:,1:orignwan(jj),jj))
!    got centers in reduced coor
     do iwan=1,orignwan(jj)
       do ii=1,3
         if(wc_rec(ii,iwan,jj)<zero) then
           shift=REAL(CEILING(ABS(wc_rec(ii,iwan,jj))),dp)
           wc_rec(ii,iwan,jj) = wc_rec(ii,iwan,jj)+shift
         end if
         if(wc_rec(ii,iwan,jj)>one) then
           shift=-REAL(INT(wc_rec(ii,iwan,jj)),dp)
           wc_rec(ii,iwan,jj) = wc_rec(ii,iwan,jj)+shift
         end if
       end do
     end do
     call xred2xcart(orignwan(jj),rprimd,wanncent(:,1:orignwan(jj),jj), &
&     wc_rec(:,1:orignwan(jj),jj))
   end do

!  ====================================================================

   write(std_out,*) ch10,'Wannier centres translated to unit cell and spr:',ch10
   do jj=1,nsppol
     write(std_out,*) 'nsppol=',jj,ch10
     do iwan=1,orignwan(jj)
       write(std_out,*) (wanncent(ii,iwan,jj),ii=1,3), wannspr(iwan,jj)
     end do
   end do
 end if !vdw_nfrag>0

!Spliting of p-like into 2 s-like MLWFs
!Eqs. (22) and (23) of J.Chem.Phys.135:154105 (2011) [[cite:Andrinopoulos2011]]

 if ( any (npwf(:)/=0) ) then

   write(std_out,*) 'Indexes of p-like MLWFs and its spin:'

   do isppol=1,nsppol
     do jj=1,npwf(isppol)

       write(std_out,*) inwan(jj,isppol),isppol

       wanncent(1:2,orignwan(isppol)+jj,isppol) = wanncent(1:2,inwan(jj,isppol),isppol)

       wanncent(3,orignwan(isppol)+jj,isppol) = wanncent(3,inwan(jj,isppol),isppol)  &
&       + 15.d0*wannspr(inwan(jj,isppol),isppol) / (eight*sqrt(30.d0))

       wanncent(3,inwan(jj,isppol),isppol) = wanncent(3,inwan(jj,isppol),isppol)  &
&       - 15.d0*wannspr(inwan(jj,isppol),isppol) / (eight*sqrt(30.d0))

       wannspr(orignwan(isppol)+jj,isppol) = seven*wannspr(inwan(jj,isppol),isppol) / (eight*sqrt2)

       wannspr(inwan(jj,isppol),isppol) = seven*wannspr(inwan(jj,isppol),isppol) / (eight*sqrt2)

       newocc_wan(orignwan(isppol)+jj,isppol) = tdocc_wan(inwan(jj,isppol),isppol) / two

       newocc_wan(inwan(jj,isppol),isppol) = tdocc_wan(inwan(jj,isppol),isppol) / two

     end do
   end do

   write(std_out,*) ch10,'Wannier centres and spreads after splitting of p-like MLWFs:',ch10
   do isppol=1,nsppol
     write(std_out,*) 'nsppol=',isppol,ch10
     do iwan=1,nwan(isppol)
       write(std_out,*) (wanncent(jj,iwan,isppol),jj=1,3), wannspr(iwan,isppol)
     end do
   end do

 end if ! any(npwf(:)/=0)

!Asign each MLWFs to one fragment, the same as their nearest atom:

 call order_wannier(mwan,natom,nwan,nsppol,ord,vdw_typfrag,wanncent,xcart)

 write(std_out,*) ch10,'Wannier centres and fragments',ch10
 do ll=1,abs(vdw_nfrag)
   write(std_out,*) 'MLWF centers in fragment',ll,ch10
   do jj=1,nsppol
     do iwan=1,nwan(jj)
       if (ord(iwan,jj)==ll) then
         write(std_out,*) 'X', (Bohr_Ang*wanncent(ii,iwan,jj),ii=1,3),iwan,jj
       end if
     end do
   end do
 end do

 write(std_out,*) ch10,'Occupation Matrix diagonal terms:',ch10
 do ll=1,abs(vdw_nfrag)
   write(std_out,*) 'For MLWF centers in fragment',ll,ch10
   do jj=1,nsppol
     do iwan=1,nwan(jj)
       if (ord(iwan,jj)==ll) then
         write(std_out,*) newocc_wan(iwan,jj),ch10
       end if
     end do
   end do
 end do

!Amalgamation of close MLWFs, see J.Chem.Phys.135:154105 (2011) [[cite:Andrinopoulos2011]]

 if (all(npwf(:)==0).and.vdw_xc/=14) then !amalgamation is done only if no p-like

   mwan_half=mwan/2
   ABI_MALLOC(amagr,(mwan,nsppol,mwan_half))
   ABI_MALLOC(nw,(nsppol,mwan_half))
   nw=0

   call amalgam(amagr,ngr,nsppol,nw,mwan,ord,nwan,vdw_nfrag,wanncent,wannspr)

!  Calculating amalgamated centres, spreads and occupancies if any:

   if( any(nw(:,:) /= 0) ) then

     ABI_MALLOC(amawf,(3,ngr))
     ABI_MALLOC(amaspr,(ngr))
     ABI_MALLOC(amaocc,(ngr))

     amawf(:,:) = 0
     amaspr(:) = 0
     amaocc(:) = 0

     do igr = 1 , ngr
       do isppol =  1 , nsppol
         do ii = 1 , nw(isppol,igr)

           amawf(:,igr) =  amawf(:,igr) + wanncent(:,amagr(ii,isppol,igr),isppol)
           amaspr(igr)  =  amaspr(igr) + wannspr(amagr(ii,isppol,igr),isppol)
           amaocc(igr)  =  amaocc(igr) + newocc_wan(amagr(ii,isppol,igr),isppol)

         end do
       end do

       amawf(:,igr) = amawf(:,igr) / real(sum(nw(1:nsppol,igr)),dp )
       amaspr(igr)  = amaspr(igr) / real(sum(nw(1:nsppol,igr)),dp )

     end do

     write(std_out,*) ch10,'Amalgamated MLWFs Centres, Spreads and Occupancies:',ch10
     do igr = 1 , ngr
       write(std_out,*) (amawf(ii,igr),ii=1,3),amaspr(igr),amaocc(igr)
     end do

!    Redefining centres, spreads and occps arrays:
     ABI_MALLOC(tmp_nwan,(nsppol))

     tmp_nwan(:) = nwan(:) - sum(nw(:,1:ngr))
     tmp_mwan = maxval(tmp_nwan(:))

     ABI_MALLOC(tmp_cent,(3,tmp_mwan,nsppol))
     ABI_MALLOC(tmp_spr,(tmp_mwan,nsppol))
     ABI_MALLOC(tmp_occ,(tmp_mwan,nsppol))

     tmp_cent(:,:,:) = zero
     tmp_spr(:,:) = zero
     tmp_occ(:,:) = zero

     do isppol = 1 , nsppol
       ii = 0
       do iwan = 1 , nwan(isppol)

         if ( any(amagr(:,isppol,:) == iwan) ) cycle

         ii = ii + 1
         tmp_cent(:,ii,isppol) = wanncent(:,iwan,isppol)
         tmp_spr(ii,isppol) = wannspr(iwan,isppol)
         tmp_occ(ii,isppol) = newocc_wan(iwan,isppol)

       end do
     end do

!    Redefining wanncent, wannspr, newocc_wan:
!    Even if amalgamation occurs with MLWFs of different spins
!    the new WF are gathered with isppol=1 functions...

     nwan(1) = nwan(1) - sum(nw(1,1:ngr)) + ngr

     if (nsppol == 2) then
       nwan(2) = nwan(2) - sum(nw(2,1:ngr))
     end if

     mwan = maxval(nwan(:))

     do isppol = 1 , nsppol
       do iwan = 1 , tmp_nwan(isppol)

         wanncent(:,iwan,isppol) = tmp_cent(:,iwan,isppol)
         wannspr(iwan,isppol) = tmp_spr(iwan,isppol)
         newocc_wan(iwan,isppol) = tmp_occ(iwan,isppol)

       end do
     end do

     do igr = 1 , ngr

       wanncent(:,tmp_nwan(1)+igr,1) = amawf(:,igr)
       wannspr(tmp_nwan(1)+igr,1) = amaspr(igr)
       newocc_wan(tmp_nwan(1)+igr,1) = amaocc(igr)

     end do

!    Ordering again:
!    Asign each MLWFs to one fragment, the same as their nearest atom:

     call order_wannier(mwan,natom,nwan,nsppol,ord,vdw_typfrag,wanncent,xcart)


     write(std_out,*) ch10,'Full set of Wannier functions and spreads'
     write(std_out,*) 'after both splitting of p-like WFs and amalgamation',ch10

     do ll=1,abs(vdw_nfrag)
       write(std_out,*) 'MLWF centers and spreads in fragment',ll,ch10
       do jj=1,nsppol
         do iwan=1,nwan(jj)
           if (ord(iwan,jj)==ll) then
             write(std_out,*) 'X', (Bohr_Ang*wanncent(ii,iwan,jj),ii=1,3),Bohr_Ang*wannspr(iwan,jj)
           end if
         end do
       end do
     end do

   end if ! any(nw(:,:) /= 0)
 end if ! all((npwf(:)==0).and.vdw_xc/=14)

!vdW-WF VERSION 1

 if(vdw_xc==10) then

   ABI_MALLOC(dcenters,(3,mwan,nsppol))
   ABI_MALLOC(rc,(mwan,nsppol))
   ABI_MALLOC(rv,(mwan,nsppol))
!  Calculate intermediate quantities
   do jj=1,nsppol
     do iwan=1, nwan(jj)
       rc(iwan,jj)= three*(0.769d0+half*dlog(wannspr(iwan,jj)))
!      rv(iwan,jj)= (1.475d0-half_sqrt3*dlog(wannspr(iwan,jj)))*wannspr(iwan,jj)
!      r_v suggested in JPhysChemA 113:5224 [[cite:Silvestrelli2009]]
       rv(iwan,jj)= (rc(iwan,jj)*wannspr(iwan,jj))/sqrt3 
     end do
   end do
   corrvdw=0.0d0  !Initializing the vdW correction energy.

   do ii=1,nsppol
     do jj=1,nsppol
       do iwan=1,nwan(ii)
         do jwan=1,nwan(jj)

           call getFu(wannspr(iwan,ii),wannspr(jwan,jj),rc(iwan,ii),rc(jwan,jj),&
&           newocc_wan(iwan,ii),newocc_wan(jwan,jj),fu)

           csix(iwan,jwan,ii,jj)=( ( ((wannspr(iwan,ii))**1.5d0)*&
&           (wannspr(jwan,jj)**three))/(two*(three**1.25d0) ) )*fu

         end do
       end do
     end do
   end do

!  if (nsppol == 1) then
!  csix(:,:,:,:)=sqrt2*csix(:,:,:,:)  !For non spin polarized systems
!  end if


!  DEBUG
   write(std_out,*) ch10,'C6ij coefficients matrix',ch10
   do ii=1,nsppol
     do jj=1,nsppol
       do iwan=1,nwan(ii)
         write(std_out,*) (csix(iwan,jwan,ii,jj),jwan=1,nwan(jj)),ch10
       end do
     end do
   end do
!  END DEBUG

!  test   k=0
   do ii=1,nsppol
     do iwan=1,nwan(ii)
       do inx=-abs(vdw_supercell(1)),abs(vdw_supercell(1))
         do iny=-abs(vdw_supercell(2)),abs(vdw_supercell(2))
           do inz=-abs(vdw_supercell(3)),abs(vdw_supercell(3))
             do jj=1,nsppol
               do jwan=1,nwan(jj)

                 if(inx==0.and.iny==0.and.inz==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
!                This avoids intrafragment vdW interactions.
                 if(vdw_supercell(1)<=0.and.inx==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
                 if(vdw_supercell(2)<=0.and.iny==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
                 if(vdw_supercell(3)<=0.and.inz==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
!                Last three conditions allow proper treatment of layered systems.

                 dcenters(:,jwan,jj) = (real(inx,dp))*rprimd(:,1)+(real(iny,dp))*rprimd(:,2)+&
&                 (real(inz,dp))*rprimd(:,3)+wanncent(:,jwan,jj)
                 rij=sqrt((dcenters(1,jwan,jj)-wanncent(1,iwan,ii))**2+&
&                 (dcenters(2,jwan,jj)-wanncent(2,iwan,ii))**2+&
&                 (dcenters(3,jwan,jj)-wanncent(3,iwan,ii))**2)

                 fij=one/(one+exp(-a*(rij/(rv(iwan,ii)+rv(jwan,jj))-one))) !Damping function.

                 corrvdw = corrvdw - csix(iwan,jwan,ii,jj)*fij/(two*(rij**6)) !making the sum of eq(4) of
!                JPhysChemA 113:5224-5234 [[cite:Silvestrelli2009]]. Each term is divided by two because
!                we are counting twice within the unit cell, also the
!                interactions with neighbor cells are properly acounted for in
!                this way.

!                write(std_out,*) 'i=',iwan, 'j=',jwan, 'C6ij=', csix(iwan,jwan)
!                write(std_out,*) 'inx=',inx, "iny=",iny, "inz=",inz, "Evdw=",&
!                & -(csix(iwan,jwan)*fij/(two*rij**6))*Ha_ev*ten**3
!                write(std_out,*) 'rnl=',rnl
               end do
             end do
           end do
         end do
       end do
     end do
   end do

   ABI_FREE(dcenters)
   ABI_FREE(rc)
   ABI_FREE(rv)

   write(message, '(2a,i2,2a,f12.6,2a,f12.6,a)' )ch10,&
&   ' vdw_xc : ',10,ch10,&
&   ' van der Waals correction(Ha):',   corrvdw,ch10,&
&   ' van der Waals correction(eV):',   corrvdw*Ha_ev,ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

 end if

!vdW-WF VERSION 2: Phys. Rev. B. 85:073101 (2012) [[cite:Ambrosetti2012]]

 if (vdw_xc==11) then

   ABI_MALLOC(dcenters,(3,mwan,nsppol))
   ABI_MALLOC(rv,(mwan,nsppol))
   ABI_MALLOC(xi,(mwan,nsppol))

!  Calculate intermediate quantities
   do jj=1,nsppol
     do iwan=1, nwan(jj)
       rv(iwan,jj)= ( (1.20d0/Bohr_Ang)*wannspr(iwan,jj) )/sqrt3
       write(std_out,*) 'rv(iwan,jj)=',rv(iwan,jj),ch10
     end do
   end do

!  C6 coefficients between WF
   csix(:,:,:,:) = 0.0d0
   corrvdw = 0.0d0

   call ovlp_wann(mwan,nwan,nsppol,ord,wanncent,wannspr,xi)

!  DEBUG
   write(std_out,*)ch10,'xi(iwan,isspol)=',ch10
   do jj=1,nsppol
     write(std_out,*) (xi(iwan,jj),iwan=1,nwan(jj))
   end do
!  END DEBUG

   do ii=1,nsppol
     do jj=1,nsppol
       do iwan=1,nwan(ii)
         do jwan=1,nwan(jj)

           csix(iwan,jwan,ii,jj)=onehalf*( (wannspr(iwan,ii)*wannspr(jwan,jj))**three )*&
&           ((xi(iwan,ii)*xi(jwan,jj))*gama**onehalf)/( sqrt(xi(iwan,ii))*&
&           wannspr(iwan,ii)**onehalf + sqrt(xi(jwan,jj))*wannspr(jwan,jj)**onehalf )

         end do
       end do
     end do
   end do

!  if (nsppol == 1) then
!  csix(:,:,:,:)=sqrt2*csix(:,:,:,:)  !For non spin polarized systems
!  end if

!  DEBUG
   write(std_out,*) ch10,'C6ij coefficients:',ch10
   do ii=1,nsppol
     do jj=1,nsppol
       do iwan=1,nwan(ii)
         write(std_out,*) (csix(iwan,jwan,ii,jj),jwan=1,nwan(jj))
       end do
     end do
   end do
!  END DEBUG
   do ii=1,nsppol
     do iwan=1,nwan(ii)
       do inx=-abs(vdw_supercell(1)),abs(vdw_supercell(1))
         do iny=-abs(vdw_supercell(2)),abs(vdw_supercell(2))
           do inz=-abs(vdw_supercell(3)),abs(vdw_supercell(3))
             do jj=1,nsppol
               do jwan=1,nwan(jj)

                 if(inx==0.and.iny==0.and.inz==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
!                This avoids intrafragment vdW interactions.
                 if(vdw_supercell(1)<=0.and.inx==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
                 if(vdw_supercell(2)<=0.and.iny==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
                 if(vdw_supercell(3)<=0.and.inz==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
!                Last three conditions allow proper treatment of layered systems.

                 dcenters(:,jwan,jj) = (real(inx,dp))*rprimd(:,1)+(real(iny,dp))*rprimd(:,2)+&
&                 (real(inz,dp))*rprimd(:,3)+wanncent(:,jwan,jj)
                 rij=sqrt((dcenters(1,jwan,jj)-wanncent(1,iwan,ii))**2+&
&                 (dcenters(2,jwan,jj)-wanncent(2,iwan,ii))**2+&
&                 (dcenters(3,jwan,jj)-wanncent(3,iwan,ii))**2)

                 fij=one/(one+exp(-a*(rij/(rv(iwan,ii)+rv(jwan,jj))-one))) !Damping function.
!                DEBUG
!                write(std_out,*) 'f_i,j=',fij,ch10
!                END DEBUG
                 corrvdw = corrvdw - csix(iwan,jwan,ii,jj)*fij/(two*(rij**6)) !making the sum of eq(4) of
!                JPhysChemA 113:5224-5234 [[cite:Silvestrelli2009]]
               end do
             end do
           end do
         end do
       end do
     end do
   end do

   write(message, '(2a,i2,2a,f12.6,2a,f12.6,a)' )ch10,&
&   ' vdw_xc : ',11,ch10,&
&   ' van der Waals correction(Ha):',   corrvdw,ch10,&
&   ' van der Waals correction(eV):',   corrvdw*Ha_ev,ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   ABI_FREE(dcenters)
   ABI_FREE(rv)
   ABI_FREE(xi)
 end if

!vdW-WF VERSION 3 (Using the long range limit of VV10 functional)

 if(vdw_xc==12) then

   ABI_MALLOC(dcenters,(3,mwan,nsppol))
   ABI_MALLOC(rc,(mwan,nsppol))
   ABI_MALLOC(rv,(mwan,nsppol))
!  Calculate intermediate quantities
   do jj=1,nsppol
     do iwan=1, nwan(jj)
!      rc(iwan,jj)= three*(0.769d0+half*dlog(wannspr(iwan,jj))) !from Silvestrelli see above.
       rc(iwan,jj)=three*wannspr(iwan,jj) !integral cutoff
!      rv(iwan,jj)= (1.475d0-half_sqrt3*dlog(wannspr(iwan,jj)))*wannspr(iwan,jj)
       rv(iwan,jj)= wannspr(iwan,jj)*sqrt3*(0.769d0+half*dlog(wannspr(iwan,jj)))
!      r_v suggested in JPhysChemA 113:5224 [[cite:Silvestrelli2009]]
     end do
   end do
   corrvdw=0.0d0  !Initializing the vdW correction energy.

   do ii=1,nsppol
     do jj=1,nsppol
       do iwan=1,nwan(ii)
         do jwan=1,nwan(jj)

           call vv10limit(wannspr(iwan,ii),wannspr(jwan,jj),rc(iwan,ii),rc(jwan,jj),fu)

           csix(iwan,jwan,ii,jj)=(1296.0d0/( (wannspr(iwan,ii)*wannspr(jwan,jj) )**3))*fu
!          vv10limit needs revision as an error regarding occupations has been included
!          currently we are calculating it with 1 electron per MLWF and there is an error four-->two
         end do
       end do
     end do
   end do

!  if (nsppol == 1) then
!  csix(:,:,:,:)=sqrt2*csix(:,:,:,:)  !For non spin polarized systems
!  end if


!  DEBUG

   write(std_out,*) ch10,'C6ij coefficients matrix',ch10
   do ii=1,nsppol
     do jj=1,nsppol
       do iwan=1,nwan(ii)
         write(std_out,*) (csix(iwan,jwan,ii,jj),jwan=1,nwan(jj)),ch10
       end do
     end do
   end do
!  END DEBUG

!  test   k=0
   do ii=1,nsppol
     do iwan=1,nwan(ii)
       do inx=-abs(vdw_supercell(1)),abs(vdw_supercell(1))
         do iny=-abs(vdw_supercell(2)),abs(vdw_supercell(2))
           do inz=-abs(vdw_supercell(3)),abs(vdw_supercell(3))
             do jj=1,nsppol
               do jwan=1,nwan(jj)

                 if(inx==0.and.iny==0.and.inz==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
!                This avoids intrafragment vdW interactions.
                 if(vdw_supercell(1)<=0.and.inx==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
                 if(vdw_supercell(2)<=0.and.iny==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
                 if(vdw_supercell(3)<=0.and.inz==0.and.ord(jwan,jj)==ord(iwan,ii)) cycle
!                Last three conditions allow proper treatment of layered systems.

                 dcenters(:,jwan,jj) = (real(inx,dp))*rprimd(:,1)+(real(iny,dp))*rprimd(:,2)+&
&                 (real(inz,dp))*rprimd(:,3)+wanncent(:,jwan,jj)
                 rij=sqrt((dcenters(1,jwan,jj)-wanncent(1,iwan,ii))**2+&
&                 (dcenters(2,jwan,jj)-wanncent(2,iwan,ii))**2+&
&                 (dcenters(3,jwan,jj)-wanncent(3,iwan,ii))**2)

                 fij=one/(one+exp(-a*(rij/(rv(iwan,ii)+rv(jwan,jj))-one))) !Damping function.

                 corrvdw = corrvdw - csix(iwan,jwan,ii,jj)*fij/(two*(rij**6)) !making the sum of eq(4) of
!                JPhysChemA 113:5224-5234 [[cite:Silvestrelli2009]]. Each term is divided by two because
!                we are counting twice within the unit cell, also the
!                interactions with neighbor cells are properly acounted for in
!                this way.

!                write(std_out,*) 'i=',iwan, 'j=',jwan, 'C6ij=', csix(iwan,jwan)
!                write(std_out,*) 'inx=',inx, "iny=",iny, "inz=",inz, "Evdw=",&
!                & -(csix(iwan,jwan)*fij/(two*rij**6))*Ha_ev*ten**3
!                write(std_out,*) 'rnl=',rnl
               end do
             end do
           end do
         end do
       end do
     end do
   end do

   ABI_FREE(dcenters)
   ABI_FREE(rc)
   ABI_FREE(rv)

   write(message, '(2a,i2,2a,f12.6,2a,f12.6,a)' )ch10,&
&   ' vdw_xc : ',12,ch10,&
&   ' van der Waals correction(Ha):',   corrvdw,ch10,&
&   ' van der Waals correction(eV):',   corrvdw*Ha_ev,ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

 end if


!vdW-QHO-WF method.

 if(vdw_xc==14) then

! There is no need of building the full set of MLWFs corresponding to the vdw_supercell
! since the matrix elements can be computed on the fly by translating the MLWFs centers.
! The polarizability and QHO frequencies are obteined for the MLWFs in the unit cell:

   ABI_MALLOC(polar,(mwan))
   ABI_MALLOC(omega,(mwan))

   corrvdw=zero
   fu=zero

   do isppol=1,nsppol
     polar=zero
     omega=zero
     do iwan=1,nwan(isppol)
       polar(iwan)=gama1*(wannspr(iwan,isppol)**3)
! assuming Z( not zeta) is the charge of a single Wannier function, 2 for non polarized and 1 for polarized)
       omega(iwan)=sqrt(zeta*(3-nsppol)/polar(iwan))
       fu=fu+omega(iwan)
     end do
!DEBUG
     write(std_out,*) 'Unit cell non interacting QHO energy:',ch10
     write(std_out,*) (1.5d0)*fu,ch10
!ENDDEBUG

!  Total number of unit cells considered:
     nc=(2*abs(vdw_supercell(1))+1)*(2*abs(vdw_supercell(2))+1)*(2*abs(vdw_supercell(3))+1)
!DEBUG
     write(std_out,*) 'Evaluation of vdW energy from ',nc,' unit cells.',ch10

     write(std_out,*) 'VdW supercell non interacting QHO energy:',ch10
     fu=nc*fu
     write(std_out,*) (1.5d0)*fu,ch10
!ENDDEBUG
     ABI_MALLOC(c_QHO,(3*mwan*nc,3*mwan*nc))
     ABI_MALLOC(Tij_dip,(3*mwan*nc,3*mwan*nc))
     ABI_MALLOC(dcenters,(3,mwan,nsppol))

     c_QHO(:,:)=zero
     Tij_dip(:,:)=zero
     inx=0
     iny=0

   !writing matrix diagonal terms
     do ll=1,nc
       do iwan=1,nwan(isppol)
         do ii=1,3
           inx=inx+1
           c_QHO(inx,inx)=omega(iwan)*omega(iwan)
         end do
       end do
     end do

   !writing terms for interactions from each cell to the central unit cell
   ! icx, icy, icz labels the cells in the supercell defined by vdw_supercell
   ! while inx and iny label QHO matrix elements.
   ! iny should start from (nc/2)*3*mwan + 1 --> r=r_central-r_cell, and r_cell=displaced positions
   ! inx starts from 1 up to (nc/2)*3*mwan +1, this includes central cell intra-interactions.

     inx=0
     do icz=-abs(vdw_supercell(3)),0
       if (icz==0) then
         mm=0
       else
         mm=abs(vdw_supercell(2))
       end if
       do icy=-abs(vdw_supercell(2)),mm
         if (icy==0) then
           ll=0
         else
           ll=abs(vdw_supercell(1))
         end if
         do icx=-abs(vdw_supercell(1)),ll
           do iwan=1,nwan(isppol)
             do ii=1,3
               inx=inx+1
       ! loop over the MLWFs in the 'central' unit cell:
               iny=(nc/2)*3*mwan
               do jwan=1,nwan(isppol)
                 do jj=1,3
                   iny=iny+1

                   if (inx==iny) cycle !in order to avoid digonal terms which  were already computed

                   dcenters(:,iwan,isppol) = (real(icx,dp))*rprimd(:,1)+(real(icy,dp))*rprimd(:,2)+&
&                   (real(icz,dp))*rprimd(:,3)+wanncent(:,iwan,isppol)

                   rij_c = -dcenters(:,iwan,isppol)+wanncent(:,jwan,isppol)
                   rij = dnrm2(3,rij_c,1)
!                  rij=sqrt(dot_product(rij_c,rij_c))
                   if (rij==zero) cycle
!DEBUG
!    write(std_out,*) 'rij=',rij,' inx=',inx,' iny=',iny, ch10
!ENDDEBUG
! This corresponds to beta*sigma_ij in the original paper:
                   fij=beta*sqrt(wannspr(iwan,isppol)*wannspr(iwan,isppol)+wannspr(jwan,isppol)*wannspr(jwan,isppol))
                   erfValue = abi_derf(rij/fij)

                   if (ii==jj) then
                     ll=1
                   else
                     ll=0
                   end if ! ii==jj

                   Tij_dip(inx,iny)=-((3*rij_c(ii)*rij_c(jj)-rij*rij*ll)/rij**5)* &
&                   (erfValue-two*rij*exp(-((rij/fij)**2))/(sqrt(pi)*fij)) + &
&                   2.*two*rij_c(ii)*rij_c(jj)*exp(-((rij/fij)**2))/(sqrt(pi)*fij*fij*fij*rij*rij)

                   c_QHO(inx,iny)=omega(iwan)*omega(jwan)*sqrt(polar(iwan)*polar(jwan))*Tij_dip(inx,iny)

                 end do !jj=1,3
               end do  !jwan=1,nwan
             end do   !ii=1,3
           end do    !iwan=1,nwan
         end do     !icx=-abs(vdw_supercell(1)),ll
       end do      !icy=-abs(vdw_supercell(2)),mm
     end do       !icz=-abs(vdw_supercell(3)),0


   !writing terms for interactions from the central unit cell to each cell
   ! icx, icy, icz labels the cells in the supercell defined by vdw_supercell
   ! while inx and iny label QHO matrix elements.
   ! inx should start from (nc/2)*3*mwan + 1 --> r=-r_central+r_cell, and r_cell=displaced positions
   ! iny starts from  (nc/2)*3*mwan+3*man+1   to avoid central cell intra interactions which were built
   ! before

     iny=(nc/2)*3*mwan+3*mwan

     do icz=0,abs(vdw_supercell(3))
       if (icz==0) then
         mm=1
       else
         mm=-abs(vdw_supercell(2))
       end if
       do icy=mm,abs(vdw_supercell(2))
         if (icy==0) then
           ll=1
         else
           ll=-abs(vdw_supercell(1))
         end if
         do icx=ll,abs(vdw_supercell(1))
           do iwan=1,nwan(isppol)
             do ii=1,3
               iny=iny+1
       ! loop over the MLWFs in the 'central' unit cell:
               inx=(nc/2)*3*mwan
               do jwan=1,nwan(isppol)
                 do jj=1,3
                   inx=inx+1

                   if (inx==iny) cycle !in order to avoid digonal terms which  were already computed


                   dcenters(:,iwan,isppol) = (real(icx,dp))*rprimd(:,1)+(real(icy,dp))*rprimd(:,2)+&
&                   (real(icz,dp))*rprimd(:,3)+wanncent(:,iwan,isppol)

                   rij_c = dcenters(:,iwan,isppol)-wanncent(:,jwan,isppol)
                   rij = dnrm2(3,rij_c,1)
!                  rij=sqrt(dot_product(rij_c,rij_c))
                   if(rij==zero) cycle
!DEBUG
!    write(std_out,*) 'rij=',rij,' inx=',inx,' iny=',iny, ch10
!ENDDEBUG
! This corresponds to beta*sigma_ij in the original paper:
                   fij=beta*sqrt(wannspr(iwan,isppol)*wannspr(iwan,isppol)+wannspr(jwan,isppol)*wannspr(jwan,isppol))
                   erfValue = abi_derf(rij/fij)

                   if (ii==jj) then
                     ll=1
                   else
                     ll=0
                   end if ! ii==jj

                   Tij_dip(inx,iny)=-((3*rij_c(ii)*rij_c(jj)-rij*rij*ll)/rij**5)* &
&                   (erfValue-two*rij*exp(-((rij/fij)**2))/(sqrt(pi)*fij)) + &
&                   2.*two*rij_c(ii)*rij_c(jj)*exp(-((rij/fij)**2))/(sqrt(pi)*fij*fij*fij*rij*rij)

                   c_QHO(inx,iny)=omega(iwan)*omega(jwan)*sqrt(polar(iwan)*polar(jwan))*Tij_dip(inx,iny)

                 end do !jj=1,3
               end do  !jwan=1,nwan
             end do   !ii=1,3
           end do    !iwan=1,nwan
         end do     !icx=-abs(vdw_supercell(1)),ll
       end do      !icy=-abs(vdw_supercell(2)),mm
     end do       !icz=-abs(vdw_supercell(3)),0


! Here we diagonalize the matrix c_QHO and the eigenvalues come back in vector eigv
     ABI_MALLOC(matrx,((3*mwan*nc*(3*mwan*nc+1))/2))
     ABI_MALLOC(eigv,(3*mwan*nc))
     ABI_MALLOC(eigvec,(3*mwan*nc,3*mwan*nc))
     ABI_MALLOC(zhpev1,(3*2*mwan*nc-1))
     ABI_MALLOC(zhpev2,(3*3*mwan*nc-2))
     matrx(:)=cmplx(zero,zero)
     do jj=1,3*mwan*nc
       do ii=1,jj
         matrx(ii+(jj-1)*jj/2)=cmplx(c_QHO(ii,jj),0.0d0)
       end do
     end do

!DEBUG
!   write(std_out,*) 'Printing the real part of elements in array matrx:',ch10
!   do jj=1,3*mwan*nc*(3*mwan*nc+1)/2
!    write(std_out,*) real(matrx(jj))
!   enddo
!ENDDEBUG
     call ZHPEV ('N','U',3*mwan*nc,matrx,eigv,eigvec,3*mwan*nc,zhpev1,zhpev2,ier)
!DEBUG
     write(std_out,*) 'Last argument of ZHPEV: ier=',ch10
     write(std_out,*) ier,ch10
     write(std_out,*) 'List of c_QHO eigenvaules:',ch10
     do ll=1,3*mwan*nc
       write(std_out,*) eigv(ll)
     end do
!ENDDEBUG
     if(ier/=0) then !vz_d
       ABI_ERROR('zhpev fails!') !vz_d
     end if !vz_d

     ABI_FREE(matrx)
     ABI_FREE(eigvec)
     ABI_FREE(zhpev1)
     ABI_FREE(zhpev2)

     do ii=1,3*mwan*nc  !3*nwan(isppol)
       corrvdw=corrvdw+sqrt(eigv(ii))
     end do

   end do  ! end isppol


   corrvdw=0.5*corrvdw

!DEBUG
   write(std_out,*) 'Half the sum of interacting matrix eigenvalues square roots:',ch10
   write(std_out,*) corrvdw,ch10
!ENDDEBUG



   corrvdw=corrvdw-1.5d0*fu

   ABI_FREE(c_QHO)
   ABI_FREE(Tij_dip)
   ABI_FREE(dcenters)
   ABI_FREE(eigv)
   ABI_FREE(polar)
   ABI_FREE(omega)

   write(message, '(2a,i2,2a,f12.6,2a,f12.6,a)' )ch10,&
&   ' vdw_xc : ',14,ch10,&
&   ' van der Waals correction(Ha):',   corrvdw,ch10,&
&   ' van der Waals correction(eV):',   corrvdw*Ha_ev,ch10
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

 end if ! vdw-QHO

 if(allocated(ord))then
   ABI_FREE(ord)
 end if
 if(allocated(wanncent))then
   ABI_FREE(wanncent)
 end if
 if(allocated(wannspr))then
   ABI_FREE(wannspr)
 end if
 if(allocated(newocc_wan))then
   ABI_FREE(newocc_wan)
 end if
 if(allocated(npwf))then
   ABI_FREE(npwf)
 end if
 if (allocated(inwan))then
   ABI_FREE(inwan)
 end if
 if(allocated(nw))then
   ABI_FREE(nw)
 end if
 if(allocated(nwan))then
   ABI_FREE(nwan)
 end if
 ABI_FREE(wc_rec)
 if(allocated(amagr))then
   ABI_FREE(amagr)
 end if
 if(allocated(amawf))then
   ABI_FREE(amawf)
 end if
 if(allocated(Tij_dip))then
   ABI_FREE(Tij_dip)
 end if
 if(allocated(c_QHO))then
   ABI_FREE(c_QHO)
 end if
 if(allocated(amaspr))then
   ABI_FREE(amaspr)
 end if
 if(allocated(amaocc))then
   ABI_FREE(amaocc)
 end if
 if(allocated(tmp_cent))then
   ABI_FREE(tmp_cent)
 end if
 if(allocated(tmp_nwan))then
   ABI_FREE(tmp_nwan)
 end if
 if(allocated(tmp_spr))then
   ABI_FREE(tmp_spr)
 end if
 if(allocated(tmp_occ))then
   ABI_FREE(tmp_occ)
 end if


end subroutine evdw_wannier
!!***

!!****f* ABINIT/getFu
!! NAME
!! getFu
!!
!! FUNCTION
!!  Performs double integral needed to evaluate C6
!!  coefficients. Eq. (9) in J.Phys.Chem. 113:5224 [[cite:Silvestrelli2009]]
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_evdw_wannier
!!
!! CHILDREN
!!
!! SOURCE

 subroutine getFu(sn,sl,rn,rl,occn,occl,fu) ! sn-->spread(n), sl-->spread(l), rn --> rc(n), rl --> rc(l)

 implicit none
 real(dp),intent(in)::sn,sl,rn,rl,occn,occl
 real(dp),intent(out)::fu
 !local variables
 integer::nx,ny,ix,iy
 real(dp)::deltax,deltay
 real(dp)::beta,xc,yc,y,x
 real(dp),allocatable::arg1(:),res1(:),arg2(:),res2(:)

! *************************************************************************

 ny=100
 nx=100
 beta=(sn/sl)**(1.5d0)
 xc=rn
 yc=rl
 deltax=xc/(real(nx,dp)-1.d0)
 deltay=yc/(real(ny,dp)-1.d0)

 ABI_MALLOC(arg1,(ny))
 ABI_MALLOC(res1,(ny))
 ABI_MALLOC(arg2,(nx))
 ABI_MALLOC(res2,(nx))

 do ix=1,nx

   x=deltax*(real(ix,dp)-1.d0)

   do iy=1,ny
     y=deltay*(real(iy,dp)-1.d0)
     arg1(iy)=( (y**2.d0)*exp(-y) )/( (exp(-x)/(beta*sqrt(occn))) + exp(-y)/(sqrt(occl)) )
   end do

   call simpson_int(ny,deltay,arg1,res1)
   arg2(ix)=(x**2.d0)*exp(-x)*res1(ny)

 end do

 call simpson_int(nx,deltax,arg2,res2)

 Fu = res2(nx)

 ABI_FREE(arg1)
 ABI_FREE(res1)
 ABI_FREE(arg2)
 ABI_FREE(res2)
end subroutine getFu
!!***

!!****f* ABINIT/order_wannier
!! NAME
!! order_wannier
!!
!! FUNCTION
!!  Assign each MLWF with a corresponding fragment of atoms, according
!!  to vdw_typfrag array. Assignation is done by evaluating the distance
!!  from each MLWF center to the unit cell atoms. MLWFs belong to the
!!  same fragment as their nearest atom.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_evdw_wannier
!!
!! CHILDREN
!!
!! SOURCE
 subroutine order_wannier(mwan,natom,nwan,nsppol,ord,vdw_typfrag,wanncent,xcart)

   implicit none
!Arguments
   integer, intent(in)    :: mwan,natom,nsppol,nwan(nsppol),vdw_typfrag(natom) !vz_d
   integer, intent(inout) :: ord(mwan,nsppol)
   real(dp),intent(in)    :: wanncent(3,mwan,nsppol),xcart(3,natom)
!Local variables
   integer :: ii,jj,ll
   real(dp):: dis,dnrm2,mindi
   real(dp), allocatable :: tmp(:)
! *************************************************************************

 ABI_MALLOC(tmp,(3))

 do ll=1,nsppol
   do ii=1,nwan(ll)
     tmp(:) = wanncent(:,ii,ll) - xcart(:,1)
     mindi = dnrm2(3,tmp,1)
!     mindi=sqrt( dot_product(wanncent(:,ii,ll),wanncent(:,ii,ll))+dot_product(xcart(:,1),xcart(:,1))&
!&     -2*(dot_product(wanncent(:,ii,ll),xcart(:,1))) )
     ord(ii,ll)=vdw_typfrag(1)
     do jj=2,natom
       tmp(:) = wanncent(:,ii,ll) - xcart(:,jj)
       dis = dnrm2(3,tmp,1)
!       dis=sqrt( dot_product(wanncent(:,ii,ll),wanncent(:,ii,ll))+dot_product(xcart(:,jj),xcart(:,jj))&
!&       -2*(dot_product(wanncent(:,ii,ll),xcart(:,jj))) )
       if(dis<=mindi) then
         mindi=dis
         ord(ii,ll)=vdw_typfrag(jj)
       end if
     end do
   end do
 end do

 ABI_FREE(tmp)

 end subroutine order_wannier
!!***

!!****f* ABINIT/ovlp_wann
!! NAME
!! ovlp_wann
!!
!! FUNCTION
!!  Evaluate volumen reduction of MLWFs
!!  due to intrafragment overlapping
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_evdw_wannier
!!
!! CHILDREN
!!
!! SOURCE
 subroutine ovlp_wann(mwan,nwan,nsppol,ord,wanncent,wannspr,xi)

   implicit none
!Arguments
   integer, intent(in)  :: mwan,nsppol,nwan(nsppol),ord(mwan,nsppol) !vz_d
   real(dp),intent(in)  :: wanncent(3,mwan,nsppol),wannspr(mwan,nsppol)
   real(dp), intent(out)  :: xi(mwan,nsppol)
!Local variables
   integer :: ii,iwan,ix,iy,iz,jj,jwan,neigh,steps
   real(dp):: dis,disi,discent,veff,vfree,dnrm2
   integer, allocatable :: intsec(:,:,:,:)
   real(dp), allocatable :: rpoint(:), tmp(:)
   real(dp), parameter :: delt = 0.05d0 !Bohr, spatial mesh (1D) step
! *************************************************************************

 ABI_MALLOC(intsec,(mwan,nsppol,mwan,nsppol))
 ABI_MALLOC(rpoint,(3))
 ABI_MALLOC(tmp,(3))
 intsec(:,:,:,:) = 0
 xi(:,:) = 0.0d0
!detecting WF intersecting neighbors
 do ii=1,nsppol
   do iwan=1,nwan(ii)
     do jj=1,nsppol
       do jwan=1,nwan(jj)
         dis = 0.0d0
         if (ord(jwan,jj)==ord(iwan,ii)) then


           tmp(:) = wanncent(:,iwan,ii) - wanncent(:,jwan,jj)
           dis =  dnrm2(3,tmp,1)
!           dis=sqrt(  dot_product(wanncent(:,iwan,ii),wanncent(:,iwan,ii))+&
!&           dot_product(wanncent(:,jwan,jj),wanncent(:,jwan,jj))&
!&           -2*( dot_product(wanncent(:,iwan,ii),wanncent(:,jwan,jj)) )  )

           if ( ii == jj ) then
             if ( dis<=(wannspr(iwan,ii)+wannspr(jwan,jj)).and.iwan/=jwan ) then
               intsec(iwan,ii,jwan,jj) = 1
             end if
           end if
           if ( ii /= jj) then
             if ( dis<=(wannspr(iwan,ii)+wannspr(jwan,jj)) ) then
               intsec(iwan,ii,jwan,jj) = 1
             end if
           end if

         end if
       end do
     end do
   end do
 end do

!DEBUG
 write(std_out,*) 'intsec(iwan,ii,jwan,jj)=',ch10
 do ii=1,nsppol
   do iwan=1,nwan(ii)
     do jj=1,nsppol
       write(std_out,*) (intsec(iwan,ii,jwan,jj),jwan=1,nwan(jj)),ch10
     end do
   end do
 end do
!END DEBUG
!Determining both free and effective volumes.
!Eqs (6) and (7) in PRB 85:073101 [[cite:Ambrosetti2012]].
!Creation of grids around each WF centre.
!Calculation of intersection volumes.
 do ii = 1,nsppol
   do iwan = 1,nwan(ii)
!    Spatial meshes and volume parameters
     steps=NINT(wannspr(iwan,ii)/delt)
     vfree = 0
     veff = 0
     rpoint(:) = 0.0d0
     do iz=-steps,steps
       do iy=-steps,steps
         do ix=-steps,steps
           neigh = 0
           rpoint(1) = wanncent(1,iwan,ii) + ix*delt
           rpoint(2) = wanncent(2,iwan,ii) + iy*delt
           rpoint(3) = wanncent(3,iwan,ii) + iz*delt

           tmp(:) = wanncent(:,iwan,ii) - rpoint(:)
           discent = dnrm2(3,tmp,1)

!           discent = sqrt( dot_product(wanncent(:,iwan,ii),wanncent(:,iwan,ii))&
!&           +dot_product( rpoint(:),rpoint(:) )&
!&           -2*( dot_product( wanncent(:,iwan,ii),rpoint(:) ) ) )

           if (discent > wannspr(iwan,ii)) cycle
           if (discent <= wannspr(iwan,ii)) then

             neigh = 1
             do jj = 1,nsppol
               do jwan = 1,nwan(jj)
                 if ( intsec(iwan,ii,jwan,jj) == 0 ) cycle
                 if ( intsec(iwan,ii,jwan,jj) == 1 ) then

                   tmp(:) = rpoint(:) - wanncent(:,jwan,jj)
                   disi = dnrm2(3,tmp,1)
!                   disi = sqrt( dot_product(rpoint(:),rpoint(:))&
!&                   +dot_product( wanncent(:,jwan,jj),wanncent(:,jwan,jj) )&
!&                   -2*( dot_product(rpoint(:),wanncent(:,jwan,jj)) ) )
                   if (disi <= wannspr(jwan,jj)) then
                     neigh = neigh + 1
                   end if
                 end if
               end do
             end do
             if (nsppol==1) then
               veff = veff + 1/(real(neigh,dp)**2)
             end if
             if (nsppol==2) then
               veff = veff + 1/real(neigh,dp)
             end if
             vfree = vfree + 1/real(neigh,dp)
           end if
         end do
       end do
     end do
!    write(std_out,*) 'iwan=',iwan,'ii=',ii,ch10
!    write(std_out,*) 'vfree=',vfree,'neigh=',neigh,'veff=',veff,ch10
     xi(iwan,ii) = veff/vfree
!    write(std_out,*) 'xi(iwan,ii)=',xi(iwan,ii),ch10
   end do
 end do

 ABI_FREE(intsec)
 ABI_FREE(rpoint)
 ABI_FREE(tmp)

 end subroutine ovlp_wann
!!***

!!****f* ABINIT/vv10limit
!! NAME
!! vv10limit
!!
!! FUNCTION
!!  Performs double integral needed to evaluate C6
!!  coefficients from the long range limit of VV10
!!  functional (Phys. Rev. A. 81:062708 (2010)) [[cite:Vydrov2010]]
!!  as expressed in terms of MLWFs.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_evdw_wannier
!!
!! CHILDREN
!!
!! SOURCE

 subroutine vv10limit(sn,sl,rn,rl,fu) ! sn-->spread(n), sl-->spread(l), rn --> rc(n), rl --> rc(l)

 implicit none
 real(dp),intent(in)::sn,sl,rn,rl
 real(dp),intent(out)::fu
 !local variables
 integer::nx,ny,ix,iy
 real(dp)::deltax,deltay,pown,powl
 real(dp)::xc,yc,y,x,wgn,wgl,wox,woy
 real(dp),parameter :: cons = 0.0093d0 !related to local band gap model, VV10
 real(dp),allocatable::arg1(:),res1(:),arg2(:),res2(:)
! *************************************************************************

 ny=1000
 nx=1000

 xc=rn
 yc=rl
 deltax=xc/(real(nx,dp)-1.d0)
 deltay=yc/(real(ny,dp)-1.d0)

 ABI_MALLOC(arg1,(ny))
 ABI_MALLOC(res1,(ny))
 ABI_MALLOC(arg2,(nx))
 ABI_MALLOC(res2,(nx))

 wgn = cons*( (18.0d0/(sn*sqrt3**three))**4 )
 pown = two*sqrt3/sn
 wgl = cons*( (18.0d0/(sl*sqrt3**three))**4 )
 powl = two*sqrt3/sl

 do ix=1,nx

   x = deltax*(real(ix,dp)-1.d0)
   wox = sqrt(wgn + (four*pown/sn**two)*exp(-pown*x))

   do iy=1,ny

     y = deltay*(real(iy,dp)-1.d0)
     woy = sqrt(wgl + (four*powl/sl**two)*exp(-powl*y))

     arg1(iy)=( (y**two)*exp(-powl*y) )/( woy*(wox+woy) )

   end do

   call simpson_int(ny,deltay,arg1,res1)
   arg2(ix)=(x**two)*exp(-pown*x)*res1(ny)/wox

 end do

 call simpson_int(nx,deltax,arg2,res2)

 fu = res2(nx)

!DEBUG
 write(std_out,*) ch10,'Int argument',ch10
 do ix=1,nx
   write(std_out,*) deltax*(real(ix,dp)-1.d0), arg2(ix)
 end do
!END DEBUG

 ABI_FREE(arg1)
 ABI_FREE(res1)
 ABI_FREE(arg2)
 ABI_FREE(res2)
end subroutine vv10limit
!!***

!!****f* ABINIT/amalgam
!! NAME
!! amalgam
!!
!! FUNCTION
!!  Amalgamates MLWFs, which are close enough,
!!  into one MLWF as suggested in J.Chem.Phys.135:154105 (2011) [[cite:Andrinopoulos2011]]
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_evdw_wannier
!!
!! CHILDREN
!!
!! SOURCE

 subroutine amalgam(amagr,ngr,nsppol,nw,mwan,ord,nwan,vdw_nfrag,wanncent,wannspr)

 implicit none
 !Arguments
 integer,intent(in) :: nsppol,mwan,vdw_nfrag
 integer,intent(in) :: ord(mwan,nsppol),nwan(nsppol)
 real(dp),intent(in):: wanncent(3,mwan,nsppol),wannspr(mwan,nsppol)
 integer,intent(out):: ngr
 integer,intent(out):: nw(nsppol,mwan/2),amagr(mwan,nsppol,mwan/2)
 !local variables
 integer :: dimen,ii,igr,isppol,iw,iwan,jj,jsppol,jwan,ll
 real(dp):: dis, dnrm2
 real(dp),allocatable :: tmp(:)
! *************************************************************************

 ABI_MALLOC(tmp,(3))

!Selecting pairs of MLWFs satisfying amalgamation criteria
 write(std_out,*) 'Searching for MLWFs close enough to amalgamate...',ch10

 dimen = iabs(vdw_nfrag)

!Grouping MLWFs and amalgamation

 ngr = 0
 amagr(:,:,:) = 0
 nw(:,:) = 0

 do ll = 1 , dimen
   do isppol = 1 , nsppol
     jsppol = isppol
     do iwan = 2 , nwan(isppol)
       do jwan = 1 , iwan-1

         if (ord(iwan,isppol)==ll .and. ord(jwan,jsppol)==ll ) then

           tmp(:) = wanncent(:,iwan,isppol) - wanncent(:,jwan,jsppol)
           dis = dnrm2(3,tmp,1)

!           dis=sqrt( dot_product(wanncent(:,iwan,isppol),wanncent(:,iwan,isppol)) &
!&           + dot_product(wanncent(:,jwan,jsppol),wanncent(:,jwan,jsppol))&
!&           - 2*(dot_product(wanncent(:,iwan,isppol),wanncent(:,jwan,jsppol))) )

           if ( dis <= (wannspr(iwan,isppol) + wannspr(jwan,jsppol)) / three ) then

             if ( all(amagr(:,isppol,:) /= iwan) .and. &
&             all(amagr(:,jsppol,:) /= jwan) ) then

               ngr = ngr + 1
               amagr(1,isppol,ngr) = jwan
               amagr(2,jsppol,ngr) = iwan
               nw(isppol,ngr) = 2
               cycle

             end if

             if  ( any(amagr(:,isppol,:) == iwan) .and. &
&             any(amagr(:,jsppol,:) == jwan) ) cycle

             do igr = 1 , mwan/2
               do iw = 1 , mwan

                 if ( amagr(iw,isppol,igr) ==  jwan .and. &
&                 all(amagr(:,isppol,igr) /= iwan) ) then
                   nw(isppol,igr) = nw(isppol,igr) + 1
                   amagr(nw(isppol,igr),isppol,igr) = iwan
                   cycle
                 end if

                 if ( amagr(iw,isppol,igr) ==  iwan .and. &
&                 all(amagr(:,isppol,igr) /= jwan) ) then
                   nw(isppol,igr) = nw(isppol,igr) + 1
                   amagr(nw(isppol,igr),isppol,igr) = jwan
                   cycle
                 end if

               end do
             end do

           end if  !if dis < (wannspr(iwan,isppol) + wannspr(jwan,jsppol))/three
         end if  !if (ord(iwan,isppol)==ll .and. ord(jwan,jsppol)==ll )
       end do  !jwan
     end do  !iwan
   end do  !isppol


   if (nsppol == 2) then
     isppol = 1
     jsppol = 2
     do iwan = 1 , nwan(isppol)
       do jwan = 1 , nwan(jsppol)

         if (ord(iwan,isppol)==ll .and. ord(jwan,jsppol)==ll ) then

           tmp(:) =  wanncent(:,iwan,isppol) - wanncent(:,jwan,jsppol)
           dis = dnrm2(3,tmp,1)

!           dis=sqrt( dot_product(wanncent(:,iwan,isppol),wanncent(:,iwan,isppol)) &
!&           + dot_product(wanncent(:,jwan,jsppol),wanncent(:,jwan,jsppol))&
!&           - 2*(dot_product(wanncent(:,iwan,isppol),wanncent(:,jwan,jsppol))) )

           if ( dis <= (wannspr(iwan,isppol) + wannspr(jwan,jsppol)) / three ) then

             if ( all(amagr(:,isppol,:) /= iwan) .and. &
&             all(amagr(:,jsppol,:) /= jwan) ) then

               ngr = ngr + 1
               amagr(1,isppol,ngr) = iwan
               amagr(1,jsppol,ngr) = jwan
               nw(isppol,ngr) = nw(isppol,ngr) + 1
               nw(jsppol,ngr) = nw(jsppol,ngr) + 1
               cycle

             end if

             if  ( any(amagr(:,isppol,:) == iwan) .and. &
&             any(amagr(:,jsppol,:) == jwan) ) cycle

             do igr = 1 , mwan/2
               do iw = 1 , mwan

                 if ( amagr(iw,jsppol,igr) ==  jwan .and. &
&                 all(amagr(:,isppol,igr) /= iwan) ) then
                   nw(isppol,igr) = nw(isppol,igr) + 1
                   amagr(nw(isppol,igr),isppol,igr) = iwan
                   cycle
                 end if

                 if ( amagr(iw,isppol,igr) ==  iwan .and. &
&                 all(amagr(:,jsppol,igr) /= jwan) ) then
                   nw(jsppol,igr) = nw(jsppol,igr) + 1
                   amagr(nw(jsppol,igr),jsppol,igr) = jwan
                   cycle
                 end if

               end do
             end do

           end if

         end if

       end do
     end do
   end if !if (nsppol == 2)

 end do !ll

 write(std_out,*) 'Number of amalgamation groups:',ngr,ch10
 if(ngr/=0)then
   do ii = 1 , ngr
     do isppol = 1 ,nsppol
       write(std_out,*) 'Number of MLWFs in group',ii,':',nw(isppol,ii),ch10
       write(std_out,*) 'MLWFs in group',ii,': WFindex,spin,group ',ch10
       do jj = 1, nw(isppol,ii)
         write(std_out,*) amagr(jj,isppol,ii),isppol,ii,ch10
       end do
     end do
   end do
 end if

!DEBUG
!write(std_out,*)' amalgam : exit '
!write(std_out,*)' nw =',nw
!call flush
!ENDDEBUG

 ABI_FREE(tmp)
 end subroutine amalgam
!!***

end module m_evdw_wannier
!!***
