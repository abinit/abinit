!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_qparticles
!! NAME
!!  m_qparticles
!!
!! FUNCTION
!!  This module contains tools for the IO of the QP file and other procedures
!!  related to the calculation of the quasiparticle amplitudes represented in terms
!!  of KS states.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2019 ABINIT group (FB, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_qparticles

 use defs_basis
 use m_abicore
 use m_hdr
 use m_errors
 use m_nctk
 use m_distribfft


 use defs_datatypes,   only : pseudopotential_type, ebands_t
 use defs_abitypes,    only : MPI_type
 use m_io_tools,       only : open_file, file_exists, isncfile
 use m_fstrings,       only : int2char10, itoa, sjoin
 use m_numeric_tools,  only : linfit, c2r, set2unit, interpol3d, rhophi
 use m_gwdefs,         only : sigparams_t
 use m_crystal,        only : crystal_t
 use m_bz_mesh,        only : kmesh_t
 use m_ebands,         only : get_valence_idx
 use m_sigma,          only : sigma_t
 use m_pawtab,         only : pawtab_type
 use m_pawrhoij,       only : pawrhoij_type, pawrhoij_alloc, pawrhoij_io, pawrhoij_inquire_dim
 use m_fourier_interpol,only : fourier_interpol

 implicit none

 private

 public :: wrqps             ! Write a QPS file.
 public :: rdqps             ! Read a QPS file.
 public :: show_QP           ! Report the components of a QP amplitude in terms of KS eigenstates.
 public :: rdgw              ! Read GW corrections from an external file.
 public :: updt_m_lda_to_qp  ! Updates the matrix of unitary transformation from lda to qp states.

CONTAINS  !=======================================================================================
!!***

!!****f* m_qparticles/wrqps
!! NAME
!! wrqps
!!
!! FUNCTION
!!  Write the _QPS file containing information on the quasi-particles energies and wavefunctions.
!!
!! INPUTS
!!  fname=The name of the file
!!  Sigp<sigparams_t>=Parameters characterizing the self-energy calculation.
!!     %nsppol=1 for unpolarized, 2 for spin-polarized
!!     %nbnds=number of bands used for sigma
!!  Sr<sigma_t>=Structure containing the results of the sigma run.
!!     %en_qp_diago(nbnds,nibz,nsppol)= NEW quasi-particle energies
!!     %eigvec_qp(nbnds,nbnds,nibz,nsppol)= NEW QP amplitudes in the KS basis set
!!      obtained by diagonalizing H0 + Herm(Sigma).
!!  m_lda_to_qp(nbnds,nbnds,nibz,nsppol)= expansion of the OLD QP amplitudes in terms of KS wavefunctions
!!  Kmesh<kmesh_t>=information on the k-point sampling.
!!     %nibz=number of irreducible k-points
!!     %ibz(3,kibz)=reduced coordinates of the irreducible k-points
!!  nfftot=Total number of FFT points for density
!!  ngfftf(18)=Info on the FFT mesh for the density.
!!  nscf=Number of self consistent cycles performed
!!  nspden=number of spin-density components
!!  Cryst<crystal_t>=Structure defining the crystal structure.
!!  Psps<type(pseudopotential_type)>=variables related to pseudopotentials.
!!  Pawrhoij(Cryst%natom*Psps%usepaw)<type(pawrhoij_type)>= rhoij datastructure.
!!  BSt<ebands_t>=Structure containing the band structure energies (only used is nscf==-1)
!!
!! OUTPUT
!!  Only writing
!!
!! NOTES
!!  Old QPS fileformat:
!!   |
!!   | No. of QPSCF cycles already performed.
!!   | No. of k-points in the IBZ.
!!   | Total number of bands used to construct the Green's function (nbnds)
!!   | nsppol
!!   | For each spin and k-point in the IBZ:
!!   |   Reduced coordinates of the k-point.
!!   |   for each band:
!!   |     QP energies obtained by diagonalizing the QPSCGW Hamiltonian.
!!   |     <\psi_{ib,k,s}^{KS}|\psi_{jb,k,s}^{QP}>$, ib=1,nbnds
!!   | FFT dimensions of the fine grid
!!   | QP density in real space.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine wrqps(fname,Sigp,Cryst,Kmesh,Psps,Pawtab,Pawrhoij,nspden,nscf,nfftot,ngfftf,Sr,Bst,m_lda_to_qp,rho_qp)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot,nscf,nspden
 character(len=*),intent(in) :: fname
 type(kmesh_t),intent(in) :: Kmesh
 type(ebands_t),intent(in) :: BSt
 type(sigparams_t),intent(in) :: Sigp
 type(sigma_t),intent(in) :: Sr
 type(crystal_t),intent(in) :: Cryst
 type(Pseudopotential_type),intent(in) :: Psps
!arrays
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(in) :: rho_qp(nfftot,nspden)
 complex(dpc),intent(in) :: m_lda_to_qp(Sigp%nbnds,Sigp%nbnds,Kmesh%nibz,Sigp%nsppol)
 type(Pawrhoij_type),intent(inout) :: Pawrhoij(Cryst%natom*Psps%usepaw)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: ib,ik,is,unqps,iatom,itypat
 character(len=500) :: msg
!arrays
 integer,allocatable :: nlmn_type(:)
 complex(dpc),allocatable :: mtmp(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 if (nscf >= 0) then
   write(msg,'(3a)')ch10,' writing QP data on file : ',TRIM(fname)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
 end if

 if (open_file(fname,msg,newunit=unqps,form='formatted',status='unknown') /= 0) then
   MSG_ERROR(msg)
 end if

 write(unqps,*)nscf+1
 write(unqps,*)Kmesh%nibz
 write(unqps,*)Sigp%nbnds
 write(unqps,*)Sigp%nsppol

 ABI_MALLOC(mtmp,(Sigp%nbnds,Sigp%nbnds))

 if (nscf>=0) then
   ! Write the new m_lda_to_qp on file.
   do is=1,Sigp%nsppol
     do ik=1,Kmesh%nibz
       write(unqps,*)Kmesh%ibz(:,ik)
       do ib=1,Sigp%nbnds
         write(unqps,*)Sr%en_qp_diago(ib,ik,is)
         write(unqps,*)m_lda_to_qp(:,ib,ik,is)
       end do
     end do
   end do
 else if (nscf==-1) then
   ! Write fake QPS file with KS band structure (Mainly used for G0W)
   call set2unit(mtmp)
   do is=1,Sigp%nsppol
     do ik=1,Kmesh%nibz
       write(unqps,*)Kmesh%ibz(:,ik)
       do ib=1,Sigp%nbnds
         write(unqps,*)BSt%eig(ib,ik,is)
         write(unqps,*)mtmp(:,ib)
       end do
     end do
   end do
 else
   MSG_ERROR(sjoin("Wrong nscf ",itoa(nscf)))
 end if

 ABI_FREE(mtmp)

 write(msg,'(a,f9.4)')' (wrqps) planewave contribution to nelect: ',SUM(rho_qp(:,1))*Cryst%ucvol/nfftot
 call wrtout(std_out,msg,'COLL')
 if (nspden == 4) then
   write(msg,'(a,3f9.4)')' mx, my, mz: ',&
     SUM(rho_qp(:,2))*Cryst%ucvol/nfftot,SUM(rho_qp(:,3))*Cryst%ucvol/nfftot,SUM(rho_qp(:,4))*Cryst%ucvol/nfftot
   call wrtout(std_out,msg,'COLL')
 end if

 ! Write FFT dimensions and QP density
 write(unqps,*)ngfftf(1:3)
 write(unqps,*)rho_qp(:,:)

 if (Psps%usepaw==1) then
   ! Write QP rhoij to be used for on-site density mixing.
   ABI_MALLOC(nlmn_type,(Cryst%ntypat))
   do itypat=1,Cryst%ntypat
     nlmn_type(itypat)=Pawtab(itypat)%lmn_size
   end do

   write(unqps,*) Cryst%natom, Cryst%ntypat
   write(unqps,*) (Cryst%typat(iatom), iatom=1,Cryst%natom)
   write(unqps,*) (nlmn_type(itypat), itypat=1,Cryst%ntypat)
   write(unqps,*) Pawrhoij(1)%nsppol, Pawrhoij(1)%nspden

   call pawrhoij_io(pawrhoij,unqps,Sigp%nsppol,Sigp%nspinor,nspden,nlmn_type,Cryst%typat,&
&                HDR_LATEST_HEADFORM,"Write",form="formatted")
   ABI_FREE(nlmn_type)
 end if

 close(unqps)

 DBG_EXIT("COLL")

end subroutine wrqps
!!***

!----------------------------------------------------------------------

!!****f* m_qparticles/rdqps
!! NAME
!! rdqps
!!
!! FUNCTION
!!  Read a _QPS file containing the QP energies of the previous iteration, the coefficients
!!  defining the QP amplitudes in terms of the KS basis set and the QP density for mixing.
!!
!! INPUTS
!!  nfftot=Total number of FFT points for density
!!  ngfftf(18)=Info on the FFT mesh for the density.
!!  nspden=Number of SPin-DENsity components.
!!  usepaw=1 if we are using PAW.
!!  fname=Name of the file
!!  dimrho=1 if density has to be read, 0 otherwise
!!  BSt<ebands_t>=Structure containing the initial band structure.
!!  ucvol=Volume of the unit cell
!!
!! OUTPUT
!!  nbsc=number of bands used to describe the QP amplitudes
!!  nscf=number of iterations that have been performed (==0 if we start from a KS calculation)
!!  m_lda_to_qp(mband,mband,nibz,nsppol)=matrix giving the decomposition of the QP
!!   wavefunction in the mainfold generated by the KS wavefunctions
!!   (i.e. $ m_lda_to_qp(ib,jb,k,s) := <\psi_{ib,k,s}^{KS}|\psi_{jb,k,s}^{QP}>$
!!  rhor_out(nfftot,nspden)=quasiparticle density
!!
!! SIDE EFFECTS
!!  BSt<ebands_t>=Structure containing the initial band structure.
!!     %en_qp(mband,nkpt,nsppol)=QP energies at iteration nscf
!!
!! TODO
!!  The value of nspden is not reported in the QPS file thus we have a possible undetected error.
!!
!! PARENTS
!!      bethe_salpeter,mlwfovlp_qp,screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine rdqps(BSt,fname,usepaw,nspden,dimrho,nscf,&
& nfftot,ngfftf,ucvol,Cryst,Pawtab,MPI_enreg,nbsc,m_lda_to_qp,rhor_out,Pawrhoij)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftot,nspden,usepaw,dimrho
 integer,intent(out) :: nbsc,nscf
 real(dp),intent(in) :: ucvol
 character(len=*),intent(in) :: fname
 type(crystal_t),intent(in) :: Cryst
 type(ebands_t),intent(inout) :: BSt
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(out) :: rhor_out(nfftot,nspden*dimrho)
 complex(dpc),intent(out) :: m_lda_to_qp(BSt%mband,BSt%mband,BSt%nkpt,BSt%nsppol)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat*usepaw)
 type(Pawrhoij_type),intent(inout) :: Pawrhoij(Cryst%natom*usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: ib,ii,ik,isppol,nbandR,nkibzR,nsppolR,unqps,my_rank,ispden
 integer :: ifft,n1,n2,n3,ir1,ir2,ir3,ios
 integer :: cplex_fft,optin,optout,nfft_found
 integer :: iatom,natomR,nspdenR,ntypatR,itypat
 real(dp) :: uerr,nelect_qps,ratio
 logical,parameter :: use_FFT_interpolation=.TRUE.
 logical :: ltest
 character(len=500) :: msg
!arrays
 integer :: ngfft_found(18)
 integer,allocatable :: nlmn_type(:),typatR(:)
 real(dp) :: kibz(3),rr(3),rhogdum(1,1)
 real(dp),allocatable :: en_tmp(:)
 real(dp),allocatable :: rhor_tmp(:,:)
 complex(dpc),allocatable :: mtmp(:,:),utest(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 ABI_CHECK(ALL(BSt%nband==BSt%nband(1)),"No. of bands must be constant")
 ABI_CHECK(dimrho==0.or.dimrho==1,'dimrho must be 0 or 1')

 ! This does not work in parallel !!?
 !% my_rank = xmpi_comm_rank(MPI_enreg%spaceComm)
 my_rank = MPI_enreg%me_kpt

 ! Check whether file exists or not.
 write(msg,'(5a)')ch10,&
  ' rdqps: reading QP wavefunctions of the previous step ',ch10,&
  '        looking for file ',TRIM(fname)
 call wrtout([std_out, ab_out], msg)

 if (.not.file_exists(fname)) then
   write(msg,'(2a)')' file not found, 1st iteration initialized with KS eigenelements ',ch10
   call wrtout([std_out, ab_out], msg)
   nscf=0; RETURN
 end if

 if (.not.isncfile(fname)) then
   if (open_file(fname,msg,newunit=unqps,form='formatted',status='unknown') /= 0) then
     MSG_ERROR(msg)
   end if

   ! TODO the _QPS file should contain additional information
   read(unqps,*)nscf
   write(msg,'(a,i4,a)')' Number of iteration(s) already performed: ',nscf,ch10
   call wrtout([std_out, ab_out], msg)

   read(unqps,*)nkibzR
   if (nkibzR/=BSt%nkpt) then
     write(msg,'(2(a,i0))')'Wrong number of k-points; Expected: ',BSt%nkpt,', Found: ',nkibzR
     MSG_ERROR(msg)
   end if

   read(unqps,*)nbandR
   nbsc=MIN(nbandR,BSt%mband)

   if (nbsc/=BSt%mband) then
     write(msg,'(3a,i4,a,i4)')&
      'QPS file contains less bands than that used in the present calculation ',ch10,&
      'Required: ',BSt%mband,', Found: ',nbandR
     MSG_WARNING(msg)
   end if

   if (nbsc/=nbandR) then
     write(msg,'(3a,i4,a)')&
      'The QPS file contains more bands than that used in the present calculation ',ch10,&
      'only the first ',nbandR,' bands will be read'
     MSG_COMMENT(msg)
   end if

   ABI_MALLOC(mtmp,(nbandR,nbandR))
   ABI_MALLOC(en_tmp,(nbandR))
   read(unqps,*)nsppolR

   ABI_CHECK(nsppolR==BSt%nsppol,"QPS generated with different nsppol")

   ! Read energies and transformation for each k-point and spin.
   ! TODO: The format of the QPS file must be standardized !
   ! For example we might add the occupation numbers.
   do isppol=1,BSt%nsppol
     do ik=1,BSt%nkpt
       read(unqps,*)kibz(:)
       write(msg,'(a,i5,a,3(f6.3,1x),4x,a,i2)')' Reading ik ',ik,')  k = ',kibz(:),' is = ',isppol
       call wrtout(std_out,msg,'COLL')
       ltest=(ALL(ABS(kibz(:)-BSt%kptns(:,ik))<0.001))
       ABI_CHECK(ltest,'Wrong k-point read')
       do ib=1,nbandR
         read(unqps,*)en_tmp(ib)
         read(unqps,*)mtmp(:,ib)
       end do

       ! Store transformation and update energies.
       m_lda_to_qp(1:nbsc,1:nbsc,ik,isppol)=mtmp(1:nbsc,1:nbsc)
       BSt%eig(1:nbsc,ik,isppol)=en_tmp(1:nbsc)

       ! Chech if matrix is unitary.
       ABI_MALLOC(utest,(nbsc,nbsc))
       utest(:,:) = TRANSPOSE(mtmp(1:nbsc,1:nbsc)) !this is just for the buggy gfortran
       utest(:,:) = MATMUL(CONJG(utest),mtmp(1:nbsc,1:nbsc))
       do ii=1,nbsc
         utest(ii,ii)=utest(ii,ii)-one
       end do
       uerr=MAXVAL(ABS(utest))
       if (uerr>tol6) then
         write(msg,'(a,es16.8)')' KS -> QP matrix is not unitary, MAX error = ',uerr
         MSG_WARNING(msg)
       end if
       ABI_FREE(utest)
     end do !ik
   end do !isppol

   ABI_FREE(mtmp)
   ABI_FREE(en_tmp)

   ! Read the QP density.
   ! The two FFT grids might differ. In case perform an FFT interpolation to have rhor on the input mesh.
   if (dimrho==1) then
     read(unqps,*)n1,n2,n3

     if (all(ngfftf(1:3)== [n1, n2, n3]) ) then
       read(unqps,*)rhor_out(:,:)
     else
       write(msg,'(2a,a,5(i3,a),i3)')&
        'FFT meshes differ. Performing Fourier interpolation. ',ch10,&
        'Found: ',n1,' x',n2,' x',n3,'; Expected: ',ngfftf(1),' x',ngfftf(2),' x',ngfftf(3)
       MSG_COMMENT(msg)

       ABI_MALLOC(rhor_tmp,(n1*n2*n3,nspden))
       read(unqps,*)rhor_tmp(:,:)

       if (use_FFT_interpolation) then
         ngfft_found(1:3)=(/n1,n2,n3/)
         ngfft_found(4)=2*(ngfft_found(1)/2)+1 ! 4:18 are not used, anyway!
         ngfft_found(5)=2*(ngfft_found(2)/2)+1
         ngfft_found(6)=ngfft_found(3)
         ngfft_found(7:18)=ngfftf(7:18)
         nfft_found=PRODUCT(ngfft_found(1:3)) !no FFT para

         cplex_fft =1 ! Real quantities.
         optin     =0 ! Input is taken from rhor.
         optout    =0 ! Output is only in real space.
         call destroy_distribfft(MPI_enreg%distribfft)
         call init_distribfft(MPI_enreg%distribfft,'c',MPI_enreg%nproc_fft,ngfftf(2),ngfftf(3))
         call init_distribfft(MPI_enreg%distribfft,'f',MPI_enreg%nproc_fft,ngfft_found(2),ngfft_found(3))

         call fourier_interpol(cplex_fft,nspden,optin,optout,nfft_found,ngfft_found,nfftot,ngfftf,&
           MPI_enreg,rhor_tmp,rhor_out,rhogdum,rhogdum)

       else
         ! Linear interpolation.
         do ispden=1,nspden
           do ir3=0,ngfftf(3)-1
             rr(3)=DBLE(ir3)/n3
             do ir2=0,ngfftf(2)-1
               rr(2)=DBLE(ir2)/n2
               do ir1=0,ngfftf(1)-1
                 rr(1)=DBLE(ir1)/n1
                 ifft = 1 +ir1 +ir2*ngfftf(1) +ir3*ngfftf(1)*ngfftf(2)
                 rhor_out(ifft,ispden) = interpol3d(rr,n1,n2,n3,rhor_tmp(:,ispden))
               end do
             end do
           end do
         end do
       end if

       ABI_FREE(rhor_tmp)
     end if

     ! Test the normalization of the QPS density.
     ! There might be errors due to the interpolation or the truncation of the G basis set
     ! Density will be renormalized in the caller since for PAW we still have to add the onsite contribution.
     if (usepaw==0) then
       nelect_qps=SUM(rhor_out(:,1))*ucvol/nfftot; ratio=BSt%nelect/nelect_qps
       write(msg,'(3(a,f9.4))')&
         ' Number of electrons calculated using the QPS density = ',nelect_qps,' Expected = ',BSt%nelect,' ratio = ',ratio
       call wrtout(std_out, msg)
       !!rhor_out(:,:)=ratio*rhor_out(:,:)
     end if

     if (usepaw==1) then
       ! Write QP_rhoij for on-site density mixing.
       read(unqps,*,iostat=ios)natomR,ntypatR
       if (ios/=0) then
         msg="Old version of QPS file found. DO NOT USE rhoqpmix for this run."
         MSG_WARNING(msg)
         call wrtout(ab_out,msg,"COLL")
         ! Init dummy rhoij just to avoid problems in sigma when rhoij is freed.
         call pawrhoij_inquire_dim(nspden_rhoij=nspdenR, nspden=nspden)
         call pawrhoij_alloc(Pawrhoij,1,nspdenR,BSt%nspinor,BSt%nsppol,Cryst%typat,pawtab=Pawtab)
         close(unqps)
         RETURN
       end if

       ABI_CHECK(natomR ==Cryst%natom, "mismatch in natom")
       ABI_CHECK(ntypatR==Cryst%ntypat,"mismatch in ntypat")
       ABI_MALLOC(nlmn_type,(ntypatR))
       ABI_MALLOC(typatR,(ntypatR))

       read(unqps,*)(typatR(iatom), iatom=1,natomR)
       ABI_CHECK(ALL(Cryst%typat==typatR),"mismatch in typat")

       read(unqps,*)(nlmn_type(itypat), itypat=1,ntypatR)
       do itypat =1,Cryst%ntypat
         if (nlmn_type(itypat)/=Pawtab(itypat)%lmn_size) then
           MSG_ERROR("mismatch in nlmn_type, check QPS file")
         end if
       end do

       read(unqps,*) nsppolR,nspdenR
       ABI_CHECK(nsppolR==BSt%nsppol,"mismatch in nsppol")
       ABI_CHECK(nspdenR==nspden    ,"mismatch in nspden")

       call pawrhoij_io(pawrhoij,unqps,BSt%nsppol,BSt%nspinor,nspden,nlmn_type,Cryst%typat,&
                        HDR_LATEST_HEADFORM,"Read",form="formatted")
       !% call pawrhoij_io(pawrhoij,std_out,BSt%nsppol,BSt%nspinor,nspden,nlmn_type,Cryst%typat,HDR_LATEST_HEADFORM,"Echo")

       ABI_FREE(nlmn_type)
       ABI_FREE(typatR)
     end if ! usepaw

   end if !dimrho=1

   close(unqps)

 else
   MSG_ERROR("netdf format not implemented")
 end if

 DBG_EXIT("COLL")

end subroutine rdqps
!!***

!----------------------------------------------------------------------

!!****f* m_qparticles/show_QP
!! NAME
!! show_QP
!!
!! FUNCTION
!! Print in a nice format (?) the expansion coefficients of the quasiparticle
!! amplitudes in terms of the KS eigenvectors
!!
!! INPUTS
!!  Bst<ebands_t>=Description of the band structure.
!!    %nsppol=1 for unpolarized, 2 for spin-polarized.
!!    %mband=Max number of bands (in GW doesn"t depend on k an spin)
!!    %nkpt=number of irreducible k-points.
!!    %eig(mband,nkpt,nsppol)= QP energies for each k-point, band and spin.
!!  m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)=matrix giving the decomposition of the QP
!!   amplitued in the mainfold generated by the KS wavefunctions
!!   (i.e $ m_lda_to_qp(ib,jb,k,s) := \langle \psi_{ib,k,s}^{KS}| \psi_{jb,k,s}^{QP}\rangle $
!!  fromb,tob=initial and final band index for QP, only states in this range are printed
!!  prtvol=Verbosity level (not used)
!!  unit=Unit number of the output file
!! tolmat[Optional]=Only components whose coefficient has modulus larger than tolmat are shown (default is 0.01)
!!
!! OUTPUT
!!  Only printing
!!
!! NOTES
!!  Only master node should call this routine.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine show_QP(Bst,m_lda_to_qp,fromb,tob,unit,prtvol,tolmat,kmask)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: fromb,tob
 integer,optional,intent(in) :: prtvol,unit
 real(dp),optional,intent(in) :: tolmat
 type(ebands_t),intent(in) :: Bst
!arrays
 logical,optional,intent(in) :: kmask(Bst%nkpt)
 complex(dpc),intent(in) :: m_lda_to_qp(Bst%mband,Bst%mband,Bst%nkpt,Bst%nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: NBRA=5
 logical,parameter :: use_rhophi=.True.
 integer :: ib_start,ib_stop,my_prtvol,counter,ib_KS,ib_QP,ikibz,isp,nspace,my_unt,nband_k
 real(dp) :: my_tolmat,rho,phi
 character(len=10) :: bks,bqp,k_tag,spin_tag
 character(len=500) :: KS_row,KS_ket,tmpstr,QP_ket
!arrays
 real(dp) :: cx(2)

! *********************************************************************

 my_unt   =std_out  ; if (PRESENT(unit  )) my_unt   =unit
 my_prtvol=0        ; if (PRESENT(prtvol)) my_prtvol=prtvol
 ib_start =1        ; if (PRESENT(fromb )) ib_start =fromb
 ib_stop  =Bst%mband; if (PRESENT(tob   )) ib_stop  =tob
 my_tolmat=0.001    ; if (PRESENT(tolmat)) my_tolmat=ABS(tolmat)

 ! I suppose nband_k is constant thus the check is done here.
 if (ib_start<=0       ) ib_start=1
 if (ib_start>Bst%mband) ib_start=Bst%mband
 if (ib_stop<=0        ) ib_stop=1
 if (ib_stop>Bst%mband ) ib_stop=Bst%mband

 ! Have to follow rules 7.f.
 write(my_unt,'(/,a,/,a,/,a,f6.3,a,/,a)')&
   ' '//REPEAT('*',76),&
&  ' ***** QP amplitudes expressed as linear combination of KS eigenstates. *****',&
&  ' ***** Only KS components whose modulus is larger than ',my_tolmat,' are shown  ***** ',&
&  ' '//REPEAT('*',76)
 if (use_rhophi) then
   write(my_unt,"(a)")"Complex coefficients given in (rho, phi) polar representation."
 else
   write(my_unt,"(a)")"Complex coefficients given in (Re, Im) representation."
 end if

 if (PRESENT(kmask)) then
   if (.not.ALL(kmask)) write(my_unt,'(/,a,i3,a)')' Only ',COUNT(kmask),' k-points are reported '
 end if

 do isp=1,Bst%nsppol
   call int2char10(isp,spin_tag)
   write(my_unt,'(/,a,i2,a,/)')' >>>>> Begin block for spin ',isp,' <<<<< '

   do ikibz=1,Bst%nkpt
     if (PRESENT(kmask)) then
       if (.not.kmask(ikibz)) CYCLE
     end if
     call int2char10(ikibz,k_tag)
     nband_k=Bst%nband(ikibz+(isp-1)*Bst%nkpt)
     write(my_unt,'(a,i4,a,3es16.8,a,f6.3,/)')' k-point: ',ikibz,') ',Bst%kptns(:,ikibz),'; wtk= ',Bst%wtk(ikibz)

     do ib_QP=ib_start,ib_stop
       call int2char10(ib_QP,bqp)
       QP_ket=' |QP: b='//TRIM(bqp)//'; s='//TRIM(spin_tag)//'> = '
       write(my_unt,'(a)')TRIM(QP_ket)
       nspace=LEN(TRIM(QP_ket))

       counter=0 ; KS_row=REPEAT('',nspace+2)
       do ib_KS=1,Bst%mband
         if (ABS(m_lda_to_qp(ib_KS,ib_QP,ikibz,isp))<my_tolmat) CYCLE
         counter=counter+1
         call int2char10(ib_KS,bks)
         write(tmpstr,'(3a)')' |',TRIM(bks),'>'

         if (use_rhophi) then
           ! coefficient as (rho, phi)
           cx(1) = real(m_lda_to_qp(ib_KS,ib_QP,ikibz,isp))
           cx(2) = aimag(m_lda_to_qp(ib_KS,ib_QP,ikibz,isp))
           call rhophi(cx, phi, rho)
           write(KS_ket,'(1x,2f7.3,a,1x)')rho, phi, TRIM(tmpstr)
         else
           ! coefficient as (Re, Im)
           write(KS_ket,'(1x,2f7.3,a,1x)')m_lda_to_qp(ib_KS,ib_QP,ikibz,isp),TRIM(tmpstr)
         end if
         KS_row=TRIM(KS_row)//TRIM(KS_ket)
         if (MOD(counter,NBRA)==0) then  ! nbra KS kets per row
           write(my_unt,'(a)')TRIM(KS_row)
           KS_row=REPEAT('',nspace+2)
         end if
       end do

       if (MOD(counter,NBRA)/=0) write(my_unt,'(a)')TRIM(KS_row) ! Last row, if any
       write(my_unt,'(a)')''
     end do !ib_QP

   end do !ikibz
 end do !isp

 write(my_unt,'(a,/)')' '//REPEAT('*',76)

end subroutine show_QP
!!***

!----------------------------------------------------------------------

!!****f* m_qparticles/rdgw
!! NAME
!! rdgw
!!
!! FUNCTION
!!  This subroutine reads the GW corrections from a _GW file.
!!
!! INPUTS
!!  [extrapolate]= if .TRUE., the routine extrapolates the
!!    GW corrections for the states that have not been explicitly evaluated (default).
!!    If .FALSE., only the GW states that have been calculated will be used to replace
!!    the input eigenvalues stored in Bst%eig
!!  Bst<ebands_t>=type describing the Band structure.
!!    %nbnds=number of bands.
!!    %nkpt=number of irred k-points.
!!    %nsppol=number of spin
!!    %kptns(3,nkpt)=irreducible k-points
!!
!! SIDE EFFECTS
!!   Bst%eig(%mband,%nkpt,%nsppol)=Overwritten with GW energies according to extrapolate flag.
!!
!! OUTPUT
!!   igwene(Bst%mband,Bst%nkpt,Bst%nsppol)= The imaginary part of the QP energies.
!!
!! PARENTS
!!      mlwfovlp_qp,screening,setup_bse,sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine rdgw(Bst,fname,igwene,extrapolate)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: fname
 logical,optional,intent(in) :: extrapolate
 type(ebands_t),intent(inout) :: Bst
!arrays
 real(dp),intent(out) :: igwene(Bst%mband,Bst%nkpt,Bst%nsppol)

!Local variables ------------------------------
!scalars
 integer :: ib,ibr,ik,ikibz,ikr,is,nn,nbandR,nkibzR,nsppolR,unt,nbv
 real(dp) :: alpha,beta,degw,egw_r,egw_i,smrt
 logical :: do_extrapolate
 character(len=500) :: msg
!arrays
 integer,allocatable :: vbik(:,:),seen(:)
 real(dp) :: kread(3)
 real(dp),allocatable :: gwcorr(:,:,:)

!************************************************************************

 call wrtout(std_out,'Reading GW corrections from file: '//TRIM(fname),'COLL')
 ABI_CHECK(ALL(Bst%nband==Bst%mband),"nband must be constant")

 if (open_file(fname,msg,newunit=unt,status='old') /=0) then
   MSG_ERROR(msg)
 end if

 read(unt,*)nkibzR,nsppolR

 ABI_CHECK(nsppolR==Bst%nsppol,"mismatch in nsppol")
 if (nkibzR/=Bst%nkpt) then
   write(msg,'(a,i4,a,i4,2a)')&
&   'Found less k-points than that required ',nkibzR,'/',Bst%nkpt,ch10,&
&   'Some k-points will be skipped. Continuing anyway '
   MSG_WARNING(msg)
 end if

 ABI_MALLOC(gwcorr,(Bst%mband,Bst%nkpt,Bst%nsppol))
 ABI_MALLOC(seen,(Bst%nkpt))
 gwcorr=zero
 igwene=zero

 do is=1,Bst%nsppol
   seen=0

   do ikr=1,nkibzR
     read(unt,*)kread(:)
     read(unt,*)nbandR
     ikibz=0
     do ik=1,Bst%nkpt
       if (ALL(ABS(kread(:)-Bst%kptns(:,ik))<0.0001)) then
         ikibz=ik
         seen(ik) = seen(ik) + 1
       end if
     end do
     do ib=1,nbandR
       read(unt,*)ibr,egw_r,degw,egw_i
       if (ibr<=Bst%mband .and. ikibz/=0) then
         gwcorr(ibr,ikibz,is)=degw/Ha_eV
         igwene(ibr,ikibz,is)=egw_i/Ha_eV
       end if
     end do
   end do

   if (ANY(seen/=1)) then
     do ik=1,Bst%nkpt
       if (seen(ik)/=1) then
         write(msg,'(a,3f8.3,a)')" k-point: ",Bst%kptns(:,ik)," not found in the GW file!"
         MSG_WARNING(msg)
       end if
     end do
   end if

 end do

 ABI_FREE(seen)
 close(unt)

 do_extrapolate=.TRUE.; if (PRESENT(extrapolate)) do_extrapolate=extrapolate

 if (.not. do_extrapolate) then ! Only the bands calculated are updated.
   Bst%eig = Bst%eig + gwcorr

 else

   if (ANY(ABS(igwene)>tol6)) then
     write(msg,'(4a)')ch10,&
&      "The GW file contains QP energies with non-zero imaginary part",ch10,&
&      "Extrapolation not coded, change the source! "
     MSG_ERROR(msg)
   end if

   ABI_MALLOC(vbik,(BSt%nkpt,BSt%nsppol))
   vbik(:,:) = get_valence_idx(BSt)

   do is=1,Bst%nsppol
     do ik=1,Bst%nkpt

      nbv=vbik(ik,is) ! Index of the (valence band| Fermi band) for each spin
      nn=Bst%mband-nbv

      do ib=nbv+1,Bst%mband
        if ( ABS(gwcorr(ib,ik,is)) < tol16) then
          nn=ib-1-nbv
          if (nn>1) then
            call wrtout(std_out,"Linear extrapolating (conduction) GW corrections beyond the read values","COLL")
            smrt=linfit(nn,Bst%eig(nbv+1:nbv+nn,ik,is),gwcorr(nbv+1:nbv+nn,ik,is),alpha,beta)
          else
            call wrtout(std_out,"Assuming constant (conduction) GW corrections beyond the read values",'COLL')
            alpha=zero
            beta =gwcorr(nbv+nn,ik,is)
          end if
          EXIT !ib loop
        end if
      end do !ib

      do ib=nbv+nn+1,Bst%mband
        gwcorr(ib,ik,is)= alpha*Bst%eig(ib,ik,is) + beta
      end do

      nn=nbv
      do ib=nbv,1,-1
        if ( ABS(gwcorr(ib,ik,is)) < tol16) then
         nn=nbv-ib
         if (nn>1) then
           call wrtout(std_out,"Linear extrapolating (valence) GW corrections beyond the read values","COLL")
           smrt=linfit(nn,Bst%eig(nbv-nn+1:nbv,ik,is),gwcorr(nbv-nn+1:nbv,ik,is),alpha,beta)
         else
           call wrtout(std_out,"Assuming constant (valence) GW corrections beyond the read values","COLL")
           alpha=zero
           beta =gwcorr(nbv,ik,is)
         end if
         EXIT !ib
        end if
      end do !ib

      do ib=1,nbv-nn
        gwcorr(ib,ik,is)=alpha*Bst%eig(ib,ik,is) + beta
      end do

     end do !ik
   end do !is

   call wrtout(std_out,' k  s     GW corrections [eV] ','COLL')
   do is=1,Bst%nsppol
     do ik=1,Bst%nkpt
       write(msg,'(i3,1x,i3,10f7.2/50(10x,10f7.2/))')ik,is,(Ha_eV*gwcorr(ib,ik,is),ib=1,Bst%mband)
       call wrtout(std_out,msg,"COLL")
     end do
   end do
   Bst%eig = Bst%eig + gwcorr
   ABI_FREE(vbik)
 end if

 call wrtout(std_out,' k   s    GW eigenvalues [eV]',"COLL")
 do is=1,Bst%nsppol
   do ik=1,Bst%nkpt
     write(std_out,'(2(i3,1x),7x,10f7.2/50(15x,10f7.2/))')ik,is,(Ha_eV*Bst%eig(ib,ik,is),ib=1,Bst%mband)
   end do
 end do

 ABI_FREE(gwcorr)

end subroutine rdgw
!!***

!----------------------------------------------------------------------

!!****f* m_qparticles/updt_m_lda_to_qp
!! NAME
!! updt_m_lda_to_qp
!!
!! FUNCTION
!! Updates the matrix containing the unitary transformation from the lda states
!! to the quasiparticle states.
!!
!! INPUTS
!!  Sigp<sigparams_t>=Parameters characterizing the self-energy calculation.
!!     %nsppol=1 for unpolarized, 2 for spin-polarized
!!     %nbnds=number of bands used for sigma
!!  Sr<sigma_t>=Structure containing the results of the sigma run.
!!     %en_qp_diago(nbnds,nibz,nsppol)= NEW quasi-particle energies
!!     %eigvec_qp(nbnds,nbnds,nibz,nsppol)= NEW QP amplitudes in the KS basis set
!!      obtained by diagonalizing H0 + Herm(Sigma).
!!  Kmesh<kmesh_t>=information on the k-point sampling.
!!     %nibz=number of irreducible k-points
!!     %ibz(3,kibz)=reduced coordinates of the irreducible k-points
!!  nscf=Number of self consistent cycles performed
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  m_lda_to_qp(nbnds,nbnds,nibz,nsppol)= overwritten with the new QP amplitudes
!!                                        in terms of KS wavefunctions
!!
!! NOTES
!!  Only master node should call this routine.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine updt_m_lda_to_qp(Sigp,Kmesh,nscf,Sr,m_lda_to_qp)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nscf
 type(kmesh_t),intent(in) :: Kmesh
 type(sigparams_t),intent(in) :: Sigp
 type(sigma_t),intent(in) :: Sr
!arrays
 complex(dpc),intent(inout) :: m_lda_to_qp(Sigp%nbnds,Sigp%nbnds,Kmesh%nibz,Sigp%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ik,is
!arrays
 complex(dpc),allocatable :: mtmp(:,:)

! *************************************************************************

 if (nscf>=0) then
   ! Calculate the new m_lda_to_qp
   ABI_MALLOC(mtmp,(Sigp%nbnds,Sigp%nbnds))
   do is=1,Sigp%nsppol
     do ik=1,Kmesh%nibz
       mtmp(:,:)=m_lda_to_qp(:,:,ik,is)
       m_lda_to_qp(:,:,ik,is)=MATMUL(mtmp(:,:),Sr%eigvec_qp(:,:,ik,is))
     end do
   end do
   ABI_FREE(mtmp)
 end if

end subroutine updt_m_lda_to_qp

!----------------------------------------------------------------------

END MODULE m_qparticles
!!***
