!{\src2tex{textfont=tt}}
!!****p* ABINIT/bsepostproc
!! NAME
!! bsepostproc
!!
!! FUNCTION
!!  Utility for post-processing Bethe-Salpeter results
!!
!! COPYRIGHT
!! Copyright (C) 2013-2018 ABINIT group (YG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  (main program)
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,close_haydock,continued_fract,destroy_mpi_enreg
!!      flush_unit,herald,initmpi_seq,open_haydock,read_dim_haydock
!!      read_haydock,timein,wrtout,xmpi_end,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 program bsepostproc

 use defs_basis
 use defs_abitypes
 use m_build_info
 use m_xmpi
 use m_haydock_io
 use m_numeric_tools

 use m_time,      only : timein
 use m_specialmsg,only : specialmsg_getcount, herald
 use m_io_tools,  only : get_unit, flush_unit
 use m_mpinfo,    only : destroy_mpi_enreg, nullify_mpi_enreg, initmpi_seq

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bsepostproc'
!End of the abilint section

 implicit none

!Arguments ----------------------------
!Local variables-----------------------
!scalars
 integer :: funt
 integer :: ios
 integer :: niter_file, io, iq
 integer :: n_all_omegas, nomega
 integer :: term
 real(dp) :: broad_in
 real(dp) :: omega_min,omega_max,delta_omega
 real(dp) :: omegaev
 real(dp) :: tcpui,twalli
 complex(dpc) :: factor
 character(len=50) :: restart_file
 character(len=500) :: frm
 character(len=24) :: codename
 character(len=50) :: output_file
 type(haydock_type) :: haydock_file
 type(MPI_type) :: mpi_enreg
!arrays
 real(dp),allocatable :: tmp_eps(:,:)
 real(dp),allocatable :: bb_file(:)
 complex(dpc),allocatable :: omega(:),green_temp(:),green(:,:)
 complex(dpc),allocatable :: aa_file(:),phi_n_file(:),phi_nm1_file(:)
 complex(dpc),allocatable :: all_omegas(:)

!*******************************************************

!Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

!Initialize MPI
 call xmpi_init()

 call timein(tcpui,twalli)

!Default for sequential use
 call initmpi_seq(mpi_enreg)

 codename='BSEPOSTPROC'//REPEAT(' ',13)
 call herald(codename,abinit_version,std_out)

 write(std_out,'(a)') "Broad_in (eV) ?"
 read(std_in,*) broad_in

 write(std_out,'(a)') "Range of frequencies (eV)"
 read(std_in,*) omega_min, omega_max, delta_omega

 write(std_out,'(a)') "Terminator ?"
 read(std_in,*) term

 write(std_out,'(a)') "Input file ?"
 read(std_in,*) restart_file

 write(std_out,'(a)') "Output file ?"
 read(std_in,*) output_file

 omega_min = omega_min/Ha_eV
 omega_max = omega_max/Ha_eV
 delta_omega = delta_omega/Ha_eV

 broad_in = broad_in/Ha_eV

 nomega = (omega_max - omega_min)/delta_omega + 1
 ABI_MALLOC(omega,(nomega))
 do io=1,nomega
   omega(io) = (omega_min + (io-1)*delta_omega)  + j_dpc*broad_in
 end do

!Create new frequencies "mirror" in negative range to add
!their contributions. Can be improved by computing only once
!zero frequency, but loosing clearness
 n_all_omegas = 2*nomega

 ABI_MALLOC(all_omegas,(n_all_omegas))
!Put all omegas with frequency > 0 in table
 all_omegas(nomega+1:n_all_omegas) = omega
!Put all omegas with frequency < 0
!Warning, the broadening must be kept positive
 all_omegas(1:nomega) = -DBLE(omega(nomega:1:-1)) &
& + j_dpc*AIMAG(omega(nomega:1:-1))

 ABI_MALLOC(green_temp,(n_all_omegas))

 call open_haydock(restart_file,haydock_file)

 call read_dim_haydock(haydock_file)

 ABI_MALLOC(green,(nomega,haydock_file%nq))

 do iq = 1, haydock_file%nq
   call read_haydock(haydock_file,haydock_file%qpoints(:,iq),aa_file,bb_file,phi_nm1_file,phi_n_file, niter_file,factor)

   call continued_fract(niter_file,term,aa_file,bb_file,n_all_omegas,all_omegas,green_temp)

!  Computing result from two ranges of frequencies
!  The real part is added, the imaginary part is substracted
   green(:,iq) = green_temp(nomega+1:n_all_omegas)+CONJG(green_temp(nomega:1:-1))

   green(:,iq) = cone+factor*green(:,iq)

 end do

 call close_haydock(haydock_file)

 ABI_MALLOC(tmp_eps,(2,haydock_file%nq))

 funt = get_unit()
 open(unit=funt,file=output_file,form="formatted",iostat=ios)

 write(funt,'(a)')"# omega [eV]    RE(eps(q=1)) IM(eps(q=1) RE(eps(q=2) ) ... "
!write(frm,*)'(f7.3,',2*BSp%nq,'es12.4)'
 write(frm,*)'(f7.3,',2*haydock_file%nq,'(1x,f9.4))'
 do io=1,nomega
   omegaev = DBLE(omega(io))*Ha_eV
   tmp_eps(1,:) = REAL (green(io,:))
   tmp_eps(2,:) = AIMAG(green(io,:))
!  where (ABS(tmp_eps) < SMALL) ! this to improve the portability of the automatic tests.
!  tmp_eps = zero
!  end where
   write(funt,frm) omegaev,(tmp_eps(:,iq), iq=1,haydock_file%nq)
 end do

 ABI_FREE(tmp_eps)

 close(funt)

 ABI_FREE(omega)
 ABI_FREE(all_omegas)
 ABI_FREE(green_temp)

 call wrtout(std_out,ch10//" Analysis completed.","COLL")

 call flush_unit(std_out)

 call destroy_mpi_enreg(mpi_enreg)
 call xmpi_end()

 end program bsepostproc
!!***
