  !!****m* ABINIT/m_multibinit_main2
  !! NAME
  !! m_multibinit_main2
  !!
  !! FUNCTION
  !! multibinit_main2 is the main function of multibinit. It runs after
  !! the filenames from files file is read. It then does everything else.
  !!
  !! Datatypes:
  !!
  !!
  !! Subroutines:
  !! multibinit_main2
  !!
  !!
  !! COPYRIGHT
  !! Copyright (C) 2001-2019 ABINIT group (hexu)
  !! This file is distributed under the terms of the
  !! GNU General Public License, see ~abinit/COPYING
  !! or http://www.gnu.org/copyleft/gpl.txt .
  !! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
  !!
  !! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
module m_multibinit_main2
  use defs_basis
  use defs_abitypes
  use m_build_info
  use m_xmpi
  use m_xomp
  use m_abicore
  use m_errors
  use m_multibinit_manager, only: mb_manager_t

  implicit none

  !!***
contains

  !!****f* m_multibinit_main2/multibinit_main2
  !!
  !! NAME
  !! multibinit_main2
  !!
  !! FUNCTION
  !! The main function 
  !!
  !! INPUTS
  !! filnam: file names from files file
  !!
  !! OUTPUT
  !!
  !! PARENTS
  !!
  !!
  !! CHILDREN
  !!
  !!
  !! SOURCE
  subroutine multibinit_main2(filnam)
    character(len=fnlen), intent(inout) :: filnam(17)
    type(mb_manager_t) :: manager
    call manager%run_all(filnam)
  end subroutine multibinit_main2
  !!***

end module m_multibinit_main2
