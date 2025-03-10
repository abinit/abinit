! **************************************************************************************************
!  Copyright (C) 2020-2023 Green-X library
!  This file is distributed under the terms of the APACHE2 License.
!
! **************************************************************************************************
/*
#ifdef HAVE_FC_LONG_LINES
#define _REGISTER_EXC(msg) call register_exc(msg, filename=__FILE__, lineno=__LINE__) 
#else
#define _REGISTER_EXC(msg) call register_exc(msg)
#endif
*/

#define _REGISTER_EXC(msg) ABI_ERROR(msg)
