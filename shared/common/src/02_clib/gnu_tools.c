/*
 * Copyright (C) 2009-2020 ABINIT group (MG)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

#include "abi_clib.h"

#ifdef HAVE_MCHECK_H
#include <mcheck.h>
#endif 

/* - Function: void mtrace (void)
 * When the mtrace function is called it looks for an environment variable named MALLOC_TRACE. 
 * This variable is supposed to contain a valid file name. The user must have write access. 
 * If the file already exists it is truncated. If the environment variable is not set or it does 
 * not name a valid file which can be opened for writing nothing is done. The behavior of malloc etc. is not changed. 
 * If the named file is successfully opened, mtrace installs special handlers for the functions malloc, realloc, and free 
 * (see Hooks for Malloc). From then on, all uses of these functions are traced and protocolled into the file. 
 * There is now of course a speed penalty for all calls to the traced functions so tracing should not be enabled during 
 * normal use. 
 * This function is a GNU extension and generally not available on other systems. The prototype can be found in mcheck.h. 
 */
void clib_mtrace (int *ierr)
{
#ifdef HAVE_MCHECK_H
  mtrace();
  *ierr=0;
#else
  *ierr=1;
#endif
}

/* - Function: void muntrace (void)
 * The muntrace function can be called after mtrace was used to enable tracing the malloc calls. 
 * If no (successful) call of mtrace was made muntrace does nothing. 
 * Otherwise it deinstalls the handlers for malloc, realloc, and free and then closes the protocol file. 
 * No calls are protocolled anymore and the program runs again at full speed. 
 * Please note that not only the application uses the traced functions, also libraries (including the C library itself) 
 * use these functions.
 * 
 * It is no good idea to call muntrace before the program terminated. The libraries are informed about 
 * the termination of the program only after the program returns from main or calls exit and so cannot 
 * free the memory they use before this time. 
 * 
 * So the best thing one can do is to call mtrace as the very first function in the program and never call muntrace. 
 * So the program traces almost all uses of the malloc functions (except those calls which are executed by 
 * constructors of the program or used libraries).
 * 
 * This function is a GNU extension and generally not available on other systems. The prototype can be found in mcheck.h.
 */
void clib_muntrace (int* ierr)
{
#ifdef HAVE_MCHECK_H
  muntrace();
  *ierr=0;
#else
  *ierr=1;
#endif
}

/*  
 * You can ask malloc to check the consistency of dynamic memory by using the mcheck function. 
 * This function is a GNU extension, declared in mcheck.h.
 * Calling mcheck tells malloc to perform occasional consistency checks. These will catch things such as 
 * writing past the end of a block that was allocated with malloc.
 * It is too late to begin allocation checking once you have allocated anything with malloc. 
 * So mcheck does nothing in that case. The function returns -1 if you call it too late, 
 * and 0 otherwise (when it is successful).
 * function prototype: mcheck (void (*abortfn) (enum mcheck_status status))
 */
void clib_mcheck (int* ierr)
{
#ifdef HAVE_MCHECK_H
  *ierr = mcheck (NULL);
#else
  *ierr = 2;
#endif
}

