/*
 * Copyright (C) 2009-2019 ABINIT group (MG)
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

/* This variable is set by xatexit if it is called.  This way, xmalloc
   doesn't drag xatexit into the link.
void (*_xexit_cleanup) ((void));
*/


/*
Terminates the program.  If any functions have been registered with
the @code{xatexit} replacement function, they will be called first.
Termination is handled via the system's normal @code{exit} call.
*/
void
xexit (int code)
{
 /* if (_xexit_cleanup != NULL)
    (*_xexit_cleanup) ();*/
  exit (code);
}
