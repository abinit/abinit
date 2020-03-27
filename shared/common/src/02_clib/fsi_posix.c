/* Fortran callable procedures compliant to the File System Interface defined in POSIX version? */
/*
 Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 02111-1307, USA.
*/

#include <sys/stat.h>
#include "abi_clib.h"
#include "xmalloc.h"


void c_mkdir(char *path, int *ierr)
{
   /* S_IRWXU read, write, execute/search by owner
   http://pubs.opengroup.org/onlinepubs/7908799/xsh/sysstat.h.html
   */
   *ierr = mkdir(path, S_IRWXU);
}
