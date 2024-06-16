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

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <errno.h>
#include "abi_clib.h"
#include "xmalloc.h"

/* Return 0 if dirpath is not an existent directory */
int _directory_exist(const char *dirpath) {
    struct stat info;

    if (stat(dirpath, &info) != 0) {
        if (errno == ENOENT) {
            return 0; /* Path does not exist */
        } else {
            return -1; /* Other error (e.g., permission denied) */
        }
    } else if (info.st_mode & S_IFDIR) {
        return 1; /* Path is a directory */

    } else {
        return 2; /* Path exists but is not a directory */
    }
}

/* Create directory dirpath if it does not exist */
void c_mkdir_if_needed(char *dirpath, int *ierr)
{
   *ierr = _directory_exist(dirpath);

   if (*ierr == 0) {
     /* S_IRWXU read, write, execute/search by owner: http://pubs.opengroup.org/onlinepubs/7908799/xsh/sysstat.h.html */
     *ierr = mkdir(dirpath, S_IRWXU);
   }
   else if (*ierr == 1) {
     *ierr = 0;
   }
}
