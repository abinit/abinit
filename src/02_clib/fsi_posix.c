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

#include "abi_clib.h"
#include "xmalloc.h"
#include "string_f.h"

/* Create a directory */
void 
FC_FUNC_(clib_mkdir,CLIB_MKDIR)
  (int *ierr, STR_F_TYPE fname STR_ARG1)
{
  struct stat buf;
  char* cname;
  TO_C_STR1(fname, cname);

  *ierr = 1;
  if (!*cname) return; 

  if (stat(cname, &buf) == 0) {
    xfree(cname);
    *ierr=2; return;
  }

  *ierr = mkdir(cname, 0775);
  if (*ierr == -1) fprintf (stderr, "Cannot create dir %s; %s\n", cname, strerror(errno) );
  xfree(cname);
}
                 
/* change to a directory */
void 
FC_FUNC_(clib_chdir,CLIB_CHDIR)
  (int *ierr, STR_F_TYPE fname STR_ARG1)
{
  char* cname;
  TO_C_STR1(fname, cname);
  *ierr = chdir(cname);
  if (*ierr == -1) fprintf (stderr, "Couldn't change to dir %s; %s\n",cname, strerror (errno));
  xfree(cname);
}

/* Rename a directory */
void 
FC_FUNC_(clib_rename,CLIB_RENAME)
  (int *ierr, STR_F_TYPE from_fname, STR_F_TYPE to_fname STR_ARG2)
{
  char* from_cname;
  char* to_cname;

  TO_C_STR1(from_fname, from_cname);
  TO_C_STR1(to_fname, to_cname);
  *ierr = rename(from_cname, to_cname);

  if (*ierr == -1) fprintf (stderr, "Couldn't rename file dir %s; %s\n",from_cname, strerror (errno));
  xfree(from_cname);
  xfree(to_cname);
}

/* Delete a directory */
void 
FC_FUNC_(clib_remove, CLIB_REMOVE)
  (int *ierr, STR_F_TYPE fname STR_ARG1)
{
  char* cname;
  TO_C_STR1(fname, cname);
  *ierr = remove(cname);
  if (*ierr == -1) fprintf (stderr, "Couldn't remove file %s; %s\n",cname, strerror (errno));
  xfree(cname);

}

/* GNU's getcwd (NULL, 0) using the standard behavior of getcwd */
static char* gnu_getcwd(void);

static char*
gnu_getcwd (void)
{
  size_t size = 100;
  while (1){
    char* buffer = (char*) xmalloc (size);
    if ( (char*) getcwd (buffer, size) == buffer) return buffer;
    xfree (buffer);
    if (errno != ERANGE) return NULL;
    size *= 2;
  }
}

/* Return the name of the current working directory */
void 
FC_FUNC_(clib_getcwd, CLIB_GETCWD)
  (int* ierr, STR_F_TYPE fname STR_ARG1)
{
  int ln; 
  char* cname = gnu_getcwd();
  ln = strlen(cname);
  *ierr = 0;

  if (ln > l1){
    fprintf(stderr, " Dir name %s will be truncated\n", cname);
    *ierr= 1;
  }
  TO_F_STR1(cname, fname);
  xfree(cname);
}

static char* gnu_gethostname(void);

static char*
gnu_gethostname (void)
{
  size_t size = 100;
  while (1){
    char *buffer = (char *) xmalloc(size);
    if ( gethostname(buffer, size) == 0) return buffer;
    xfree (buffer);
    if (errno != ENAMETOOLONG) return NULL;
    size *= 2;
  }
}

/* Returns hostname */
void 
FC_FUNC_(clib_gethname, CLIB_GETHNAME)
  (int* ierr, STR_F_TYPE fname STR_ARG1)
{
  int ln;
  char* cname = gnu_gethostname(); 

  *ierr = 0;
  if (cname == NULL){
    *ierr = ERANGE; /*check this*/
    return;
  }

  ln = strlen(cname);
  if (ln > l1){
    *ierr = ERANGE; /*check this*/
    fprintf(stderr, "hostname %s will be truncated\n", cname);
  }
  
  TO_F_STR1(cname, fname);
  xfree(cname);
}
