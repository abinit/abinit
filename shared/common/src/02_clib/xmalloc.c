/* xmalloc.c -- safe versions of malloc, realloc, free */
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
#include "xmalloc.h"

static long int COUNT_MALLOC = 0, COUNT_FREE = 0;
static size_t BYTES_ALLOCATED = 0;
/*
 * Counting bytes deallocated by free requires some non-standard API to access the status of the mallocator
*/

static void memory_error_and_abort (const char *fname)
{
  fprintf (stderr, "%s: out of virtual memory\n", fname);
#ifdef HAVE_ABORT
  abort();
#else
  exit(EXIT_FAILURE);
#endif
}

/* Return a pointer to free()able block of memory large enough
   to hold BYTES number of bytes. If the memory cannot be allocated,
   print an error message and abort. */

void* xmalloc(size_t bytes)
{
  void* temp;

  temp = malloc (bytes);
  if (temp == NULL)
    memory_error_and_abort("xmalloc");

#ifdef HAVE_MEM_PROFILING
  COUNT_MALLOC += 1;
  BYTES_ALLOCATED += bytes;
#endif

  return temp;
}

/*
void* xrealloc (void *pointer, size_t bytes)
{
  void* temp;
  temp = pointer ? realloc (pointer, bytes) : malloc (bytes);

  if (temp == NULL)
    memory_error_and_abort ("xrealloc");
  return (temp);
}
*/

/* Free pointer returned by xmalloc */

void xfree (void* ptr)
{
  if (ptr){
    free(ptr);
#ifdef HAVE_MEM_PROFILING
    COUNT_FREE += 1;
#endif
  }
}
