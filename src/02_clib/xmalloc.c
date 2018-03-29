/* xmalloc.c -- safe versions of malloc, realloc, free */
/*
 * Copyright (C) 2009-2018 ABINIT group (MG)
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

static void
memory_error_and_abort (const char *fname)
{
  fprintf (stderr, "%s: out of virtual memory\n", fname);
#ifdef HAVE_ABORT
  abort();
#else
  exit(EXIT_FAILURE);
#endif
}

/* Return a pointer to free()able block of memory large enough
   to hold BYTES number of bytes.  If the memory cannot be allocated,
   print an error message and abort. */
void*
xmalloc(size_t bytes)
{
  void* temp;

  temp = malloc (bytes);
  if (temp == NULL)
    memory_error_and_abort ("xmalloc");
  return (temp);
}

void*
xrealloc (void *pointer, size_t bytes)
{
  void* temp;

  temp = pointer ? realloc (pointer, bytes) : malloc (bytes);

  if (temp == NULL)
    memory_error_and_abort ("xrealloc");
  return (temp);
}

/* Use this as the function to call when adding unwind protects so we
   don't need to know what free() returns. */
void
xfree (void* ptr)
{
  if (ptr)
    free (ptr);
}

#if 0
typedef void* voidptr_t;

/* Variadic function for freeing several pointers. Use at your own risk */
void
vxfree (const int n, ...)
{
  int i;
  void* ptr;
  va_list vl;

  va_start(vl, n);

  for (i=0; i<n; i++){
    ptr = va_arg(vl, voidptr_t );
    xfree(ptr);
  }

  va_end(vl);
}
#endif
