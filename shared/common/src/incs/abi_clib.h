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


/* ABINIT internal header file */
#ifndef _ABINIT_CLIB_H
#define _ABINIT_CLIB_H

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_STDIO_H
#include <stdio.h>
#endif

#ifdef HAVE_STDLIB_H
#include <stdlib.h>		/* size_t */
#endif

#ifdef HAVE_STDARG_H
#include <stdarg.h>		/* va_list */
#endif

#ifdef HAVE_STDDEF_H
#include <stddef.h>             /* ptrdiff_t */
#endif

#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#ifdef HAVE_STDINT_H
#include <stdint.h>             /* uintptr_t, maybe */
#endif

#ifdef HAVE_INTTYPES_H
#include <inttypes.h>           /* uintptr_t, maybe. C99 */
#endif

#ifdef HAVE_STRING_H
#include <string.h>
#endif

#ifdef HAVE_ERRNO_H
#include <errno.h>
#endif

#ifdef HAVE_MALLOC_H
#include <malloc.h>
#elif HAVE_SYS_MALLOC_H
#include <sys/malloc.h>
#endif

#ifdef HAVE_MATH_H
#include <math.h>
#endif

#ifdef HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#define STRINGIZEx(x) #x
#define STRINGIZE(x) STRINGIZEx(x)
#define CONCAT(prefix, name) prefix ## name

/* determine precision */ 
#if defined(ABINIT_SINGLE)
  typedef float R;
#elif defined(ABINIT_LDOUBLE)
  typedef long double R;
#else
  typedef double R;
#endif

/* dummy use of unused parameters to silence compiler warnings */
#define UNUSED(x) (void)x

#if 0
/* inline version */
#define IABS(x) (((x) < 0) ? (0 - (x)) : (x))

/*-----------------------------------------------------------------------*/
/* assert.c: */
void assertion_failed (const char *s, const int line, const char *file);

/* always check */
#define ASSERT(ex)						 \
      (void) ( (ex) || (assertion_failed (#ex, __LINE__, __FILE__), 0) )

#ifdef DEBUG_MODE      /* check only if debug enabled */
#define DBG_ASSERT(ex)						 \
      (void)( (ex) || (assertion_failed (#ex, __LINE__, __FILE__), 0) )
#else
#define DBG_ASSERT(ex) /* nothing */
#endif

extern void debug(const char *format, ...);
#define DBG debug
#endif

/*-----------------------------------------------------------------------*/
/* xmalloc.c: */

void* xmalloc(size_t bytes);
/* void* xrealloc(void* pointer, size_t bytes); */
void  xfree(void *ptr);

/*-----------------------------------------------------------------------*/
/* xexit.c */
void xexit(int code);

#endif /* _ABINIT_CLIB_H */
