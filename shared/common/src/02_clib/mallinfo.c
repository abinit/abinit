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

void clib_mallinfo
  (long int *arena, long int *hblkhd, long int *usmblks, long int *fsmblks, long int *uordblks, long int *fordblks)
{
#ifdef HAVE_MALLINFO
  struct mallinfo info = mallinfo();
  *arena    = info.arena;
  *hblkhd   = info.hblkhd;
  *usmblks  = info.usmblks;
  *fsmblks  = info.fsmblks;
  *uordblks = info.uordblks;
  *fordblks = info.fordblks;
#else
  *arena    = -1.0; 
  *hblkhd   = -1.0;
  *usmblks  = -1.0;
  *fsmblks  = -1.0;
  *uordblks = -1.0;
  *fordblks = -1.0;
#endif
}


/* SVID2/XPG mallinfo structure */
#if 0
struct mallinfo {
  int arena;    /* non-mmapped space allocated from system */
  int ordblks;  /* number of free chunks */
  int smblks;   /* number of fastbin blocks */
  int hblks;    /* number of mmapped regions */
  int hblkhd;   /* space in mmapped regions */
  int usmblks;  /* maximum total allocated space */
  int fsmblks;  /* space available in freed fastbin blocks */
  int uordblks; /* total allocated space */
  int fordblks; /* total free space */
  int keepcost; /* top-most, releasable (via malloc_trim) space */
};
#endif
