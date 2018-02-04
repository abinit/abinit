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

/* Returns in algn the alignment of a pointer modulo nbytes NONPORTABLE */
void 
FC_FUNC_(clib_alignment_of,clib_ALIGNMENT_OF)
  (void *p, int *nbytes, int *algn)
{
  if (*nbytes){
    *algn = (int)(((uintptr_t) p) % *nbytes);
     /*printf("address = %p; uintptr = %d\n", p, (uintptr_t) p);*/
  }
  else 
    *algn = 0;
}
