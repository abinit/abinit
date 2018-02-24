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
#include "string_f.h"

/* 
 * This function causes any buffered output on stream to be delivered to the file. 
 * If stream is a null pointer, then fflush causes buffered output on all open 
 * output streams to be flushed. It returns EOF if a write error occurs, or zero otherwise.
 */
void 
FC_FUNC_(clib_fflush,CLIB_FFLUSH)
  (int* ierr)
{
  *ierr = fflush(NULL);
}

/* This function returns a string that is the value of the environment variable fname
 * If the environment variable fname is not defined, the value is a null pointer. 
 */
void 
FC_FUNC_(clib_getenv,CLIB_GETENV)
  (int* ierr, STR_F_TYPE fname, STR_F_TYPE fvalue STR_ARG2)
{
  char *cname, *cvalue;
  TO_C_STR1(fname, cname);
  cvalue = getenv (cname);
  xfree(cname);

  *ierr = 1;
  if (cvalue && strlen(cvalue) <= l2 ){
    *ierr = 0;
    TO_F_STR2(cvalue, fvalue);
  } 
}

/*
*  Replacement for the F2003 intrinsic GET_ENVIRONMENT_VARIABLE.
*  This version can be called in Fortran code but does not support optional arguments.
*
* CALL GET_ENVIRONMENT_VARIABLE (NAME[,VALUE,LENGTH,STATUS,TRIM_NAME]) 
*   obtains the value of an environment variable.
*
*   NAME is a default character INTENT(IN) scalar that identifies the required environment
*          variable. The interpretation of case is processor dependent.
*   [VALUE] is a default character INTENT(OUT) scalar that is assigned the value of the environment variable.
*   [LENGTH]is a default integer INTENT(OUT) scalar. If the specified environment variable exists and has a value, 
*           LENGTH is set to the length (number of characters) of that value. Otherwise, LENGTH is set to 0.
*   [STATUS] is a default integer INTENT(OUT) scalar that indicates success or failure.
*   [TRIM_NAME] is a logical INTENT(IN) scalar that indicates whether trailing blanks
*               in NAME are considered significant.
*/

void 
FC_FUNC_(clib_get_environment_variable,CLIB_GET_ENVIRONMENT_VARIABLE)
  (const STR_F_TYPE fname, STR_F_TYPE fvalue, int *length, int *status, const int *trim_name STR_ARG2)
  /* name_len, value_len) */
{
  char *cname, *cvalue;
  int q_len;

  if (*trim_name)
    TO_C_STR1(fname, cname)
  else 
    TO_C_STR1_NOTRIM(fname, cname)

  cvalue = getenv(cname);
  q_len = (cvalue == NULL) ? 0 : strlen(cvalue);
  xfree(cname);

  *length = (cvalue == NULL) ? 0 : q_len;

  if (cvalue == NULL)  
    *status = 1;  /* $fname does not exist */

  else if (cvalue != NULL && l2 < q_len)  
    *status = -1; /* Fortran string is too short */

  else
    *status = 0;  /* OK */

  if (cvalue != NULL) TO_F_STR2(cvalue, fvalue)
}
