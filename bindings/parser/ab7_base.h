/** -*- C header file -*-
 *
 * Copyright (C) 2009-2018 ABINIT Group (Damien Caliste)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#ifndef AB7_BASE
#define AB7_BASE

/**
 * Ab7Error:
 * @AB7_NO_ERROR: no error.
 * @AB7_ERROR_INVARS_OBJ: wrong dataset object.
 * @AB7_ERROR_INVARS_ATT: wrong attribute in dataset.
 * @AB7_ERROR_INVARS_ID: wrong dataset index.
 * @AB7_ERROR_INVARS_SIZE: wrong size when accessing arrays.
 *
 * An error code.
 */
typedef enum
  {
    AB7_NO_ERROR,
    AB7_ERROR_OBJ,
    AB7_ERROR_ARG,
    AB7_ERROR_INVARS_ATT,
    AB7_ERROR_INVARS_ID,
    AB7_ERROR_INVARS_SIZE,
    AB7_ERROR_SYM_NOT_PRIMITIVE,
    AB7_ERROR_SYM_BRAVAIS_XRED
  } Ab7Error;

char* ab7_error_string_from_id(Ab7Error errno);

void ab7_mpi_init();

#ifndef GLIB_MAJOR_VERSION
typedef int gboolean;
#define g_malloc(A) malloc(A)
#define g_free(A) free(A)
#endif

#endif
