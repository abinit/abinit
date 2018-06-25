/** -*- C source file -*-
 *
 * Copyright (C) 2009-2018 ABINIT Group (Damien Caliste)
 *
 * This file is part of the ABINIT software package. For license information,
 * please see the COPYING file in the top-level directory of the ABINIT source
 * distribution.
 *
 */

#include "ab7_base.h"
#include <config.h>

/*
  AB7_NO_ERROR,
  AB7_ERROR_OBJ,
  AB7_ERROR_ARG,
  AB7_ERROR_INVARS_ATT,
  AB7_ERROR_INVARS_ID,
  AB7_ERROR_INVARS_SIZE,
  AB7_ERROR_SYM_NOT_PRIMITIVE,
  AB7_ERROR_SYM_BRAVAIS_XRED
*/
static char* messages[8] = {
  "No error.",
  "Wrong pointer to object.",
  "Wrong value for one argument of the calling routine.",
  "Unknown attribute or wrong attribute type from Dtset structure.",
  "Out of bounds dtset number.",
  "Wrong input size for array.",
  "The cell is not primitive.",
  "The bravais lattice has more symmetries than the lattice"
  " system found from the atom coordinates."};

char* ab7_error_string_from_id(Ab7Error errno)
{
  if (errno < 0 || errno >= 8)
    return "Unknown error id.";
  else
    return messages[errno];
}

void ABI_FC_MOD(m_xmpi,M_XMPI, xmpi_init,XMPI_INIT)();
void ab7_mpi_init()
{
  ABI_FC_MOD(m_xmpi,M_XMPI, xmpi_init,XMPI_INIT)();
}
