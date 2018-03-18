// TODO hexu: should this file be removed?
// It is not used by the fortran modules, but in C convention, it's better to keep a header file of the functions.

#include <stdio.h>
#include <string.h>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include "list.h" // TODO hexu: Note that list.h might be removed.



int xml_read_spin_system(char *fname, double *ref_energy, double *unitcell[9],
                         int *natoms, double *masses[], int *nmatoms,
                         int *index_spin[], double *gyroratios[],
                         double *damping_factors[],
                         double *positions[], double *spinat[]);

int xml_read_spin_exchange( char * fname, int *exc_nnz, int *exc_ilist[],
                            int *exc_jlist[], int *exc_Rlist[],
                            double *exc_vallist[]);

int xml_read_spin_dmi( char * fname, int *dmi_nnz, int *dmi_ilist[],
                       int *dmi_jlist[], int *dmi_Rlist[],
                       double *dmi_vallist[]);

int xml_read_spin_uni(char * fname, int *uni_nnz, int *uni_ilist[],
                        double *uni_amplitude_list[],
                        double *uni_direction_list[]);

int xml_read_spin_bilinear( char * fname, int *bi_nnz, int *bi_ilist[],
                            int *bi_jlist[], int *bi_Rlist[],
                            double *bi_vallist[]);

void xml_read_spin(char *fname, double *ref_energy, double *unitcell[9],
                   int *natoms, double *masses[], int *nmatoms,
                   int *index_spin[], double *gyroratios[],
                   double *positions[], double *spinat[],
                   // exchange
                   int *exc_nnz, int *exc_ilist[],
                   int *exc_jlist[], int *exc_Rlist[],
                   double *exc_vallist[],
                   //dmi
                   int *dmi_nnz, int *dmi_ilist[],
                   int *dmi_jlist[], int *dmi_Rlist[],
                   double *dmi_vallist[],
                   //uniaxial SIA
                   int *uni_nnz, int *uni_ilist[],
                   double *uni_amplitude_list[],
                   double *uni_direction_list[],
                   //bilinear
                   int *bi_nnz, int *bi_ilist[],
                   int *bi_jlist[], int *bi_Rlist[],
                   double *bi_vallist[]);
