// TODO: moved to effpot_xml.c, remove this file.

#include "list.h"
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <stdio.h>
#include <string.h>

void string2Array(char *input, double **farray, size_t *size) {
  char *token = strtok(input, " \n\t");
  Array tmp;
  initArray(&tmp, 3);
  while (token != NULL) {
    insertArray(&tmp, strtod(token, NULL));
    token = strtok(NULL, " ,\t\n");
  }
  copyArraytoCArray(&tmp, farray, size);
  freeArray(&tmp);
}

void string2IntArray(char *input, int **farray, size_t *size) {
  char *token = strtok(input, " \n\t");
  IntArray tmp;
  initIntArray(&tmp, 3);
  while (token != NULL) {
    insertIntArray(&tmp, strtod(token, NULL));
    token = strtok(NULL, " ,\n\t");
  }
  copyIntArrayToCIntArray(&tmp, farray, size);
  freeIntArray(&tmp);
}

int xml_read_spin_system(char *fname, double *ref_energy, double *unitcell[],
                         int *natoms, double *masses[], int *nmatoms,
                         int *index_spin[], double *gyroratios[],
                         double *damping_factors[],
                         double *positions[], double *spinat[]) {
  Array mass_array;
  Array gyroratio_array;
  Array damping_factor_array;
  IntArray index_spin_array;
  Array position_array;
  Array spinat_array;
  initArray(&mass_array, 3);
  initArray(&gyroratio_array, 3);
  initArray(&damping_factor_array, 3);
  initIntArray(&index_spin_array, 3);
  initArray(&position_array, 3);
  initArray(&spinat_array, 3);

  *natoms = 0;
  *nmatoms = 0;

  size_t size;

  xmlDocPtr doc;
  xmlNodePtr cur, cur2;
  xmlChar *key;
  doc = xmlParseFile(fname);
  if (doc == NULL) {
    fprintf(stderr, "Document parse failed. \n");
    return 1;
  }

  cur = xmlDocGetRootElement(doc);

  if (xmlStrcmp(cur->name, (const xmlChar *)("System_definition"))) {
    fprintf(stderr, "System_definition not found at the root.\n");
    return 1;
  }
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    // read energy;
    if (!xmlStrcmp(cur->name, (const xmlChar *)("energy"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      *ref_energy = strtod((const char *)key, NULL);
      printf("  energy: %lf\n", *ref_energy);
      xmlFree(key);
    }

    // read unit cell
    if (!xmlStrcmp(cur->name, (const xmlChar *)("unitcell"))) {
      printf("%s\n", cur->name);
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      printf("unitcell: %s\n", key);
      string2Array((char *)key, unitcell, &size);
      xmlFree(key);
      key = xmlGetProp(cur, BAD_CAST "units");
      printf("units: %s\n", key);
      xmlFree(key);
    }

    // read atoms
    if (!xmlStrcmp(cur->name, BAD_CAST "atom")) {
      (*natoms)++;
      printf("%s\n", cur->name);
      int ind_spin = -1;
      // mass
      key = xmlGetProp(cur, BAD_CAST "mass");

      printf("mass: %s\n", key);
      insertArray(&mass_array, strtod((const char *)key, NULL));
      xmlFree(key);


      // index_spin (optional if it's not a spin xml file. )

      key = xmlGetProp(cur, BAD_CAST "index_spin");
      if (key != NULL) {
        ind_spin = (int)strtol((const char *)key, NULL, 10);
      } else {
        ind_spin = -1;
      }

      printf("index_spin: %d\n", ind_spin);

      // gyroratio (optional)
      key = xmlGetProp(cur, BAD_CAST "gyroratio");
      printf("gyroratio: %s\n", key);
      if (key != NULL) {
        insertArray(&gyroratio_array, strtod((const char *)key, NULL));
      }
      else if(ind_spin>0) {
        insertArray(&gyroratio_array, 0.0);
      }
      xmlFree(key);

      // damping_factors(optional)
      key = xmlGetProp(cur, BAD_CAST "damping_factor");
      printf("damping_factor: %s\n", key);
      if (key != NULL) {
        insertArray(&damping_factor_array, strtod((const char *)key, NULL));
      }
      else if(ind_spin>0) {
        insertArray(&damping_factor_array, 1.0);
      }
      xmlFree(key);

      if (ind_spin > 0) {
        (*nmatoms)++;
      }
      insertIntArray(&index_spin_array, ind_spin);
      xmlFree(key);

      cur2 = cur->xmlChildrenNode;

      while (cur2 != NULL) {
        // positions
        if (!xmlStrcmp(cur2->name, BAD_CAST "position")) {
          key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
          printf("position: %s\n", key);
          double *pos;
          string2Array((char *)key, &pos, &size);
          for (size_t i = 0; i < size; ++i) {
            insertArray(&position_array, pos[i]);
          }
          xmlFree(key);
        }

        // spinat, which is optional.
        if (!xmlStrcmp(cur2->name, BAD_CAST "spinat")) {
          key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
          printf("spinat: %s\n", key);
          double *spinat_tmp;
          string2Array((char *)key, &spinat_tmp, &size);
	  for (int i=0;i<size;++i){printf("%lf ", spinat_tmp[i]);}
          if (size != 3) {
            fprintf(stderr,
                    "Error reading xml file, spinat should be a 3-vector, size is %zu", size);
          }
          for (size_t i = 0; i < size; ++i) {
            insertArray(&spinat_array, spinat_tmp[i]);
          }
          xmlFree(key);
        }

        cur2 = cur2->next;
      }

    }
    cur = cur->next;
  }
  copyArraytoCArray(&mass_array, masses, &size);
  if ((int)size != *natoms) {
    fprintf(stderr, "Number of masses not equal to number of atoms.\n");
  }
  copyIntArrayToCIntArray(&index_spin_array, index_spin, &size);
  if ((int)size != *natoms) {
    fprintf(stderr, "Number of spin_indexes not equal to number of atoms.\n");
  }
  copyArraytoCArray(&position_array, positions, &size);
  if ((int)size / 3 != *natoms) {
    fprintf(stderr, "Number of positions not equal to number of atoms.\n");
  }
  copyArraytoCArray(&spinat_array, spinat, &size);
  if ((int)size / 3 != *nmatoms) {
    fprintf(stderr, "Number of spinat not equal to number of magnetic atoms.\n");
  }
  copyArraytoCArray(&gyroratio_array, gyroratios, &size);
  if ((int)size != *nmatoms) {
    fprintf(stderr,
            "Number of gyroratios not equal to number of magnetic atoms");
  }

  copyArraytoCArray(&damping_factor_array, damping_factors, &size);
  if ((int)size != *nmatoms) {
    fprintf(stderr,
            "Number of damping_factors not equal to number of magnetic atoms");
  }
  //
  //
  printf("Here\n");
  freeIntArray(&index_spin_array);
  freeArray(&mass_array);
  freeArray(&gyroratio_array);
  freeArray(&damping_factor_array);
  freeArray(&position_array);
  freeArray(&spinat_array);
  printf("Here\n");
  return 0;
}

int xml_read_spin_exchange( char * fname, int *exc_nnz, int *exc_ilist[],
                            int *exc_jlist[], int *exc_Rlist[],
                             double *exc_vallist[]){
  *exc_nnz=0;
  IntArray i_array, j_array, R_array;
  Array val_array;
  initIntArray(&i_array, 3);
  initIntArray(&j_array, 3);
  initIntArray(&R_array, 9);
  initArray(&val_array, 9);
  int counter =0;

  xmlDocPtr doc;
  xmlNodePtr cur, cur2, cur3;
  xmlChar *key;
  doc = xmlParseFile(fname);
  if (doc == NULL) {
    fprintf(stderr, "Document parse failed. \n");
    return 1;
  }

  cur = xmlDocGetRootElement(doc);

  if (xmlStrcmp(cur->name, BAD_CAST"System_definition")) {
    fprintf(stderr, "System_definition not found at the root.\n");
    return 1;
  }
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if (!xmlStrcmp(cur->name, BAD_CAST"spin_exchange_list")) {
      cur2=cur->xmlChildrenNode;
      while (cur2 != NULL){
        if (!xmlStrcmp(cur2->name, BAD_CAST"nterms")) {
          key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
          *exc_nnz=strtol((char *)key, NULL, 10);
          xmlFree(key);
          *exc_ilist=(int *)malloc(sizeof(int)*(*exc_nnz));
          *exc_jlist=(int *)malloc(sizeof(int)*(*exc_nnz));
          *exc_Rlist=(int *)malloc(sizeof(int)*(*exc_nnz)*3);
          *exc_vallist = (double *)malloc(sizeof(double)*(*exc_nnz)*3);
        }
        if(!xmlStrcmp(cur2->name, BAD_CAST"spin_exchange_term")) {
          cur3=cur2->xmlChildrenNode;
          while(cur3!=NULL){
            if (!xmlStrcmp(cur3->name, BAD_CAST"ijR")) {
              key = xmlNodeListGetString(doc, cur3->xmlChildrenNode, 1);
              int *itmp;
              size_t size;
              string2IntArray((char *)key, &itmp, &size);
              (*exc_ilist)[counter]=itmp[0];
              (*exc_jlist)[counter]=itmp[1];
              (*exc_Rlist)[counter*3]=itmp[2];
              (*exc_Rlist)[counter*3+1]=itmp[3];
              (*exc_Rlist)[counter*3+2]=itmp[4];
              xmlFree(key);
            }
            if (!xmlStrcmp(cur3->name, BAD_CAST"data")) {
              key = xmlNodeListGetString(doc, cur3->xmlChildrenNode, 1);
              double *dtmp;
              size_t size;
              string2Array((char *)key, &dtmp, &size);
              xmlFree(key);
              for(size_t i=0; i< size; i++)
                {
                  (*exc_vallist)[counter*3+i]=dtmp[i];
                }
            }

            cur3=cur3->next;
          }
          counter++;
        }
        cur2=cur2->next;
      }
    }
    cur=cur->next;
  }
  return 0;
}


int xml_read_spin_dmi( char * fname, int *dmi_nnz, int *dmi_ilist[],
                            int *dmi_jlist[], int *dmi_Rlist[],
                             double *dmi_vallist[]){
  *dmi_nnz=0;
  IntArray i_array, j_array, R_array;
  Array val_array;
  initIntArray(&i_array, 3);
  initIntArray(&j_array, 3);
  initIntArray(&R_array, 9);
  initArray(&val_array, 9);
  int counter =0;

  xmlDocPtr doc;
  xmlNodePtr cur, cur2, cur3;
  xmlChar *key;
  doc = xmlParseFile(fname);
  if (doc == NULL) {
    fprintf(stderr, "Document parse failed. \n");
    return 1;
  }

  cur = xmlDocGetRootElement(doc);

  if (xmlStrcmp(cur->name, BAD_CAST"System_definition")) {
    fprintf(stderr, "System_definition not found at the root.\n");
    return 1;
  }
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if (!xmlStrcmp(cur->name, BAD_CAST"spin_DMI_list")) {
      cur2=cur->xmlChildrenNode;
      while (cur2 != NULL){
        if (!xmlStrcmp(cur2->name, BAD_CAST"nterms")) {
          key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
          *dmi_nnz=strtol((char *)key, NULL, 10);
          xmlFree(key);
          *dmi_ilist=(int *)malloc(sizeof(int)*(*dmi_nnz));
          *dmi_jlist=(int *)malloc(sizeof(int)*(*dmi_nnz));
          *dmi_Rlist=(int *)malloc(sizeof(int)*(*dmi_nnz)*3);
          *dmi_vallist = (double *)malloc(sizeof(double)*(*dmi_nnz)*3);
        }
        if(!xmlStrcmp(cur2->name, BAD_CAST"spin_DMI_term")) {
          cur3=cur2->xmlChildrenNode;
          while(cur3!=NULL){
            if (!xmlStrcmp(cur3->name, BAD_CAST"ijR")) {
              key = xmlNodeListGetString(doc, cur3->xmlChildrenNode, 1);
              int *itmp;
              size_t size;
              string2IntArray((char *)key, &itmp, &size);
              (*dmi_ilist)[counter]=itmp[0];
              (*dmi_jlist)[counter]=itmp[1];
              (*dmi_Rlist)[counter*3]=itmp[2];
              (*dmi_Rlist)[counter*3+1]=itmp[3];
              (*dmi_Rlist)[counter*3+2]=itmp[4];
              xmlFree(key);
            }
            if (!xmlStrcmp(cur3->name, BAD_CAST"data")) {
              key = xmlNodeListGetString(doc, cur3->xmlChildrenNode, 1);
              double *dtmp;
              size_t size;
              string2Array((char *)key, &dtmp, &size);
              xmlFree(key);
              for(size_t i=0; i< size; i++)
                {
                  (*dmi_vallist)[counter*3+i]=dtmp[i];
                }
            }

            cur3=cur3->next;
          }
          counter++;
        }
        cur2=cur2->next;
      }
    }
    cur=cur->next;
  }
  return 0;
}


int xml_read_spin_uni(char * fname, int *uni_nnz, int *uni_ilist[],
                       double *uni_amplitude_list[],
                       double *uni_direction_list[]){
  *uni_nnz=0;
  IntArray i_array;
  Array amp_array;
  Array direction_array;
  initIntArray(&i_array, 3);
  initArray(&amp_array, 3);
  initArray(&direction_array, 3);
  int counter =0;

  xmlDocPtr doc;
  xmlNodePtr cur, cur2, cur3;
  xmlChar *key;
  doc = xmlParseFile(fname);
  if (doc == NULL) {
    fprintf(stderr, "Document parse failed. \n");
    return 1;
  }

  cur = xmlDocGetRootElement(doc);

  if (xmlStrcmp(cur->name, BAD_CAST"System_definition")) {
    fprintf(stderr, "System_definition not found at the root.\n");
    return 1;
  }
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if (!xmlStrcmp(cur->name, BAD_CAST"spin_uniaxial_SIA_list")) {
      cur2=cur->xmlChildrenNode;
      while (cur2 != NULL){
        if (!xmlStrcmp(cur2->name, BAD_CAST"nterms")) {
          key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
          *uni_nnz=strtol((char *)key, NULL, 10);
          xmlFree(key);
          *uni_ilist=(int *)malloc(sizeof(int)*(*uni_nnz));
          *uni_amplitude_list = (double *)malloc(sizeof(double)*(*uni_nnz));
          *uni_direction_list = (double *)malloc(sizeof(double)*(*uni_nnz)*3);
        }
        if(!xmlStrcmp(cur2->name, BAD_CAST"spin_uniaxial_SIA_term")) {
          cur3=cur2->xmlChildrenNode;
          while(cur3!=NULL){
            if (!xmlStrcmp(cur3->name, BAD_CAST"i")) {
              key = xmlNodeListGetString(doc, cur3->xmlChildrenNode, 1);
              int *itmp;
              size_t size;
              string2IntArray((char *)key, &itmp, &size);
              (*uni_ilist)[counter]=itmp[0];
              xmlFree(key);
            }
            if (!xmlStrcmp(cur3->name, BAD_CAST"amplitude")) {
              key = xmlNodeListGetString(doc, cur3->xmlChildrenNode, 1);
              double *dtmp;
              size_t size;
              string2Array((char *)key, &dtmp, &size);
              xmlFree(key);
              for(size_t i=0; i< size; i++)
                {
                  (*uni_amplitude_list)[i]=dtmp[i];
                }
            }
            if (!xmlStrcmp(cur3->name, BAD_CAST"direction")) {
              key = xmlNodeListGetString(doc, cur3->xmlChildrenNode, 1);
              double *dtmp;
              size_t size;
              string2Array((char *)key, &dtmp, &size);
              xmlFree(key);
              for(size_t i=0; i< size; i++)
                {
                  (*uni_direction_list)[counter*3+i]=dtmp[i];
                }
            }

            cur3=cur3->next;
          }
          counter++;
        }
        cur2=cur2->next;
      }
    }
    cur=cur->next;
  }
  return 0;
}



int xml_read_spin_bilinear( char * fname, int *bi_nnz, int *bi_ilist[],
                            int *bi_jlist[], int *bi_Rlist[],
                             double *bi_vallist[]){
  *bi_nnz=0;
  IntArray i_array, j_array, R_array;
  Array val_array;
  initIntArray(&i_array, 3);
  initIntArray(&j_array, 3);
  initIntArray(&R_array, 9);
  initArray(&val_array, 27);
  int counter =0;

  xmlDocPtr doc;
  xmlNodePtr cur, cur2, cur3;
  xmlChar *key;
  doc = xmlParseFile(fname);
  if (doc == NULL) {
    fprintf(stderr, "Document parse failed. \n");
    return 1;
  }

  cur = xmlDocGetRootElement(doc);

  if (xmlStrcmp(cur->name, BAD_CAST"System_definition")) {
    fprintf(stderr, "System_definition not found at the root.\n");
    return 1;
  }
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if (!xmlStrcmp(cur->name, BAD_CAST"spin_bilinear_list")) {
      cur2=cur->xmlChildrenNode;
      while (cur2 != NULL){
        if (!xmlStrcmp(cur2->name, BAD_CAST"nterms")) {
          key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
          *bi_nnz=strtol((char *)key, NULL, 10);
          xmlFree(key);
          *bi_ilist=(int *)malloc(sizeof(int)*(*bi_nnz));
          *bi_jlist=(int *)malloc(sizeof(int)*(*bi_nnz));
          *bi_Rlist=(int *)malloc(sizeof(int)*(*bi_nnz)*3);
          *bi_vallist = (double *)malloc(sizeof(double)*(*bi_nnz)*9);
        }
        if(!xmlStrcmp(cur2->name, BAD_CAST"spin_bilinear_term")) {
          cur3=cur2->xmlChildrenNode;
          while(cur3!=NULL){
            if (!xmlStrcmp(cur3->name, BAD_CAST"ijR")) {
              key = xmlNodeListGetString(doc, cur3->xmlChildrenNode, 1);
              int *itmp;
              size_t size;
              string2IntArray((char *)key, &itmp, &size);
              (*bi_ilist)[counter]=itmp[0];
              (*bi_jlist)[counter]=itmp[1];
              (*bi_Rlist)[counter*3]=itmp[2];
              (*bi_Rlist)[counter*3+1]=itmp[3];
              (*bi_Rlist)[counter*3+2]=itmp[4];
              xmlFree(key);
            }
            if (!xmlStrcmp(cur3->name, BAD_CAST"data")) {
              key = xmlNodeListGetString(doc, cur3->xmlChildrenNode, 1);
              double *dtmp;
              size_t size;
              string2Array((char *)key, &dtmp, &size);
              xmlFree(key);
              for(size_t i=0; i< size; i++)
                {
                  (*bi_vallist)[counter*9+i]=dtmp[i];
                }
            }

            cur3=cur3->next;
          }
          counter++;
        }
        cur2=cur2->next;
      }
    }
    cur=cur->next;
  }
  return 0;
}

void xml_read_spin(char *fname, double *ref_energy, double *unitcell[9],
                   int *natoms, double *masses[], int *nmatoms,
                   int *index_spin[], double *gyroratios[], double *damping_factors[],
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
                   double *bi_vallist[]){
  printf("reading xml: system.\n");
  xml_read_spin_system(fname, ref_energy, unitcell, natoms, masses, nmatoms, index_spin, gyroratios, damping_factors, positions, spinat);

  printf("reading xml: exchange.\n");
  xml_read_spin_exchange(fname, exc_nnz, exc_ilist, exc_jlist, exc_Rlist, exc_vallist);
  printf(" %d terms readed\n", *exc_nnz);

  printf("reading xml: DMI.\n");
  xml_read_spin_dmi(fname, dmi_nnz, dmi_ilist, dmi_jlist, dmi_Rlist, dmi_vallist);
  printf(" %d terms readed.\n", *dmi_nnz);

  printf("reading xml: uniaxial single ion anisotropy.\n");
  xml_read_spin_uni(fname, uni_nnz, uni_ilist, uni_amplitude_list, uni_direction_list);
  printf(" %d terms readed.\n", *uni_nnz);

  printf("reading xml: bilinear.\n");
  xml_read_spin_bilinear(fname, bi_nnz, bi_ilist, bi_jlist, bi_Rlist, bi_vallist);
  printf(" %d terms readed.\n", *bi_nnz);

  printf("Reading xml finished!");
}


int test_read_xml() {

  char * fname="test_f.xml";
  double ref_energy;
  double *unitcell, *masses, *gyroratios, *damping_factors, *positions, *spinat;
  int natoms, nmatoms, *index_spin;
  xml_read_spin_system("test_f.xml", &ref_energy, &unitcell, &natoms, &masses,
   &nmatoms, &index_spin, &gyroratios, &damping_factors, &positions, &spinat);

  // exchange
  printf("======Exchange Terms========\n");
  int exc_nnz, *exc_ilist, *exc_jlist, *exc_Rlist;
  double *exc_vallist;
  xml_read_spin_exchange(fname, &exc_nnz, &exc_ilist, &exc_jlist, &exc_Rlist, &exc_vallist);
  for(int i=0;i<exc_nnz;i++){
    printf("%d\t%d\t%d\t%d\t%d\t : %E\t%E\t%E\n", exc_ilist[i], exc_jlist[i], exc_Rlist[3*i], exc_Rlist[3*i+1], exc_Rlist[3*i+2], exc_vallist[3*i], exc_vallist[3*i+1], exc_vallist[3*i+2]);
  }

  //dmi
  printf("======DMI Terms========\n");
  int dmi_nnz, *dmi_ilist, *dmi_jlist, *dmi_Rlist;
  double *dmi_vallist;
  xml_read_spin_dmi(fname, &dmi_nnz, &dmi_ilist, &dmi_jlist, &dmi_Rlist, &dmi_vallist);
  for(int i=0;i<dmi_nnz;i++){
    printf("%d\t%d\t%d\t%d\t%d\t :", dmi_ilist[i], dmi_jlist[i], dmi_Rlist[3*i], dmi_Rlist[3*i+1], dmi_Rlist[3*i+2]);
    for (int j=0; j < 3; ++j) {
      printf("%E\t", dmi_vallist[i*3+j]);
    }
    printf("\n");
  }

  //uniaxial SIA
  printf("======uniaxial SIA Terms========\n");
  int uni_nnz, *uni_ilist;
  double *uni_amplitude_list, *uni_direction_list;
  xml_read_spin_uni(fname, &uni_nnz, &uni_ilist, &uni_amplitude_list, &uni_direction_list);
  for(int i=0;i<uni_nnz;i++){
    printf("%d\t :", uni_ilist[i]);
    printf("%E\t:", uni_amplitude_list[i]);
    for (int j=0; j < 3; ++j) {
      printf("%E\t", uni_direction_list[i*3+j]);
    }
    printf("\n");
  }


  //bilinear
  printf("======Bilinear Terms========\n");
  int bi_nnz, *bi_ilist, *bi_jlist, *bi_Rlist;
  double *bi_vallist;
  xml_read_spin_bilinear("test_f.xml", &bi_nnz, &bi_ilist, &bi_jlist, &bi_Rlist, &bi_vallist);
  for(int i=0;i<bi_nnz;i++){
    printf("%d\t%d\t%d\t%d\t%d\t :", bi_ilist[i], bi_jlist[i], bi_Rlist[3*i], bi_Rlist[3*i+1], bi_Rlist[3*i+2]);
    for (int j=0; j < 9; ++j) {
      printf("%E\t", bi_vallist[i*9+j]);
    }
    printf("\n");
  }


  return 0;
}

//int main(){
//  test_read_xml();
//}
