/*
 * Copyright (C) 2015-2019 ABINIT group (AM)
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
 */

/* ===============================================================
 * Set of C functions interfacing the LibXML library.
 * ===============================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <config.h>
#include <sys/stat.h>

#if defined HAVE_XML

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>

static double const eV=1.6e-19;

//define type for dynamical double format array 
typedef struct {
  double *array;
  size_t used;
  size_t size;
} Array;

void initArray(Array *a, size_t initialSize) {
  a->array = (double *)malloc(initialSize * sizeof(double));
  a->used = 0;
  a->size = initialSize;
}

void insertArray(Array *a, double element) {
  if (a->used == a->size) {
    a->size *= 2;
    a->array = (double *)realloc(a->array, a->size * sizeof(double));
  }
  a->array[a->used++] = element;
}

void freeArray(Array *a) {
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
}

void copyArraytoCArray(Array *l, double **a, size_t* size){
    size_t i;
  *size=0;
  *a=(double *) malloc(sizeof(double)*l->used);
  for(i=0;i<l->used;i++){
    (*a)[i]=l->array[i];
    (*size)++;
  }
}

//define type for dynamical int type array
typedef struct {
  int *array;
  size_t used;
  size_t size;
} IntArray;

void initIntArray(IntArray *a, size_t initialSize) {
  a->array = (int *)malloc(initialSize * sizeof(int));
  a->used = 0;
  a->size = initialSize;
}

void insertIntArray(IntArray *a, int element) {
  if (a->used == a->size) {
    a->size *= 2;
    a->array = (int *)realloc(a->array, a->size * sizeof(int));
  }
  a->array[a->used++] = element;
}

void freeIntArray(IntArray *a) {
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
}

void copyIntArrayToCIntArray(IntArray *l, int **a, size_t *size){
    size_t i;
  *size=0;
  *a=(int *) malloc(sizeof(int)*l->used);
  for(i=0;i<l->used;i++){
    (*a)[i]=l->array[i];
    (*size)++;
  }
}


// string to array of double.
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

// string to array of int
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

// check if file exist
int file_exists(const char* filename){
    struct stat buffer;
    int exist = stat(filename,&buffer);
    if(exist == 0)
        return 1;
    else 
        return 0;
}


void effpot_xml_checkXML(char *filename,char *name_xml){  
  xmlDocPtr doc;
  xmlNodePtr cur;

  doc = xmlParseFile(filename);
  if (doc == NULL) printf(" error: could not parse file file.xml\n");

  cur = xmlDocGetRootElement(doc);
  printf(" Root node xml : %s\n",cur -> name);
  if (cur == NULL) {
    fprintf(stderr," The document is empty \n");
    xmlFreeDoc(doc);
    return;
  }

  if (xmlStrcmp(cur->name, (const  xmlChar *) name_xml)) {
    fprintf(stderr," Document of the wrong type, root node != %s\n",name_xml);
    xmlFreeDoc(doc);
     return;
  }
  xmlFreeDoc(doc);
}

void effpot_xml_getDimSystem(char *filename,int *natom,int *ntypat, int *nqpt, int *loc_nrpt,\
                             int *tot_nrpt){
  xmlDocPtr doc;
  int i,iatom,irpt1,irpt2,iqpt,itypat,j,present;
  xmlNodePtr cur,cur2;
  xmlChar *key,*uri;
  Array typat;

  initArray(&typat, 1);
  *natom  = 0;  *nqpt  = 0; *loc_nrpt = 0; *tot_nrpt = 0;  *ntypat = 0;
  iatom   = 0; irpt1   = 0; irpt2  = 0;  iqpt  = 0;  itypat = 0;
  present = 0;
  typat.array[0] = 0;

  doc = xmlParseFile(filename);
  if (doc == NULL) printf(" error: could not parse file file.xml\n");

  cur = xmlDocGetRootElement(doc);
  if (cur == NULL) {
    fprintf(stderr," The document is empty \n");
    xmlFreeDoc(doc);
    return;
  }
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const  xmlChar *) "atom"))) {
      iatom++;
      uri = xmlGetProp(cur, (const  xmlChar *) "mass");
      insertArray(&typat,strtod(uri,NULL)); 
      xmlFree(uri);
    } 
    if ((!xmlStrcmp(cur->name, (const  xmlChar *) "local_force_constant"))) {irpt1++;}
    if ((!xmlStrcmp(cur->name, (const  xmlChar *) "total_force_constant"))) {irpt2++;}
    if ((!xmlStrcmp(cur->name, (const  xmlChar *) "phonon"))) {
      cur2 = cur->xmlChildrenNode;
      while (cur2 != NULL) {
        if ((!xmlStrcmp(cur2->name, (const  xmlChar *) "qpoint"))) {iqpt++;}
        cur2 = cur2->next;
      }
    }
    cur = cur->next;
  }
  for(i=0;i<typat.used;i++){
    present = 0;
    for(j=i+1;j<typat.used;j++){
      if (typat.array[i] == typat.array[j] && typat.array[i] != 0){
        present = 1;
        break;
      }
    }
    if(present==0 && typat.array[i] !=0){
      itypat++;
    }
  }
  xmlFreeDoc(doc);
  freeArray(&typat);

  *natom  = iatom;
  *nqpt   = iqpt;
  *loc_nrpt   = irpt1;
  *tot_nrpt   = irpt2;
  *ntypat = itypat;
}

void effpot_xml_readSystem(char *filename,int *natom,int *ntypat,int *nrpt,int *nqpt,
                           double amu[*ntypat],double atmfrc[*nrpt][*natom][3][*natom][3],
                           int cell[*nrpt][3],double dynmat[*nqpt][*natom][3][*natom][3][2],
                           double elastic_constants[6][6],
                           double *energy,double epsilon_inf[3][3],
                           double ewald_atmfrc[*nrpt][*natom][3][*natom][3],
                           double phfrq[*nqpt][3* *natom],
                           double rprimd[3][3],double qph1l[*nqpt][3],
                           double short_atmfrc[*nrpt][*natom][3][*natom][3],
                           int typat[*natom],double xcart[*natom][3],double zeff[*natom][3][3]){
  xmlDocPtr doc;
  char *pch;
  double total_atmfrc[*nrpt][*natom][3][*natom][3];
  double local_atmfrc[*nrpt][*natom][3][*natom][3];
  int cell_local[*nrpt][3];
  int cell_total[*nrpt][3];
  int iatom,iamu,irpt1,irpt2,irpt3,iqpt,present;
  int ia,ib,mu,nu,voigt;
  int i,j;
  xmlNodePtr cur,cur2;
  xmlChar *key,*uri;

  if (*natom <= 0){ 
    printf(" error: The number of atom must be superior to zero\n");
    exit(0);
  }

  iatom   = 0;
  iamu    = 0;
  present = 0;
  irpt1 = 0;
  irpt2 = 0;
  irpt3 = 0;
  iqpt  = 0;
  voigt = 0;

  doc = xmlParseFile(filename);
  if (doc == NULL) printf(" error: could not parse file file.xml\n");

  cur = xmlDocGetRootElement(doc);
  if (cur == NULL) {
    fprintf(stderr," The document is empty \n");
    xmlFreeDoc(doc);
    return;
  }

  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const  xmlChar *) "energy"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      *energy = strtod(key, NULL);
      xmlFree(key);
    }
    else if ((!xmlStrcmp(cur->name, (const  xmlChar *) "unit_cell"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      pch = strtok(key,"\t \n");
      for(mu=0;mu<3;mu++){
        for(nu=0;nu<3;nu++){
          if (pch != NULL){
            rprimd[nu][mu]=strtod(pch,NULL);
            pch = strtok(NULL,"\t \n");
          }
        }
      }
      xmlFree(key);
    }
    else if ((!xmlStrcmp(cur->name, (const  xmlChar *) "epsilon_inf"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      pch = strtok(key,"\t \n");
      for(mu=0;mu<3;mu++){
        for(nu=0;nu<3;nu++){
          if (pch != NULL){
            epsilon_inf[mu][nu]=strtod(pch,NULL);
            pch = strtok(NULL,"\t \n");
          }
        }
      }
      xmlFree(key);
    }
    else if ((!xmlStrcmp(cur->name, (const  xmlChar *) "elastic"))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      pch = strtok(key,"\t \n");
      for(mu=0;mu<6;mu++){
        for(nu=0;nu<6;nu++){
          if (pch != NULL){
            elastic_constants[mu][nu]=strtod(pch,NULL);
            pch = strtok(NULL,"\t \n");
          }
        }
      }
      xmlFree(key);
    }
    else if ((!xmlStrcmp(cur->name, (const  xmlChar *) "atom"))) {
      uri = xmlGetProp(cur, (const  xmlChar *) "mass");
      present = 0;
      //1) fill the atomic mass unit
      for(i=0;i<*ntypat;i++){
        if(abs(amu[i]-strtod(uri,NULL))<1e-5){
        //if(amu[i]==strtod(uri,NULL)){
          present = 1;
          break;
        }
      }
      if(present==0){
        amu[iamu]=strtod(uri,NULL);
        iamu++;
      }
      // fill the typat table
      //printf("=====typat====\n");
      //printf("ntypat: %d\n", *ntypat);
      for(i=0;i<*ntypat;i++){
        if(abs(amu[i]-strtod(uri,NULL))<1e-5){
          typat[iatom]=i+1;
      // Debug typat
      //   printf("i= %d\n", i);
      //   printf("%d: %f, %f, %d\n", iatom, amu[i], strtod(uri,NULL), typat[iatom]);
        }
      }
      xmlFree(uri);
      
      cur2 = cur->xmlChildrenNode;
      while (cur2 != NULL) {
        if (iatom<=*natom) {
          if ((!xmlStrcmp(cur2->name, (const  xmlChar *) "position"))) {
            key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
            pch = strtok(key,"\t \n");
            for(mu=0;mu<3;mu++){
              if (pch != NULL){
                xcart[iatom][mu]=strtod(pch,NULL);
                pch = strtok(NULL,"\t \n");
              }
            }
          }
          if ((!xmlStrcmp(cur2->name, (const  xmlChar *) "borncharge"))) {
            key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
            pch = strtok(key,"\t \n");
            for(mu=0;mu<3;mu++){
              for(nu=0;nu<3;nu++){
                if (pch != NULL){
                  zeff[iatom][mu][nu]=strtod(pch,NULL);
                  pch = strtok(NULL,"\t \n");
                }
              }
            }
            iatom++;
          }
        }
        else{
          printf(" error: The number of atom doesn't match with the XML file\n");
          exit(0);
        } 
        cur2 = cur2->next;
      }
    }
    else if ((!xmlStrcmp(cur->name, (const  xmlChar *) "local_force_constant"))) { 
      cur2 = cur->xmlChildrenNode;
      while (cur2 != NULL) {
        if (irpt1<=*nrpt) {
          if ((!xmlStrcmp(cur2->name, (const  xmlChar *) "data"))) {
            key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
            pch = strtok(key,"\t \n");
            for(ia=0;ia<*natom;ia++){
              for(mu=0;mu<3;mu++){
                for(ib=0;ib<*natom;ib++){
                  for(nu=0;nu<3;nu++){
                    if (pch != NULL){
                      local_atmfrc[irpt1][ib][nu][ia][mu]=strtod(pch,NULL);
                      pch = strtok(NULL,"\t \n");
                    }
                  }
                }
              }
            }
          }
          if ((!xmlStrcmp(cur2->name, (const  xmlChar *) "cell"))) {
            key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
            pch = strtok(key,"\t \n");
            for(i=0;i<3;i++){
              cell_local[irpt1][i]=atoi(pch);
              pch = strtok(NULL,"\t \n");
            }
          }
        }
        else{
          printf(" error: The number of ifc doesn't match with the XML file %d %d\n",irpt1,*nrpt);
          exit(0);
        } 
        cur2 = cur2->next;
      }
      irpt1++;
    }
    else if ((!xmlStrcmp(cur->name, (const  xmlChar *) "total_force_constant"))) {      
      cur2 = cur->xmlChildrenNode;
      while (cur2 != NULL) {
        if (irpt2<=*nrpt) {
          if ((!xmlStrcmp(cur2->name, (const  xmlChar *) "data"))) {
            key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
            pch = strtok(key,"\t \n");
            for(ia=0;ia<*natom;ia++){
              for(mu=0;mu<3;mu++){
                for(ib=0;ib<*natom;ib++){
                  for(nu=0;nu<3;nu++){
                    if (pch != NULL){
                      total_atmfrc[irpt2][ib][nu][ia][mu]=strtod(pch,NULL);
                      pch = strtok(NULL,"\t \n");
                    }
                  }
                }
              }
            }
          }
          if ((!xmlStrcmp(cur2->name, (const  xmlChar *) "cell"))) {
            key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
            pch = strtok(key,"\t \n");
            for(i=0;i<3;i++){
              cell_total[irpt2][i]=atoi(pch);
              pch = strtok(NULL,"\t \n");
            }
          }
        }
        else{
          printf(" error: The number of ifc doesn't match with the XML file %d %d\n",irpt2,*nrpt);
          exit(0);
        } 
        cur2 = cur2->next;
      }
      irpt2++;
    }
    else if ((!xmlStrcmp(cur->name, (const  xmlChar *) "phonon"))){
      cur2 = cur->xmlChildrenNode;
      while (cur2 != NULL) {
        if (iqpt<=*nqpt) {
          if ((!xmlStrcmp(cur2->name, (const  xmlChar *) "qpoint"))) {
            key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
            pch = strtok(key,"\t \n");
            for(mu=0;mu<3;mu++){
              if (pch != NULL){
                qph1l[iqpt][mu]=strtod(pch,NULL);
                pch = strtok(NULL,"\t \n");
              }
            }
          }
          if ((!xmlStrcmp(cur2->name, (const  xmlChar *) "frequencies"))) {
            key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
            pch = strtok(key,"\t \n");
            for(mu=0;mu<3**natom;mu++){
              if (pch != NULL){
                phfrq[iqpt][mu]=strtod(pch,NULL);
                pch = strtok(NULL,"\t \n");
              }
            }
          }
          if ((!xmlStrcmp(cur2->name, (const  xmlChar *) "dynamical_matrix"))) {
            key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
            pch = strtok(key,"\t \n");
            for(ia=0;ia<*natom;ia++){
              for(mu=0;mu<3;mu++){
                for(ib=0;ib<*natom;ib++){
                  for(nu=0;nu<3;nu++){
                    if (pch != NULL){
                      dynmat[iqpt][ia][mu][ib][nu][0]=strtod(pch,NULL);
                      pch = strtok(NULL,"\t \n");
                    }
                  }
                }
              }
            }
          }
        }
        else{
          printf(" error: The number of qpoints doesn't match with the XML file %d %d\n",irpt2,*nrpt);
          exit(0);
        }
        cur2 = cur2->next;
      }
      iqpt++;
    }
    cur = cur->next;
  }
  xmlFreeDoc(doc);

  //  Reorder the ATMFRC
  //  Case 1: only local in the xml
  if (irpt1>0 && irpt2==0){
    for(i=0;i<irpt1;i++){
      for(j=0;j<3;j++){
        cell[i][j] = cell_local[i][j];        
      }
    }
    for(i=0;i<irpt1;i++){
      for(ia=0;ia<*natom;ia++){
        for(mu=0;mu<3;mu++){
          for(ib=0;ib<*natom;ib++){
            for(nu=0;nu<3;nu++){
              atmfrc[i][ib][nu][ia][mu]=local_atmfrc[i][ib][nu][ia][mu];
              short_atmfrc[i][ib][nu][ia][mu]=local_atmfrc[i][ib][nu][ia][mu];
              ewald_atmfrc[i][ib][nu][ia][mu]=0.0;
            }
          }
        }    
      }
    }
  //Case 2: only total in the xml
  }else if (irpt1==0 && irpt2>0){
    for(i=0;i<irpt1;i++){
      for(j=0;j<3;j++){
        cell[i][j] = cell_total[i][j];        
      }
    }
    for(i=0;i<irpt1;i++){
      for(ia=0;ia<*natom;ia++){
        for(mu=0;mu<3;mu++){
          for(ib=0;ib<*natom;ib++){
            for(nu=0;nu<3;nu++){
              atmfrc[i][ib][nu][ia][mu]=total_atmfrc[i][ib][nu][ia][mu];
              short_atmfrc[i][ib][nu][ia][mu]=0.0;
              ewald_atmfrc[i][ib][nu][ia][mu]=total_atmfrc[i][ib][nu][ia][mu];
            }
          }
        }    
      }
    }
  //Case 3: local + total in the xml
  }else if (irpt1>0 && irpt2>0){
    if (irpt1 <= irpt2){
      for(i=0;i<irpt2;i++){
        for(j=0;j<3;j++){
          cell[i][j] = cell_total[i][j];
        }
      }
      for(i=0;i<irpt2;i++){
        for(ia=0;ia<*natom;ia++){
          for(mu=0;mu<3;mu++){
            for(ib=0;ib<*natom;ib++){
              for(nu=0;nu<3;nu++){
                atmfrc[i][ib][nu][ia][mu] = total_atmfrc[i][ib][nu][ia][mu];
                ewald_atmfrc[i][ib][nu][ia][mu]= atmfrc[i][ib][nu][ia][mu];
                for(j=0;j<irpt1;j++){
                  if(cell_local[j][0] == cell[i][0] && 
                     cell_local[j][1] == cell[i][1] &&
                     cell_local[j][2] == cell[i][2] ){
                    if(ia==0 && ib==0 && mu==0 && nu==0){irpt3++;}
                    short_atmfrc[i][ib][nu][ia][mu]= local_atmfrc[j][ib][nu][ia][mu];
                    ewald_atmfrc[i][ib][nu][ia][mu] -=  short_atmfrc[i][ib][nu][ia][mu];
                  }
                }
              }    
            }
          }
        }
      }
      if(irpt3 != irpt1){
        fprintf(stdout,"\n WARNING: The number of local and total rpt are not equivalent\n");
        fprintf(stdout,"          in the XML file :%d and %d\n",irpt1,irpt3);
        fprintf(stdout,"          The missing local IFC will be set to zero\n");        
      }
    }
    else{
      fprintf(stderr,"error: Local rpt is superior to total rpt in the XML file:%d %d\n",\
              irpt1,irpt2);
      exit(0);
    }
  } 
  // TODO: hexu Temporarily disabled by hexu to make the spin only model work
  //else{
  //  fprintf(stderr,"error: Number of local and total rpt doesn't match with the XML file:%d %d\n",\
  //          irpt1,irpt2);
  //  exit(0);
  //}
}

void effpot_xml_getDimStrainCoupling(char *filename, int *nrpt,int *voigt){
  xmlDocPtr doc;
  int irpt;
  xmlNodePtr cur,cur2,cur3;
  xmlChar *uri;

  doc = xmlParseFile(filename);
  if (doc == NULL) printf(" error: could not parse file file.xml\n");

  cur = xmlDocGetRootElement(doc);
  if (cur == NULL) {
    fprintf(stderr," The document is empty \n");
    xmlFreeDoc(doc);
    return;
  }
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const  xmlChar *) "strain_coupling"))){
      irpt = 0;
      uri = xmlGetProp(cur, (const  xmlChar *) "voigt");
      if (strtod(uri,NULL) == *voigt){
        cur2 = cur->xmlChildrenNode;
        while (cur2 != NULL) {
          if (*voigt<=12) {
            if ((!xmlStrcmp(cur2->name, (const  xmlChar *) "correction_force_constant"))) {
              cur3 = cur2->xmlChildrenNode;
              while (cur3 != NULL) {
                if ((!xmlStrcmp(cur3->name, (const  xmlChar *) "data"))) {
                  irpt++;
                }
                cur3 = cur3->next;             
              }
            }
            *nrpt = irpt;
          }
          else{
            printf(" error: The number of strain doesn't match with the XML file\n");
            exit(0);
          }
          cur2 = cur2->next;
        }
      }
      xmlFree(uri);
    }
    cur = cur->next;
  }
  xmlFreeDoc(doc);  
}

void effpot_xml_readStrainCoupling(char *filename,int *natom,int *nrpt,int *voigt,
                                   double elastic3rd[6][6], double elastic_displacement[*natom][3][6],
                                   double internal_strain[*natom][3],
                                   double phonon_strain_atmfrc[*nrpt][*natom][3][*natom][3],
                                   int phonon_straincell[*nrpt][3]){
  xmlDocPtr doc;
  char *pch;
  int i,irpt,ia,ib,mu,nu;
  xmlNodePtr cur,cur2,cur3;
  xmlChar *key,*uri;

  if (*natom <= 0){ 
    printf(" error: The number of atom must be superior to zero\n");
    exit(0);
  }

  irpt    = 0;

  doc = xmlParseFile(filename);
  if (doc == NULL) printf(" error: could not parse file file.xml\n");

  cur = xmlDocGetRootElement(doc);
  if (cur == NULL) {
    fprintf(stderr," The document is empty \n");
    xmlFreeDoc(doc);
    return;
  }
 
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const  xmlChar *) "strain_coupling"))){
      irpt = 0;
      uri = xmlGetProp(cur, (const  xmlChar *) "voigt");
      if(atoi(uri) == *voigt){
        cur2 = cur->xmlChildrenNode;
        while (cur2 != NULL) {
          if (atoi(uri)<=12) {
            if ((!xmlStrcmp(cur2->name, (const  xmlChar *) "correction_force"))) {
              key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
              pch = strtok(key,"\t \n");
              for(ia=0;ia<*natom;ia++){
                for(mu=0;mu<3;mu++){
                  if (pch != NULL){
                    internal_strain[ia][mu]=strtod(pch,NULL);
                    pch = strtok(NULL,"\t \n");
                  }
                }
              }
            }
            if ((!xmlStrcmp(cur2->name, (const  xmlChar *) "elastic3rd"))) {
              key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
              pch = strtok(key,"\t \n");
              for(mu=0;mu<6;mu++){
                for(nu=0;nu<6;nu++){
                  if (pch != NULL){
                    elastic3rd[mu][nu]=strtod(pch,NULL);
                    pch = strtok(NULL,"\t \n");
                  }
                }
              }
            }
            if ((!xmlStrcmp(cur2->name, (const  xmlChar *) "correction_strain_force"))) {
              key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
              pch = strtok(key,"\t \n");
              for(ia=0;ia<*natom;ia++){
                for(mu=0;mu<3;mu++){
                  for(nu=0;nu<6;nu++){
                    if (pch != NULL){
                      elastic_displacement[ia][mu][nu]=strtod(pch,NULL);
                      pch = strtok(NULL,"\t \n");
                    }
                  }
                }
              }
            }
            if ((!xmlStrcmp(cur2->name, (const  xmlChar *) "correction_force_constant"))) {
              cur3 = cur2->xmlChildrenNode;
              while (cur3 != NULL) {
                if (irpt<=*nrpt) {
                  if ((!xmlStrcmp(cur3->name, (const  xmlChar *) "data"))) {
                    key = xmlNodeListGetString(doc, cur3->xmlChildrenNode, 1);
                    pch = strtok(key,"\t \n");
                    for(ia=0;ia<*natom;ia++){
                      for(mu=0;mu<3;mu++){
                        for(ib=0;ib<*natom;ib++){
                          for(nu=0;nu<3;nu++){
                            if (pch != NULL){
                              phonon_strain_atmfrc[irpt][ib][nu][ia][mu]=strtod(pch,NULL);
                              pch = strtok(NULL,"\t \n");
                            }
                          }
                        }
                      }
                    }
                  }
                  if ((!xmlStrcmp(cur3->name, (const  xmlChar *) "cell"))) {
                    key = xmlNodeListGetString(doc, cur3->xmlChildrenNode, 1);
                    pch = strtok(key,"\t \n");
                    for(i=0;i<3;i++){
                      phonon_straincell[irpt][i]=atoi(pch);
                      pch = strtok(NULL,"\t \n");
                    }
                    irpt++;
                  }
                }
                cur3 = cur3->next;             
              }
            }
          }
          else{
            printf(" error: The number of strain doesn't match with the XML file %d %d\n",irpt,*nrpt);
            exit(0);
          }
          cur2 = cur2->next;
        }
      }
      xmlFree(uri);
    }

    cur = cur->next;
  }
  xmlFreeDoc(doc);
}

void effpot_xml_readCoeff(char *filename,int*ncoeff,int*ndisp,int*nterm,
                          double coefficient[*ncoeff],
                          int atindx[*ndisp][2][*nterm][*ncoeff],
                          int cell[*ndisp][2][3][*nterm][*ncoeff],
                          int direction[*ndisp][*nterm][*ncoeff],
                          int power_disp[*ndisp][*nterm][*ncoeff],
                          int power_strain[*ndisp][*nterm][*ncoeff],
                          int strain[*ndisp][*nterm][*ncoeff],
                          double weight[*nterm][*ncoeff]){
  
  int i,idisp,istrain,j,iterm,icoeff;
  xmlDocPtr doc;
  char * pch;
  xmlNodePtr cur,cur2,cur3,cur4;
  xmlChar *key,*uri,*uri2;

  if (*ncoeff <= 0){ 
    printf(" error: The number of coeff must be superior to zero\n");
    exit(0);
  }

  if (*ndisp <= 0){ 
    printf(" error: The number of coeff must be superior to zero\n");
    exit(0);
  }

  if (*nterm <= 0){ 
    printf(" error: The number of coeff must be superior to zero\n");
    exit(0);
  }
 
  //Set to zero outputs
  for (icoeff=0; icoeff < *ncoeff ;icoeff++){
    coefficient[icoeff]=0;
    for (iterm=0; iterm < *nterm ;iterm++){
      weight[iterm][icoeff]=0;
      for (idisp=0; idisp < *ndisp ;idisp++){
        direction[idisp][iterm][icoeff] = 0;
        power_disp[idisp][iterm][icoeff] = 0;
        strain[idisp][iterm][icoeff] = 0;
        power_strain[idisp][iterm][icoeff] = 0;
        for (i=0;i<2;i++){
          atindx[idisp][i][iterm][icoeff] = 0;
          for (j=0;j<3;j++){
            cell[idisp][i][j][iterm][icoeff] = 0;
          }
        }      
      }
    }
  }

  idisp = 0;

  doc = xmlParseFile(filename);
  if (doc == NULL) printf(" error: could not parse file file.xml\n");
 
  cur = xmlDocGetRootElement(doc);
  if (cur == NULL) {
    fprintf(stderr," The document is empty \n");
    xmlFreeDoc(doc);
    return;
  }

  //Reset counter
  icoeff = 0; iterm = 0; idisp = 0; istrain = 0;
  cur = cur ->xmlChildrenNode;
  while (cur != NULL) {
    if (!xmlStrcmp(cur->name, (const  xmlChar *) "Heff_definition") ||
        !xmlStrcmp(cur->name, (const  xmlChar *) "Terms_definition")){
      cur = cur->xmlChildrenNode;
    }
    if (!xmlStrcmp(cur->name, (const  xmlChar *) "coefficient")){
      //Get the name of the coefficient, need to be debug..
      uri = xmlGetProp(cur, (const  xmlChar *) "text");
      //*name_coeff = &uri;
      //Get the value of the coefficient
      uri = xmlGetProp(cur, (const  xmlChar *) "value");
      if(uri != NULL) {
        coefficient[icoeff] = strtod(uri,NULL);
      }else{
        coefficient[icoeff] = 0.0;
      }
      xmlFree(uri);
      //Get the children of coeff node
      cur2 = cur->xmlChildrenNode;
      iterm=0;
      while (cur2 != NULL){
        if (!xmlStrcmp(cur2->name, (const  xmlChar *) "term")){
          //Get the weght of the term
          uri2 = xmlGetProp(cur2, (const  xmlChar *) "weight");
          weight[iterm][icoeff] = strtod(uri2,NULL);
          xmlFree(uri2);            
          //Get the children of the term
          cur3 = cur2->xmlChildrenNode;
          idisp = 0;
          istrain = 0;
          while (cur3 != NULL){
            if (!xmlStrcmp(cur3->name, (const  xmlChar *) "displacement_diff")){
              // Get the index of the atom a
              uri2 = xmlGetProp(cur3, (const  xmlChar *) "atom_a");
              atindx[idisp][0][iterm][icoeff] = strtod(uri2,NULL);
              xmlFree(uri2);
              // Get the index of the atom b
              uri2 = xmlGetProp(cur3, (const  xmlChar *) "atom_b");
              atindx[idisp][1][iterm][icoeff] = strtod(uri2,NULL);
              xmlFree(uri2);
              
              //Get the direction
              uri2 = xmlGetProp(cur3, (const  xmlChar *) "direction");
              if(strcmp(uri2,"x") == 0){direction[idisp][iterm][icoeff] = 1;}
              if(strcmp(uri2,"y") == 0){direction[idisp][iterm][icoeff] = 2;}
              if(strcmp(uri2,"z") == 0){direction[idisp][iterm][icoeff] = 3;}
              xmlFree(uri2);
              
              //Get the power
              uri2 = xmlGetProp(cur3, (const  xmlChar *) "power");
              power_disp[idisp][iterm][icoeff] = strtod(uri2,NULL);
              
              //Get the children of the displacement
              cur4 = cur3->xmlChildrenNode;
              while (cur4 != NULL){
                if (!xmlStrcmp(cur4->name, (const  xmlChar *) "cell_a")){
                  key = xmlNodeListGetString(doc, cur4->xmlChildrenNode, 1);
                  pch = strtok(key,"\t \n");  
                  for(i=0;i<3;i++){
                    if (pch != NULL){
                      cell[idisp][0][i][iterm][icoeff]=strtod(pch,NULL);
                      pch = strtok(NULL,"\t \n");
                    }
                    }
                }
                if (!xmlStrcmp(cur4->name, (const  xmlChar *) "cell_b")){
                  key = xmlNodeListGetString(doc, cur4->xmlChildrenNode, 1);
                  pch = strtok(key,"\t \n");  
                  for(i=0;i<3;i++){
                    if (pch != NULL){
                      cell[idisp][1][i][iterm][icoeff]=strtod(pch,NULL);
                      pch = strtok(NULL,"\t \n");                      
                    }
                  }
                  }
                cur4 = cur4->next;
              }
              idisp++;
            }
            if (!xmlStrcmp(cur3->name, (const  xmlChar *) "strain")){
              uri2 = xmlGetProp(cur3, (const  xmlChar *) "power");
              power_strain[istrain][iterm][icoeff] = strtod(uri2,NULL);
              xmlFree(uri2); 
              uri2 = xmlGetProp(cur3, (const  xmlChar *) "voigt");
              strain[istrain][iterm][icoeff] = strtod(uri2,NULL); 
              xmlFree(uri2); 
              istrain++;
            }
            cur3 = cur3->next;
          }
          iterm ++;
        }
        cur2 = cur2->next;
      }
      icoeff ++;
    }    
    cur = cur->next;
  }
  xmlFreeDoc(doc);
}
  
void effpot_xml_getDimCoeff(char *filename,int *ncoeff,int *nterm_max,int *ndisp_max){
  int icoeff,idisp,iterm;
  int count1,count2;
  xmlDocPtr doc;
  xmlNodePtr cur,cur2,cur3;

  icoeff = 0;
  iterm  = 0;
  idisp  = 0;
  
  doc = xmlParseFile(filename);
  if (doc == NULL) printf(" error: could not parse file file.xml\n");

  cur = xmlDocGetRootElement(doc);
  if (cur == NULL) {
    fprintf(stderr," The document is empty \n");
    xmlFreeDoc(doc);
    return;
  }
  cur = cur ->xmlChildrenNode;
  while (cur != NULL) {
    if (!xmlStrcmp(cur->name, (const  xmlChar *) "Heff_definition") ||
        !xmlStrcmp(cur->name, (const  xmlChar *) "Terms_definition")){
      cur = cur->xmlChildrenNode;
    }
    if (!xmlStrcmp(cur->name, (const  xmlChar *) "coefficient")){
      icoeff ++;
      cur2 = cur->xmlChildrenNode;
      count1 = 0;
          while (cur2 != NULL){
            if (!xmlStrcmp(cur2->name, (const  xmlChar *) "term")){
              count1 ++;
              count2 = 0;
              cur3 = cur2->xmlChildrenNode;             
              while (cur3 != NULL){
                if (!xmlStrcmp(cur3->name, (const  xmlChar *) "displacement_diff")) {count2 ++;}
                if (!xmlStrcmp(cur3->name, (const  xmlChar *) "strain")) {count2 ++;}
                cur3 = cur3->next;
              }
              if(count2 > idisp){idisp = count2;}
            }
            cur2 = cur2->next;
          }
          if(count1 > iterm){iterm = count1;}
    }
    cur = cur->next;
  }
  xmlFreeDoc(doc);
  *ncoeff = icoeff;
  *nterm_max = iterm;
  *ndisp_max = idisp;
}


void effpot_xml_getNumberKey(char *filename,char*name_key,int*result){
  int i;
  xmlDocPtr doc;
  xmlNodePtr cur;
  
  i = 0;
  doc = xmlParseFile(filename);
  if (doc == NULL) printf(" error: could not parse file file.xml\n");

  cur = xmlDocGetRootElement(doc);
  if (cur == NULL) {
    fprintf(stderr," The document is empty \n");
    xmlFreeDoc(doc);
    return;
  }
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const  xmlChar *) name_key))) {i++;}
    cur = cur->next;
  }
  xmlFreeDoc(doc);
  *result = i;
}

void effpot_xml_getValue(char *filename,char*name_key,char*result){
  xmlDocPtr doc;
  xmlNodePtr cur;
  xmlChar *key;
  
  doc = xmlParseFile(filename);
  if (doc == NULL) printf(" error: could not parse file file.xml\n");

  cur = xmlDocGetRootElement(doc);
  if (cur == NULL) {
    fprintf(stderr," The document is empty \n");
    xmlFreeDoc(doc);
    return;
  }
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const  xmlChar *) name_key))) {
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      printf(" keyword: %s\n", key);
      //      *result = key;
      xmlFree(key);
    }
    cur = cur->next;
  }
  xmlFreeDoc(doc);
}

void effpot_xml_getAttribute(char *filename,char*name_key,char*name_attributes,char*result){  
  xmlDocPtr doc;
  xmlNodePtr cur;
  xmlChar *uri;

  doc = xmlParseFile(filename);
  if (doc == NULL) printf(" error: could not parse file file.xml\n");

  cur = xmlDocGetRootElement(doc);
  if (cur == NULL) {
    fprintf(stderr," The document is empty \n");
    xmlFreeDoc(doc);
    return;
  }
  cur = cur->xmlChildrenNode;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const  xmlChar *) name_key))) {
      uri = xmlGetProp(cur, name_attributes);
      printf("uri: %s\n", uri);
      xmlFree(uri);
    }
    cur = cur->next;
  }
  xmlFreeDoc(doc);
}

/****************************************************/
/*    Below are functions to read xml for spin part */
/****************************************************/

int xml_read_spin_system(char *fname, double *ref_energy, double *unitcell[],
                         int *natoms, double *masses[], int *nspins,
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
  *nspins = 0;

  size_t size, i;
  if (file_exists(fname)==0){
  fprintf(stderr, "xml file %s does not exist. Exit!\n", fname);
  return 1;
  }

  xmlDocPtr doc;
  xmlNodePtr cur, cur2;
  xmlChar *key;
  doc = xmlParseFile(fname);
  if (doc == NULL) {
    fprintf(stderr, "Document %s parse failed. \n", fname);
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
      //printf("  energy: %lf\n", *ref_energy);
      xmlFree(key);
    }

    // read unit cell
    if (!xmlStrcmp(cur->name, (const xmlChar *)("unit_cell"))) {
      //printf("%s\n", cur->name);
      key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
      printf("unit_cell: %s\n", key);
      string2Array((char *)key, unitcell, &size);
      xmlFree(key);
      key = xmlGetProp(cur, BAD_CAST "units");
      xmlFree(key);
    }

    // read atoms
    if (!xmlStrcmp(cur->name, BAD_CAST "atom")) {
      (*natoms)++;
      //printf("%s\n", cur->name);
      int ind_spin = -1;
      // mass
      key = xmlGetProp(cur, BAD_CAST "mass");

      //printf("mass: %s\n", key);
      insertArray(&mass_array, strtod((const char *)key, NULL));
      xmlFree(key);


      // index_spin (optional if it's not a spin xml file. )

      key = xmlGetProp(cur, BAD_CAST "index_spin");
      if (key != NULL) {
        ind_spin = (int)strtol((const char *)key, NULL, 10);
      } else {
        ind_spin = -1;
      }

      //printf("index_spin: %d\n", ind_spin);

      if (ind_spin>0){
      // gyroratio (optional)
      key = xmlGetProp(cur, BAD_CAST "gyroratio");
      //printf("gyroratio: %s\n", key);
      if (key != NULL) {
        insertArray(&gyroratio_array, strtod((const char *)key, NULL));
      }
      else if(ind_spin>0) {
        insertArray(&gyroratio_array, 0.0);
      }
      xmlFree(key);

      // damping_factors(optional)
      key = xmlGetProp(cur, BAD_CAST "damping_factor");
      //printf("damping_factor: %s\n", key);
      if (key != NULL) {
        insertArray(&damping_factor_array, strtod((const char *)key, NULL));
      }
      else if(ind_spin>0) {
        insertArray(&damping_factor_array, 1.0);
      }
      xmlFree(key);

        (*nspins)++;
      }
      insertIntArray(&index_spin_array, ind_spin);

      cur2 = cur->xmlChildrenNode;

      while (cur2 != NULL) {
        // positions
        if (!xmlStrcmp(cur2->name, BAD_CAST "position")) {
          key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
          double *pos;
          string2Array((char *)key, &pos, &size);
          for (i = 0; i < size; ++i) {
            insertArray(&position_array, pos[i]);
          }
          xmlFree(key);
        }

        // spinat, which is optional.
        if (!xmlStrcmp(cur2->name, BAD_CAST "spinat")) {
          key = xmlNodeListGetString(doc, cur2->xmlChildrenNode, 1);
          //printf("spinat: %s\n", key);
          double *spinat_tmp;
          string2Array((char *)key, &spinat_tmp, &size);
	  //for (int i=0;i<size;++i){printf("%lf ", spinat_tmp[i]);}
          if (size != 3) {
            fprintf(stderr,
                    "Error reading xml file, spinat should be a 3-vector, size is %zu", size);
          }
          for (i = 0; i < size; ++i) {
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

  //if ((int)size / 3 != *nspins) {
  //  fprintf(stderr, "Number of spinat not equal to number of magnetic atoms.\n");
  //}
  copyArraytoCArray(&gyroratio_array, gyroratios, &size);
  //if ((int)size != *nspins) {
  //  fprintf(stderr,
  //          "Number of gyroratios not equal to number of magnetic atoms");
  //}

  copyArraytoCArray(&damping_factor_array, damping_factors, &size);
  //if ((int)size != *nspins) {
  //  fprintf(stderr,
  //          "Number of damping_factors not equal to number of magnetic atoms");
  //}
  //
  //
  freeIntArray(&index_spin_array);
  freeArray(&mass_array);
  freeArray(&gyroratio_array);
  freeArray(&damping_factor_array);
  freeArray(&position_array);
  freeArray(&spinat_array);
  return 0;
}

int xml_read_spin_exchange( char * fname, int *exc_nnz, int *exc_ilist[],
                            int *exc_jlist[], int *exc_Rlist[],
                             double *exc_vallist[]){
  *exc_nnz=0;
  size_t i;
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
              for(i=0; i< size; i++)
                {
                  (*exc_vallist)[counter*3+i]=dtmp[i]*eV;
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
  size_t i;

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
              for(i=0; i< size; i++)
                {
                  (*dmi_vallist)[counter*3+i]=dtmp[i]*eV;
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
  size_t i;

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
              for(i=0; i< size; i++)
                {
                  (*uni_amplitude_list)[i]=dtmp[i]*eV;
                }
            }
            if (!xmlStrcmp(cur3->name, BAD_CAST"direction")) {
              key = xmlNodeListGetString(doc, cur3->xmlChildrenNode, 1);
              double *dtmp;
              size_t size;
              string2Array((char *)key, &dtmp, &size);
              xmlFree(key);
              for(i=0; i< size; i++)
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
  size_t i;

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
              for(i=0; i< size; i++)
                {
                  (*bi_vallist)[counter*9+i]=dtmp[i]*eV;
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
                   int *natoms, double *masses[], int *nspins,
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
  xml_read_spin_system(fname, ref_energy, unitcell, natoms, masses, nspins, index_spin, gyroratios, damping_factors, positions, spinat);

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

void xml_free_spin(char *fname, double *ref_energy, double *unitcell[9],
                   int *natoms, double *masses[], int *nspins,
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
                   double *bi_vallist[])
{
  free(*unitcell);
  free(*masses);
  free(*index_spin);
  free(*gyroratios);
  free(*damping_factors);
  free(*positions);
  free(*spinat);

  unitcell=NULL;
  masses=NULL;
  index_spin=NULL;
  gyroratios=NULL;
  damping_factors=NULL;
  positions=NULL;
  spinat=NULL;

  if (*exc_nnz !=0){
    free(*exc_ilist);
    free(*exc_jlist);
    free(*exc_Rlist);
    free(*exc_vallist);
    exc_ilist=NULL;
    exc_jlist=NULL;
    exc_Rlist=NULL;
    exc_vallist=NULL;
  }

  if (*dmi_nnz !=0){
    free(*dmi_ilist);
    free(*dmi_jlist);
    free(*dmi_Rlist);
    free(*dmi_vallist);
    dmi_ilist=NULL;
    dmi_jlist=NULL;
    dmi_Rlist=NULL;
    dmi_vallist=NULL;

  }

  if (*uni_nnz!=0){
    free(*uni_ilist);
    free(*uni_amplitude_list);
    free(*uni_direction_list);
    uni_ilist=NULL;
    uni_amplitude_list=NULL;
    uni_direction_list=NULL;
  }

  if (*bi_nnz !=0){
    free(*bi_ilist);
    free(*bi_jlist);
    free(*bi_Rlist);
    free(*bi_vallist);
    bi_ilist=NULL;
    bi_jlist=NULL;
    bi_Rlist=NULL;
    bi_vallist=NULL;
  }
}

// This function is for testing.
// TODO hexu: to be removed
int test_read_xml() {

  char * fname="test_f.xml";
  double ref_energy;
  double *unitcell, *masses, *gyroratios, *damping_factors, *positions, *spinat;
  int natoms, nspins, *index_spin;
  int i, j;
  xml_read_spin_system("test_f.xml", &ref_energy, &unitcell, &natoms, &masses,
   &nspins, &index_spin, &gyroratios, &damping_factors, &positions, &spinat);

  // exchange
  printf("======Exchange Terms========\n");
  int exc_nnz, *exc_ilist, *exc_jlist, *exc_Rlist;
  double *exc_vallist;
  xml_read_spin_exchange(fname, &exc_nnz, &exc_ilist, &exc_jlist, &exc_Rlist, &exc_vallist);
  for(i=0;i<exc_nnz;i++){
    printf("%d\t%d\t%d\t%d\t%d\t : %E\t%E\t%E\n", exc_ilist[i], exc_jlist[i], exc_Rlist[3*i], exc_Rlist[3*i+1], exc_Rlist[3*i+2], exc_vallist[3*i], exc_vallist[3*i+1], exc_vallist[3*i+2]);
  }

  //dmi
  printf("======DMI Terms========\n");
  int dmi_nnz, *dmi_ilist, *dmi_jlist, *dmi_Rlist;
  double *dmi_vallist;
  xml_read_spin_dmi(fname, &dmi_nnz, &dmi_ilist, &dmi_jlist, &dmi_Rlist, &dmi_vallist);
  for(i=0;i<dmi_nnz;i++){
    printf("%d\t%d\t%d\t%d\t%d\t :", dmi_ilist[i], dmi_jlist[i], dmi_Rlist[3*i], dmi_Rlist[3*i+1], dmi_Rlist[3*i+2]);
    for (j=0; j < 3; ++j) {
      printf("%E\t", dmi_vallist[i*3+j]);
    }
    printf("\n");
  }

  //uniaxial SIA
  printf("======uniaxial SIA Terms========\n");
  int uni_nnz, *uni_ilist;
  double *uni_amplitude_list, *uni_direction_list;
  xml_read_spin_uni(fname, &uni_nnz, &uni_ilist, &uni_amplitude_list, &uni_direction_list);
  for(i=0;i<uni_nnz;i++){
    printf("%d\t :", uni_ilist[i]);
    printf("%E\t:", uni_amplitude_list[i]);
    for (j=0; j < 3; ++j) {
      printf("%E\t", uni_direction_list[i*3+j]);
    }
    printf("\n");
  }


  //bilinear
  printf("======Bilinear Terms========\n");
  int bi_nnz, *bi_ilist, *bi_jlist, *bi_Rlist;
  double *bi_vallist;
  xml_read_spin_bilinear("test_f.xml", &bi_nnz, &bi_ilist, &bi_jlist, &bi_Rlist, &bi_vallist);
  for(i=0;i<bi_nnz;i++){
    printf("%d\t%d\t%d\t%d\t%d\t :", bi_ilist[i], bi_jlist[i], bi_Rlist[3*i], bi_Rlist[3*i+1], bi_Rlist[3*i+2]);
    for (j=0; j < 9; ++j) {
      printf("%E\t", bi_vallist[i*9+j]);
    }
    printf("\n");
  }


  return 0;
}

#else
int xml_read_spin_system(char *fname, double *ref_energy, double *unitcell[],
                         int *natoms, double *masses[], int *nspins,
                         int *index_spin[], double *gyroratios[],
                         double *damping_factors[],
                         double *positions[], double *spinat[]) 
{
	fprintf(stderr, "Cannot read xml file. Please install abinit with libxml enabled.\n");
	exit(1);
	return 1;
}



void xml_read_spin(char *fname, double *ref_energy, double *unitcell[9],
                   int *natoms, double *masses[], int *nspins,
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
                   double *bi_vallist[])
{
	fprintf(stderr, "Cannot read xml file. Please install abinit with libxml enabled.\n");
	exit(1);
}


void xml_free_spin(char *fname, double *ref_energy, double *unitcell[9],
                   int *natoms, double *masses[], int *nspins,
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
                   double *bi_vallist[])
{
	fprintf(stderr, "Cannot read xml file. Please install abinit with libxml enabled.\n");
	exit(1);
}


#endif
