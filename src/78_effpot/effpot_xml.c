/*
 * Copyright (C) 2015-2017 ABINIT group (AM)
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

#if defined HAVE_LIBXML

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>

//define type for dynamical array 
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


void effpot_xml_checkXML(char *filename,char *name_xml){  
  xmlDocPtr doc;
  xmlNodePtr cur;
  xmlNodeSetPtr nodeset;

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

  present  = 0;
  iatom    = 0;
  irpt1     = 0;
  irpt2     = 0;
  iqpt    = 0;
  itypat   = 0;
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
 
  freeArray(&typat);

  *natom  = iatom;
  *nqpt   = iqpt;
  *loc_nrpt   = irpt1;
  *tot_nrpt   = irpt2;
  *ntypat = itypat;
}

void effpot_xml_readSystem(char *filename,int *natom,int *ntypat,int *nrpt,int *nqpt,
                           double amu[*ntypat],double atmfrc[*nrpt][*natom][3][*natom][3][2],
                           int cell[*nrpt][3],double dynmat[*nqpt][*natom][3][*natom][3][2],
                           double elastic_constants[6][6],
                           double *energy,double epsilon_inf[3][3],
                           double ewald_atmfrc[*nrpt][*natom][3][*natom][3][2],
                           double phfrq[*nqpt][3* *natom],
                           double rprimd[3][3],double qph1l[*nqpt][3],
                           double short_atmfrc[*nrpt][*natom][3][*natom][3][2],
                           int typat[*natom],double xcart[*natom][3],double zeff[*natom][3][3]){
  xmlDocPtr doc;
  char *pch;
  double total_atmfrc[*nrpt][*natom][3][*natom][3][2];
  double local_atmfrc[*nrpt][*natom][3][*natom][3][2];
  int cell_local[*nrpt][3];
  int cell_total[*nrpt][3];
  int iatom,iamu,irpt1,irpt2,irpt3,iqpt,present;
  int ia,ib,mu,nu,voigt;
  int i,j;
  xmlNodePtr cur,cur2,cur3;
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
      for(i=0;i<=*ntypat;i++){
        if(amu[i]==strtod(uri,NULL)){
          present = 1;
          break;
        }
      }
      if(present==0){
        amu[iamu]=strtod(uri,NULL);
        iamu++;
      }
      // fill the typat table
      for(i=0;i<=*ntypat;i++){
        if(amu[i]==strtod(uri,NULL)){
          typat[iatom]=i+1;
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
                      local_atmfrc[irpt1][ib][nu][ia][mu][0]=strtod(pch,NULL);
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
                      total_atmfrc[irpt2][ib][nu][ia][mu][0]=strtod(pch,NULL);
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

  //Reorder the ATMFRC
  //Case 1: only local in the xml
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
              atmfrc[i][ib][nu][ia][mu][0]=local_atmfrc[i][ib][nu][ia][mu][0];
              short_atmfrc[i][ib][nu][ia][mu][0]=local_atmfrc[i][ib][nu][ia][mu][0];
              ewald_atmfrc[i][ib][nu][ia][mu][0]=0.0;
              //Set imaginary part to 0
              short_atmfrc[i][ib][nu][ia][mu][1]= 0.0;
              atmfrc[i][ib][nu][ia][mu][1] = 0.0;      
              ewald_atmfrc[i][ib][nu][ia][mu][1]= 0.0;
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
              atmfrc[i][ib][nu][ia][mu][0]=total_atmfrc[i][ib][nu][ia][mu][0];
              short_atmfrc[i][ib][nu][ia][mu][0]=0.0;
              ewald_atmfrc[i][ib][nu][ia][mu][0]=total_atmfrc[i][ib][nu][ia][mu][0];
              //Set imaginary part to 0
              short_atmfrc[i][ib][nu][ia][mu][1]= 0.0;
              atmfrc[i][ib][nu][ia][mu][1] = 0.0;      
              ewald_atmfrc[i][ib][nu][ia][mu][1]= 0.0;

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
                atmfrc[i][ib][nu][ia][mu][0] = total_atmfrc[i][ib][nu][ia][mu][0];
                ewald_atmfrc[i][ib][nu][ia][mu][0]= atmfrc[i][ib][nu][ia][mu][0]-
                                                    short_atmfrc[i][ib][nu][ia][mu][0];
                //Set imaginary part to 0
                atmfrc[i][ib][nu][ia][mu][1] = 0.0; 
                ewald_atmfrc[i][ib][nu][ia][mu][1]= 0.0;
                for(j=0;j<irpt1;j++){
                  if(cell_local[j][0] == cell[i][0] && 
                     cell_local[j][1] == cell[i][1] &&
                     cell_local[j][2] == cell[i][2] ){
                    if(ia==0 && ib==0 && mu==0 && nu==0){irpt3++;}
                    short_atmfrc[i][ib][nu][ia][mu][0]= local_atmfrc[j][ib][nu][ia][mu][0];
                    short_atmfrc[i][ib][nu][ia][mu][1]= 0.0;
                  }
                }
              }    
            }
          }
        }
      }
      if(irpt3 != irpt1){
        fprintf(stdout,"\n WARNING: The number of local and total rpt are not equivalent\n");
        fprintf(stdout,"          in the XML file :%d %d\n",irpt1,irpt3);
        fprintf(stdout,"          The missing local IFC will be set to zero\n");        
      }
    }
    else{
      fprintf(stderr,"error: Local rpt is superior to total rpt in the XML file:%d %d\n",\
              irpt1,irpt2);
      exit(0);
    }
  }else{
    fprintf(stderr,"error: Number of local and total rpt doesn't match with the XML file:%d %d\n",\
            irpt1,irpt2);
    exit(0);
  }
}


void effpot_xml_getDimStrainCoupling(char *filename, int *nrpt,int *voigt){
  xmlDocPtr doc;
  int i,irpt;
  xmlNodePtr cur,cur2,cur3;
  xmlChar *key,*uri;

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
}

void effpot_xml_readStrainCoupling(char *filename,int *natom,int *nrpt,int *voigt,
                                   double elastic3rd[6][6], double elastic_displacement[*natom][3][6],
                                   double internal_strain[*natom][3],
                                   double phonon_strain_atmfrc[*nrpt][*natom][3][*natom][3],
                                   int phonon_strain_cell[*nrpt][3]){
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
                      phonon_strain_cell[irpt][i]=atoi(pch);
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
                          int power[*ndisp][*nterm][*ncoeff],double weight[*nterm][*ncoeff]){
  int i,idisp,j,iterm,jterm,icoeff;
  xmlDocPtr doc;
  char * pch;
  xmlNodePtr cur,cur2,cur3,cur4,cur5;
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
        power[idisp][iterm][icoeff] = 0;
        for (i=0;i<2;i++){
          atindx[idisp][i][iterm][icoeff] = 0;
          for (j=0;j<3;j++){
            cell[idisp][i][j][iterm][icoeff] = 0;
          }
        }      
      }
    }
  }

  jterm = 0;
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
  icoeff = 0; iterm = 0; idisp = 0;
  cur = cur ->xmlChildrenNode;
  while (cur != NULL) {
    if (!xmlStrcmp(cur->name, (const  xmlChar *) "Heff_definition") ||
        !xmlStrcmp(cur->name, (const  xmlChar *) "Terms_definition")) {
      cur2 = cur->xmlChildrenNode;
      while (cur2 != NULL){
        if (!xmlStrcmp(cur2->name, (const  xmlChar *) "coefficient")) {
          //Get the name of the coefficient, need to be debug..
          uri = xmlGetProp(cur2, (const  xmlChar *) "text");
          //*name_coeff = &uri;
          //Get the value of the coefficient
          uri = xmlGetProp(cur2, (const  xmlChar *) "value");
          if(uri != NULL) {
            coefficient[icoeff] = strtod(uri,NULL);
          }else{
            coefficient[icoeff] = 0.0;
          }        
          xmlFree(uri);
          //Get the children of coeff node
          cur3 = cur2->xmlChildrenNode;
          iterm=0;
          while (cur3 != NULL){
            if (!xmlStrcmp(cur3->name, (const  xmlChar *) "term")){
             //Get the weght of the term
              uri2 = xmlGetProp(cur3, (const  xmlChar *) "weight");
              weight[iterm][icoeff] = strtod(uri2,NULL);
              xmlFree(uri2);            
             //Get the children of the term
              cur4 = cur3->xmlChildrenNode;
              idisp = 0;
              while (cur4 != NULL){
                if (!xmlStrcmp(cur4->name, (const  xmlChar *) "displacement_diff")){
                // Get the index of the atom a
                  uri2 = xmlGetProp(cur4, (const  xmlChar *) "atom_a");
                  atindx[idisp][0][iterm][icoeff] = strtod(uri2,NULL);
                  xmlFree(uri2);
                  // Get the index of the atom b
                  uri2 = xmlGetProp(cur4, (const  xmlChar *) "atom_b");
                  atindx[idisp][1][iterm][icoeff] = strtod(uri2,NULL);
                  xmlFree(uri2);

                  //Get the direction
                  uri2 = xmlGetProp(cur4, (const  xmlChar *) "direction");
                  if(strcmp(uri2,"x") == 0){direction[idisp][iterm][icoeff] = 1;}
                  if(strcmp(uri2,"y") == 0){direction[idisp][iterm][icoeff] = 2;}
                  if(strcmp(uri2,"z") == 0){direction[idisp][iterm][icoeff] = 3;}
                  xmlFree(uri2);

                  //Get the power
                  uri2 = xmlGetProp(cur4, (const  xmlChar *) "power");
                  power[idisp][iterm][icoeff] = strtod(uri2,NULL);
                  
                  //Get the children of the displacement
                  cur5 = cur4->xmlChildrenNode;
                  while (cur5 != NULL){
                    if (!xmlStrcmp(cur5->name, (const  xmlChar *) "cell_a")){
                      key = xmlNodeListGetString(doc, cur5->xmlChildrenNode, 1);
                      pch = strtok(key,"\t \n");  
                      for(i=0;i<3;i++){
                        if (pch != NULL){
                          cell[idisp][0][i][iterm][icoeff]=strtod(pch,NULL);
                          pch = strtok(NULL,"\t \n");
                        }
                      }
                    }
                    if (!xmlStrcmp(cur5->name, (const  xmlChar *) "cell_b")){
                      key = xmlNodeListGetString(doc, cur5->xmlChildrenNode, 1);
                      pch = strtok(key,"\t \n");  
                      for(i=0;i<3;i++){
                        if (pch != NULL){
                          cell[idisp][1][i][iterm][icoeff]=strtod(pch,NULL);
                          pch = strtok(NULL,"\t \n");                      
                        }
                      }
                    }
                    cur5 = cur5->next;
                  }
                  idisp++;
                }
                if (!xmlStrcmp(cur4->name, (const  xmlChar *) "strain")){
                  uri2 = xmlGetProp(cur4, (const  xmlChar *) "power");
                  power[idisp][iterm][icoeff] = strtod(uri2,NULL);
                  xmlFree(uri2); 
                  uri2 = xmlGetProp(cur4, (const  xmlChar *) "voigt");
                  direction[idisp][iterm][icoeff] = -1 *  strtod(uri2,NULL); 
                  xmlFree(uri2); 
                 //Set to -1 the useless quantitiers for strain                       
                  for(i=0;i<2;i++){
                    atindx[idisp][i][iterm][icoeff]  = -1 ;
                    for(j=0;j<3;j++){
                      cell[idisp][i][j][iterm][icoeff]= -1;
                    }
                  }
                  idisp++;
                }
                cur4 = cur4->next;
              }
              iterm ++;
            }
            cur3 = cur3->next;
          }
          icoeff ++;
        }
        cur2 = cur2->next;
      }
    }
    cur = cur ->next;
  }
  xmlFreeDoc(doc);
}
  
void effpot_xml_getDimCoeff(char *filename,int*ncoeff,char *nterm_max,int*ndisp_max){
  int icoeff,idisp,iterm;
  int count1,count2;
  xmlDocPtr doc;
  char * pch;
  xmlNodePtr cur,cur2,cur3,cur4;

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
        !xmlStrcmp(cur->name, (const  xmlChar *) "Terms_definition")) {
      cur2 = cur->xmlChildrenNode;
      while (cur2 != NULL){
        if (!xmlStrcmp(cur2->name, (const  xmlChar *) "coefficient")){
          icoeff ++;
          cur3 = cur2->xmlChildrenNode;
          count1 = 0;
          while (cur3 != NULL){
            if (!xmlStrcmp(cur3->name, (const  xmlChar *) "term")){
              count1 ++;
              count2 = 0;
              cur4 = cur3->xmlChildrenNode;             
              while (cur4 != NULL){
                if (!xmlStrcmp(cur4->name, (const  xmlChar *) "displacement_diff")) {count2 ++;}
                if (!xmlStrcmp(cur4->name, (const  xmlChar *) "strain")) {count2 ++;}
                cur4 = cur4->next;
              }
              if(count2 > idisp){idisp = count2;}
            }
            cur3 = cur3->next;
          }
          if(count1 > iterm){iterm = count1;}
        }
        cur2 = cur2->next;        
      }
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
  xmlChar *key, *uri;

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

#endif
