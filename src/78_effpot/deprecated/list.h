// TODO this file can be removed. It is used in reading xml effective potential.
// TODO hexu: should we move this to somewhere and make it more general & useful

#include <stdio.h>
#include <stdlib.h>

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

void copyArraytoCArray(Array *l, double **a, size_t *size){
  *size=l->used;
  *a=(double *) malloc(sizeof(double)*l->used);
  for(size_t i=0;i<l->used;i++){
    (*a)[i]=l->array[i];
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
  *size=l->used;
  *a=(int *) malloc(sizeof(int)*l->used);
  for(size_t i=0;i<l->used;i++){
    (*a)[i]=l->array[i];
  }
}
