/*
 * Author: Theo Cavignac
 * Implementation of a simple pair list structure.
 * Unless you really know what you are doing you probably should
 * look at the Fortran interface instead: m_pair_list.F90
 *
 * Possible improvement:
 * - Add new type to be stored
 * - Optimise memory allocation
 */
#include <stdlib.h>

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#define FALSE 0
#define TRUE 1

/* type codes */
#define TC_EMPTY -2
#define TC_NOTFOUND -1
#define TC_INT 0
#define TC_REAL 1
#define TC_STRING 2

typedef uint8_t bool;

typedef union {
  int i;
  double r;
  char* s;
/*
  int* iarr;
  double* rarr;
*/
} value;

typedef struct Pair {
  int8_t type_code;
  char* key;
  struct Pair* next;
  value val;
} pair_t;

typedef struct {
  pair_t* first;
  pair_t* cursor;
  int length;
} pair_list;

/* Lazy string comparision
 */
static bool str_eq(char* s1, char* s2){
  while(*s1 && *s2){
    if(*s1 != *s2){
      return FALSE;
    }
    s1++; s2++;
  }
  return (*s1 | *s2) == 0;
}

/* Fortran style to C style string
 * return a new allocated pointer to the C string
 */
static char* ftoc_str(char* fstr, int length){
  int i;
  char* cstr;
  cstr = malloc((length+1)*sizeof(char));
  for(i = 0; i < length; i++){
    cstr[i] = fstr[i];
  }
  cstr[length] = '\x00';
  return cstr;
}

/* C style to Fortran style string
 * fstr must be allocated and of size length
 */
static void ctof_str(char* fstr, char* cstr, int length){
  int c = 0;
  while(cstr[c] != '\x00' && c < length){
    fstr[c] = cstr[c];
    c++;
  }
  while(c < length){
    fstr[c++] = ' ';
  }
}

/* Look for a given key, if found return FALSE and have
 * selected point to the pair, else return TRUE, allocate a new pair
 * and have selected point to it
 */
static bool get_or_create(pair_list* pl, char* ckey, pair_t** selected){
  if(pl->first){
    pair_t* prev = NULL;
    pair_t* pair = pl->first;
    pair_t* new_pair = malloc(sizeof(pair_t));
    while(pair){
      if(str_eq(ckey, pair->key)){
        *selected = pair;
        return FALSE;
      } else {
        prev = pair;
        pair = pair->next;
      }
    }
    new_pair->key = ckey;
    new_pair->next = NULL;
    prev->next = new_pair;
    *selected = new_pair;
    return TRUE;
  } else {  /* first element of the list */
    pair_t* new_pair = malloc(sizeof(pair_t));
    new_pair->key = ckey;
    new_pair->next = NULL;
    pl->first = new_pair;
    pl->cursor = new_pair;
    *selected = new_pair;
    return TRUE;
  }
}

/* free a pair after freeing the next one
 */
static void pair_free(pair_t* p){
  if(p){
    pair_free(p->next);
    free(p->key);
    if(p->type_code == TC_STRING){
      free(p->val.s);
    }
    free(p);
  }
}


/* Visible from fortran */

/* set an integer value
 */
void pair_list_seti(pair_list* l, char* fkey, int* i, int* len){
  pair_t* pair = NULL;
  char* ckey = ftoc_str(fkey, *len);
  bool is_new = get_or_create(l, ckey, &pair);
  if(!is_new){
    free(ckey);
  } else {
    if(pair->type_code == TC_STRING){
      free(pair->val.s);
    }
  }
  l->length += is_new;
  pair->type_code = TC_INT;
  pair->val.i = *i;
}

/* set a real (double) value
 */
void pair_list_setr(pair_list* l, char* fkey, double* r, int* len){
  pair_t* pair = NULL;
  char* ckey = ftoc_str(fkey, *len);
  bool is_new = get_or_create(l, ckey, &pair);
  if(!is_new){
    free(ckey);
  } else {
    if(pair->type_code == TC_STRING){
      free(pair->val.s);
    }
  }
  l->length += is_new;
  pair->type_code = TC_REAL;
  pair->val.r = *r;
}

/* set a string value
 */
void pair_list_sets(pair_list* l, char* fkey, char* s, int* len, int* len_s){
  pair_t* pair = NULL;
  char* ckey = ftoc_str(fkey, *len);
  bool is_new = get_or_create(l, ckey, &pair);
  if(!is_new){
    free(ckey);
  } else {
    if(pair->type_code == TC_STRING){
      free(pair->val.s);
    }
  }
  l->length += is_new;
  pair->type_code = TC_STRING;
  pair->val.s = ftoc_str(s, *len_s);
}

/* get a value from a key
 */
void pair_list_get_(pair_list* l, char* fkey, int* type_code, int*i, double* r, char* s, int* len, int* len_s){
  if(!l->first){
    /* list is empty */
    *type_code = TC_EMPTY;
    return;
  } else {
    char* ckey = ftoc_str(fkey, *len);
    pair_t* pair = l->first;
    while(pair){
      if(str_eq(pair->key, ckey)){
        *type_code = pair->type_code;
        switch(pair->type_code){
          case TC_REAL:
            *r = pair->val.r;
            break;
          case TC_INT:
            *i = pair->val.i;
            break;
          case TC_STRING:
            ctof_str(s, pair->val.s, *len_s);
            break;
        }
        break;
      } else {
        pair = pair->next;
      }
    }
    if(!pair){
      /* key not found */
      *type_code = TC_NOTFOUND;
    }
    free(ckey);
  }
}

/* move the cursor forward in the chain
 */
void pair_list_next(pair_list* pl){
  pl->cursor = pl->cursor->next;
}

/* free the whole chained list
 */
void pair_list_free(pair_list* pl){
  pair_free(pl->first);
  pl->first = NULL;
  pl->cursor = NULL;
  pl->length = 0;
}

/* Return the pair pointed by the cursor
 */
void pair_list_look_(pair_list* pl, char* fkey, int* type_code, int* i, double* r, char* s, int* len, int* len_s){
  pair_t* p = pl->cursor;
  if(p){
    *type_code = p->type_code;
    switch(p->type_code){
      case TC_REAL:
        *r = p->val.r;
        break;
      case TC_INT:
        *i = p->val.i;
        break;
      case TC_STRING:
        ctof_str(s, p->val.s, *len_s);
        break;
    }
    ctof_str(fkey, p->key, *len);
  } else {
    /* reached end of list */
    *type_code = TC_EMPTY;
  }
}
