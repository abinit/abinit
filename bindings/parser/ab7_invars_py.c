#include <config.h>

#include <Python.h>

#include <structmember.h>
#include <string.h>
#include <stdlib.h>

/** #define HAVE_NUMPY **/

#ifdef HAVE_NUMPY
#include <numpy/arrayobject.h>
#else
#include <numarray/arrayobject.h>
#endif

#include <pthread.h>

#include "ab7_invars.h"

#ifdef DEBUG_MODE
#define DBG_printf if (DEBUG_MODE) printf
#else
#define DBG_printf if (0) printf
#endif

#include "ab7_invars_py.h"


/* proptotypes for the replacement functions */
static char *
my_strndup (char const *s, size_t n);

static size_t
my_strnlen (const char *string, size_t maxlen);


/* The abinit module provides:
   - a new class dtsets, with an array of datasets.
*/

static PyObject* dtsets_get_ids(PyObject *self, PyObject *args)
{
  Py_INCREF(DictIds);
  return DictIds;
}

typedef struct Dtsets
{
  PyObject_HEAD
  Ab7Invars *dt;
  PyObject *filename;
} Dtsets;

static void dtsets_free(Dtsets* self)
{
  DBG_printf("Deallocate a Dtsets object %p (%p).\n", self, self->dt);
  if (self->dt)
    {
      DBG_printf(" | %d.\n", *self->dt);
      ab7_invars_free(self->dt);
    }
  Py_XDECREF(self->filename);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* dtsets_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  Dtsets *self;

  self = (Dtsets*)type->tp_alloc(type, 0);
  if (self != NULL)
    {
      DBG_printf("New dtset object %p.\n", self);
      self->dt = (Ab7Invars*)0;
      self->filename = PyString_FromString("");
      if (self->filename == NULL)
	{
	  Py_DECREF(self);
	  return NULL;
	}
    }

  return (PyObject *)self;
}

static void* _parse_ab(void *filename)
{
  DBG_printf("argument is '%s'.\n", (char*)filename);
  return ab7_invars_new_from_file((char*)filename);
}

/* These are static variables, used by the second thread to
   communicate with the main thread. */
static int _parse_error;
static char* _parse_message;

static int dtsets_init(Dtsets *self, PyObject *args, PyObject *kwds)
{
  const char* filename;
  static char *kwlist[] = {"filename", NULL};
  pthread_t parse_ab;
  Ab7Invars *dt;

  if (! PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, 
				    &filename))
    return -1;


  Py_BEGIN_ALLOW_THREADS
    pthread_create(&parse_ab, NULL, _parse_ab, (void*)filename);
  pthread_join(parse_ab, (void**)&dt);
  DBG_printf("Parsing thread return %d.\n", _parse_error);
  Py_END_ALLOW_THREADS
    
    if (dt == PTHREAD_CANCELED)
      {
	PyErr_SetString(PyExc_RuntimeError, "Abinit thread cancelled.");
	return -1;
      }

  switch (_parse_error)
    {
    case (1):
      PyErr_SetString(PyExc_SyntaxError, _parse_message);
      free(_parse_message);
      return -1;
    case (2):
      PyErr_Format(PyExc_TypeError, "File '%s' is not an ABINIT file.", filename);
      return -1;
    case (3):
      PyErr_SetString(PyExc_SystemExit, "Abinit called leave_new() routine.");
      return -1;
    default:
      break;
    }

  self->dt = dt;
  DBG_printf("Parse OK -> %d.\n", *self->dt);

  Py_XDECREF(self->filename);
  self->filename = PyString_FromString(filename);
  
  return 0;
}

static PyObject* dtsets_get_ndtset(Dtsets *self)
{
  int ndtsets;
  Ab7Error error;
  PyObject *o;

  if (self->dt == NULL)
    {
      PyErr_SetString(PyExc_AttributeError, "dt");
      return NULL;
    }

  DBG_printf("dtsets.get_ndtset called dt = %d.\n", *self->dt);

  error = ab7_invars_get_ndtset(self->dt, &ndtsets);
  if (error != AB7_NO_ERROR)
    {
      PyErr_SetString(PyExc_RuntimeError, "Can't get number of datasets.");
      return NULL;
    }
  DBG_printf(" | %d\n", ndtsets);
  
  o = PyInt_FromLong(ndtsets);
  Py_INCREF(o);
  return (PyObject*)o;
}

static PyObject* dtsets_get(Dtsets *self, PyObject *args)
{
  int idtset, id, type, rdInt, i, tmp;
  size_t n;
  double rdDbl;
  Ab7Error error;
  PyObject *tuple, *o;
  const char *strId;
  int ndims, dims[7];
#ifdef HAVE_NUMPY
  npy_intp pydims[7];
#endif

  if (self->dt == NULL)
    {
      PyErr_SetString(PyExc_AttributeError, "dt");
      return NULL;
    }

  DBG_printf("dtsets.get called dt = %d.\n", *self->dt);

  /* Read the string id and get it from the dictionnary. */
  if (!PyArg_ParseTuple(args, "si", &strId, &idtset))
    return NULL;

  tuple = PyDict_GetItemString(DictIds, strId);
  if (!tuple)
    {
      PyErr_BadArgument();
      return NULL;
    }
  
  o = PyTuple_GET_ITEM(tuple, 0);
  id = (int) ((PyIntObject*)o)->ob_ival;
  o = PyTuple_GET_ITEM(tuple, 1);
  type = (int) ((PyIntObject*)o)->ob_ival;
  DBG_printf("Get value for '%s' (%d), type %d.\n", strId, id, type);

  switch (type)
    {
    case (_INT_SCALAR):
      error = ab7_invars_get_integer(self->dt, id, idtset, &rdInt);
      if (error != AB7_NO_ERROR)
	{
	  PyErr_Format(PyExc_RuntimeError, "Can't get integer value from '%s'.", strId);
	  return NULL;
	}
      o = PyInt_FromLong(rdInt);
      Py_INCREF(o);
      return (PyObject*)o;
    case (_INT_ARRAY):
    case (_DOUBLE_ARRAY):
      /* We get the shape. */
      error = ab7_invars_get_shape(self->dt, &ndims, dims, id, idtset);
      DBG_printf("Shape is %d - %d %d %d %d %d %d %d\n",
		 ndims, dims[0], dims[1], dims[2], dims[3],
		 dims[4], dims[5], dims[6]);
      if (error != AB7_NO_ERROR)
	{
	  PyErr_SetString(PyExc_RuntimeError, "Can't get shape of array.");
	  return NULL;
	}
      n = 1;
      for (i = 0; i < ndims; i++) n *= dims[i];
      DBG_printf("Total size is %d.\n", (int)n);
      if (n > 0)
	{
#ifdef HAVE_NUMPY
	  /* Copy and reverse dims into pydims. */
	  for (i = 0; i < ndims; i++)
	    pydims[i] = (npy_intp)dims[ndims - i - 1];
#else
	  /* Reverse only. */
	  for (i = 0; i < ndims / 2; i++)
	    {tmp = dims[ndims - i - 1]; dims[ndims - i - 1] = dims[i]; dims[i] = tmp;}
#endif
	  if (type == _INT_ARRAY)
	    {
#ifdef HAVE_NUMPY
	      o = PyArray_SimpleNew(ndims, pydims, PyArray_INT);
#else
	      o = PyArray_FromDims(ndims, dims, PyArray_INT);
#endif
	      DBG_printf("PyArray allocation OK (%d).\n", (int)n);
	      error = ab7_invars_get_integer_array(self->dt,
#ifdef HAVE_NUMPY
						   (int*)(PyArray_DATA(o)),
#else
						   (int*)((PyArrayObject*)o)->data,
#endif
						   n, id, idtset);
	      DBG_printf("Read %d.\n", error);
	    }
	  else
	    {
#ifdef HAVE_NUMPY
	      o = PyArray_SimpleNew(ndims, pydims, PyArray_DOUBLE);
#else
	      o = PyArray_FromDims(ndims, dims, PyArray_DOUBLE);
#endif
	      DBG_printf("PyArray allocation OK (%d).\n", (int)n);
	      error = ab7_invars_get_real_array(self->dt,
#ifdef HAVE_NUMPY
						(double*)(PyArray_DATA(o)),
#else
						(double*)((PyArrayObject*)o)->data,
#endif
						n, id, idtset);
	      DBG_printf("Read %d.\n", error);
	    }
	  if (error != AB7_NO_ERROR)
	    {
	      Py_DECREF(o);
	      PyErr_SetString(PyExc_RuntimeError, "Can't get values of array.");
	      return NULL;
	    }      
	}
      else
	{
	  DBG_printf("Create null dimension array.\n");
	  dims[0] = 0;
	  if (type == _INT_ARRAY)
	    o = PyArray_FromDims(1, dims, PyArray_INT);
	  else
	    o = PyArray_FromDims(1, dims, PyArray_DOUBLE);
	  DBG_printf("OK.\n");
	}
      Py_INCREF(o);
      return (PyObject*)o;
    case (_DOUBLE_SCALAR):
      error = ab7_invars_get_real(self->dt, id, idtset, &rdDbl);
      if (error != AB7_NO_ERROR)
	{
	  PyErr_SetString(PyExc_NotImplementedError, "Can't get number of datasets.");
	  return NULL;
	}
      o = PyFloat_FromDouble(rdDbl);
      Py_INCREF(o);
      return (PyObject*)o;
    }
  
  Py_RETURN_NONE;
}

/* We define here all the methods of the module. */
static PyMethodDef ModuleMethods[] =
  {
    {"get_ids", (PyCFunction)dtsets_get_ids, METH_NOARGS,
     "Return a dictionary with the list of available attributes."},
    {NULL}
  };
/* We define here all the methods of the dtsets class. */
static PyMethodDef DtsetsMethods[] =
  {
    {"get_ndtset", (PyCFunction)dtsets_get_ndtset, METH_NOARGS,
     "Return the number of parsed datasets (not counting the default values)."},
    {"get", (PyCFunction)dtsets_get, METH_VARARGS,
     "Return the values of an attribute."},
    {NULL, NULL, 0, NULL}
  };
static PyMemberDef DtsetsMembers[] =
  {
/*     {"dt", T_OBJECT_EX, offsetof(Dtsets, dt), READONLY, */
/*      "Internal identifier for Fortran calls."}, */
    {"filename", T_OBJECT_EX, offsetof(Dtsets, filename), READONLY,
     "The path to the read file."},
    {NULL, 0, 0, 0, NULL}
  };
static PyTypeObject DtsetsType = {
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "ab7_invars.Dtsets",           /*tp_name*/
  sizeof(Dtsets),            /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)dtsets_free,   /*tp_dealloc*/
  0,                         /*tp_print*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  0,                         /*tp_repr*/
  0,                         /*tp_as_number*/
  0,                         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  0,                         /*tp_hash */
  0,                         /*tp_call*/
  0,                         /*tp_str*/
  0,                         /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  "Dtsets objects",          /* tp_doc */
  0,                         /* tp_traverse */
  0,                         /* tp_clear */
  0,                         /* tp_richcompare */
  0,                         /* tp_weaklistoffset */
  0,                         /* tp_iter */
  0,                         /* tp_iternext */
  DtsetsMethods,             /* tp_methods */
  DtsetsMembers,             /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)dtsets_init,     /* tp_init */
  0,                         /* tp_alloc */
  dtsets_new,                /* tp_new */
};

/* The name initab7_invars is mandatory, do not change it! */
PyMODINIT_FUNC initab7_invars(void)
{
  PyObject *module;

  if (PyType_Ready(&DtsetsType) < 0)
    return;

  ab7_mpi_init();

  module = Py_InitModule("ab7_invars", ModuleMethods);

  if (module == NULL)
    return;

  /* Use Numpy. */
  import_array()

  Py_INCREF(&DtsetsType);
  PyModule_AddObject(module, "Dtsets", (PyObject *)&DtsetsType);

  _init_dict_ids();
}

#ifndef FC_FUNC
#define FC_FUNC(A,B) A##_
#endif
void FC_FUNC(wrtout, WRTOUT)(int *unit,char message[500])
{
  char *buf, *ptError, *ptInvars0, *ptInstrng, *ptSize, *pt;

  buf = my_strndup(message, 500);
  pt = buf + 498;
/*   while (buf[0] == ' ') buf++; */
  while (pt[0] == ' ' || pt[0] == '\0') pt--;
  pt[1] = '\0';

  /* We analyse buf. If, it contains an error, we test if it is
     about natom in inarvs0. If so, the file is not a valid ABINIT
     file. On the contrary, we get the message and raise an error. */
  ptError = strstr(buf, "ERROR");
  ptInvars0 = strstr(buf, "Input natom must be defined");
  ptInstrng = strstr(buf, "The occurence of a tab");
  ptSize    = strstr(buf, "The size of your input file");
  if (ptError && !ptInvars0 && !ptInstrng && !ptSize)
    {
      _parse_error = 1; /* Abinit file with errors. */
      _parse_message = buf;
      pthread_exit(NULL);
    }
  else if (ptError && (ptInvars0 || ptInstrng || ptSize))
    {
      _parse_error = 2; /* Not an abinit file. */
      free(buf);
      pthread_exit(NULL);
    }
  else
    _parse_error = 0;
  free(buf);
}
void FC_FUNC(leave_new, LEAVE_NEW)()
{
  _parse_error = 3;
  pthread_exit(NULL);
}
void FC_FUNC(timab, TIMAB)(int *nn, int *option, double tottim[2])
{
}
void FC_FUNC(psp_from_data, PSP_FROM_DATA)()
{
}

/* A replacement function, for systems that lack strndup.

   Copyright (C) 1996, 1997, 1998, 2001, 2002, 2003, 2005, 2006, 2007, 2009,
   2010 Free Software Foundation, Inc.
*/

char *
my_strndup (char const *s, size_t n)
{
  size_t len = my_strnlen (s, n);
  char *new = malloc (len + 1);

  if (new == NULL)
    return NULL;

  new[len] = '\0';
  return memcpy (new, s, len);
}


/* A replacement function, for systems that lack strnlen.
   Copyright (C) 2005, 2006, 2007, 2009, 2010 Free Software Foundation, Inc.
   Written by Simon Josefsson.*/

/* Find the length of STRING, but scan at most MAXLEN characters.
   If no '\0' terminator is found in that many characters, return MAXLEN.  */

size_t
my_strnlen (const char *string, size_t maxlen)
{
  const char *end = memchr (string, '\0', maxlen);
  return end ? (size_t) (end - string) : maxlen;
}
