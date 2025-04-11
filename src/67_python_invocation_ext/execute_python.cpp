
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

/********************
 * Modified from the file execute_python.c of
 * Olivier Parcollet's project: execute_python.
 * See https://github.com/parcollet/execute_python
 * ******************/


#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>

#include <stdexcept>
#include <string>
#include <fstream>

#include <iostream>
#include <fstream>
#include <iomanip>

#if defined HAVE_MPI
using namespace std;
#include <mpi.h>

#include "execute_python.hpp"

//--------------------
// The function to run Python
//--------------------

#define AS_STRING(X) AS_STRING2(X)
#define AS_STRING2(X) #X

// loads the function F from the shared lib
#define LOAD(F, R, ...)                                 \
 R (*F)(__VA_ARGS__);                                   \
 *(void**)(&F) = dlsym(libpython_handle, AS_STRING(F)); \
 if ((error = dlerror()) != NULL) {                     \
  fprintf(stderr, "%s\n", error);                       \
  return 1;                                             \
 }

static void* libpython_handle = NULL;

//--------------------

extern "C" {
int init_python_interpreter(const char* python_so) {
 char* error;
 libpython_handle = dlopen(python_so, RTLD_GLOBAL | RTLD_LAZY);
 if (!libpython_handle) {
  fprintf(stderr, "Can not find Python !\n%s\n", dlerror());
  return 1;
 }
 dlerror();                     // clear any existing error
 LOAD(Py_Initialize, void, );   // loads the functions that we will need
 (*Py_Initialize)();            // initialize the interpreter
 return 0;
}
}

// //--------------------
//
// // int init_python_interpreter_from_env(const char* env_var) {
// //  char *python_so = getenv(env_var);
// //  if (!python_so) {
// //   fprintf(stderr, "Can not find the environment variable %s\n", env_var);
// //   return 1;
// //  }
// //  return init_python_interpreter(python_so);
// // }

//--------------------

extern "C" {
int execute_python_file(const char* filename) {
 char* error;
 if (!libpython_handle) {
  fprintf(stderr, "Python is not initialized. You forgot to call init_python_interpreter !\n");
  return 1;
 }

  // check Python is running
  LOAD(Py_IsInitialized, int, );
  if (!(*Py_IsInitialized)()) {
   fprintf(stderr, "Python Interpreter failed to initialize\n");
   return 1;
 }

 // open the script file, report error, and run it in the interpreter
 FILE* file = fopen(filename, "r");
 if (!file) {
  fprintf(stderr, "file %s not found \n", filename);
  return 1;
 }

 // LOAD(PyRun_SimpleString, int, const char *);
 LOAD(PyRun_SimpleFile, int, FILE*, const char*);
 int result = PyRun_SimpleFile(file, filename);
 if (result == -1) {
  fprintf(stderr, "return of PyRun_SimpleFile %d \n", result);
  return 1;
 }

 return 0;
}
}

//--------------------

extern "C" {
int close_python_interpreter() {
 fprintf(stdout, "In close_python_interpreter.");
 char* error;

 // close the interpreter
 LOAD(Py_Finalize, void, );
 (*Py_Finalize)();

 // close the shared lib
 dlclose(libpython_handle);

 return 0;
}
}

#endif
