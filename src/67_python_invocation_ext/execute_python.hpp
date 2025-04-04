
// #if defined HAVE_CONFIG_H
// #include "config.h"
// #endif

/********************
 * Modified from the file execute_python.h of
 * Olivier Parcollet's project: execute_python.
 * See https://github.com/parcollet/execute_python
 * ******************/

// #if defined HAVE_MPI
// #include <mpi.h>

// starts the interpreter
// the location of the PythonLibrary must be in the var env_var

// int init_python_interpreter_from_env(const char* env_var);

extern "C" void test_python();

// states the interpreter
extern "C" int init_python_interpreter(const char* python_so);

// // executes code in the interpreter
extern "C" int execute_python_file(const char* filename);
//
// // call only ONCE at the very end !
extern "C" int close_python_interpreter();

// #endif
