
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

/********************
 * Modified from the file execute_python.h of
 * Olivier Parcollet's project: execute_python.
 * See https://github.com/parcollet/execute_python
 * ******************/

#if defined HAVE_MPI
// starts the interpreter
// the location of the PythonLibrary must be in the var env_var
int init_python_interpreter_from_env(const char* env_var);

// states the interpreter
int init_python_interpreter(const char* python_so);

// executes code in the interpreter
int execute_python_file(const char* filename);

// call only ONCE at the very end !
int close_python_interpreter();

extern "C"{
    void invoke_python_triqs(int rank, char* filapp_in);
}
#else
extern "C"{
    void invoke_python_triqs(int rank, char* filapp_in);
}
#endif
