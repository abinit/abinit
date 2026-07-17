As of now, GPU support using OpenMP offload has been
validated with the following compilers:
- NVHPC
- Cray

It seems you use a different compiler.
While ABINIT may compile successfuly, issues may
arise during at execution, depending on OpenMP support
level in your compiler.
