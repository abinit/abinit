#include <mpi.h>
#include <netcdf.h>
#include <stdio.h>


int main(int argc, char* argv[])
{

  MPI_Init(NULL,NULL);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Info info = MPI_INFO_NULL;

  int ncid;

  int ierr = nc_create_par("conftest.nc", NC_NETCDF4, comm, info, &ncid);
  //int ierr = nc_create_par("conftest.nc", NC_CLASSIC_MODEL, comm, info, &ncid);
  printf("  nc_create_par : ierr= %d\n",ierr);
  MPI_Finalize();

  if (ierr == NC_ENOPAR)
    printf("Check how your netcdf library was build; apparenly it was not built with parallel I/O features\n");

  if(ierr != 0) return 1;

}
