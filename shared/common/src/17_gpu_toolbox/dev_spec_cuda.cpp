/* dev_spec_cuda.cpp*/

/*
 * Copyright (C) 2008-2024 ABINIT Group (MMancini,FDahm)
 * this file is distributed under the terms of the
 * gnu general public license, see ~abinit/COPYING
 * or http://www.gnu.org/copyleft/gpl.txt.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <abi_gpu_header_common.h>
#include <cuda_runtime_api.h>
#include "cuda_api_error_check.h"

static int version_2_cores(int major, int minor);

/*=========================================================================*/
/*________________________ GPU_function called by HOST_____________________*/
/*=========================================================================*/
// display CUDA device info
static void prt_dev_info()
{
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  for (int dev = 0; dev < deviceCount; ++dev)
    {
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, dev);
      int NProcs=deviceProp.multiProcessorCount;
      int NCores=version_2_cores(deviceProp.major, deviceProp.minor);
      printf("\n___________________________________________________________________\n");
      printf(  "__________  Graphic Card Properties  ______________________________\n");
      printf("\n  Device %d: \"%s\"\n", dev, deviceProp.name);
      printf("  Revision number:                               %d.%d\n", deviceProp.major,deviceProp.minor);
      printf("  Total amount of global memory:                 %3.1f Mbytes\n", deviceProp.totalGlobalMem/1048576.);
      printf("  Clock rate:                                    %3.1f GHz\n", deviceProp.clockRate/1000000.);
      printf("  Number of processors/cores:                    %d/%d\n", NProcs,NCores);
      if (NCores<0) {
        printf("  Max GFLOPS:                                    undefined (add new def. in version_2_cores function)\n");
      } else {
        printf("  Max GFLOPS:                                    %d GFP\n", NCores*deviceProp.multiProcessorCount * deviceProp.clockRate/1000000);
      }
      printf("  Total amount of constant memory:               %d bytes\n",(int) deviceProp.totalConstMem);
      printf("  Total amount of shared memory per block:       %d bytes\n",(int) deviceProp.sharedMemPerBlock);
      printf("  Total number of registers available per block: %d\n", deviceProp.regsPerBlock);
      printf("___________________________________________________________________\n");
      fflush(stdout);
      if( (int) deviceProp.totalConstMem<0) break;
      //if(deviceProp.major==9999){printf("EXIT: PROBLEM WITH AVAILABLE DEVICES \n");exit(0);}
    }
}


// Explicit Cuda Error ---------------------
void check_err(int line )
{
  /* cuda check errors */
  cudaError_t cudaError;
  cudaError = cudaGetLastError();
  if(cudaError != cudaSuccess)
    {
      fprintf(stderr, "CUDA Runtime API Error reported : %s %d\n", cudaGetErrorString(cudaError),line);
      exit(EXIT_FAILURE);
    }
  return;
}


// Gives the number of GPU devices ---------
extern "C"
void get_gpu_ndev_(int* ndevice)
{
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  *ndevice = deviceCount;

  return;
}

// Gives the max memory available for a GPU device ---------
extern "C"
void get_gpu_max_mem_(int* device, float* max_mem)
{
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, *device);
  *max_mem = (float) deviceProp.totalGlobalMem;
   return;
}

// Gives currently free memory available for current GPU device ---------
extern "C"
void gpu_get_free_mem_cpp(size_t* free_mem)
{
   size_t max_mem;
   cudaMemGetInfo(&max_mem, free_mem);
   return;
}

// Set the device if it exists   -----------------
extern "C"
void set_dev_(int* gpudevice)
{
 if(*gpudevice >-1){
   cudaError_t cudaError;
   int deviceCount;
   cudaGetDeviceCount(&deviceCount);
   if(deviceCount>*gpudevice){
     cudaSetDevice(*gpudevice);
     cudaError = cudaGetLastError();
     if(cudaError != cudaSuccess){
       fprintf(stderr, "CUDA Runtime API Error reported : %s\n", cudaGetErrorString(cudaError));
       fflush(stderr);
       exit(1);
     }
   }
   else *gpudevice=-1;
 }
  return;
}


// Unset the devices  -----------------
extern "C"
void unset_dev_()
{
#if defined HAVE_GPU_CUDA3
  cudaThreadExit();
#else
  cudaDeviceReset();
#endif
  return;
}

// Synchronize device (makes the CPU waits the GPU to finish all running kernels)
// this is required when using mamanged memory in order to reuse safely on CPU data
// that were processed / modified by the GPU
extern "C"
void gpu_device_synchronize_cpp()
{
  cudaError_t cudaError = cudaDeviceSynchronize();
  if(cudaError != cudaSuccess)
    {
      fprintf(stderr, "CUDA Runtime API Error reported : %s when trying to call cudaDeviceSynchronize\n", cudaGetErrorString(cudaError));
      fflush(stderr);
    }
  return;
}

//
extern "C"
void gpu_get_device_cpp(int *deviceId)
{
  CHECK_CUDA_ERROR( cudaGetDevice(deviceId) );

  return;
}

//
extern "C"
void gpu_data_prefetch_async_cpp(const void* devPtr, size_t count, int deviceId)
{

  CHECK_CUDA_ERROR( cudaMemPrefetchAsync(devPtr, count, deviceId) );

  return;
}

//
extern "C"
void gpu_memory_advise_cpp(const void* devPtr, size_t count, cudaMemoryAdvise advice, int deviceId)
{

  CHECK_CUDA_ERROR( cudaMemAdvise(devPtr, count, advice, deviceId) );

  return;
}

// Get context  -----------------------
extern "C"
void check_context_(int *res,char *message)
{
  *res=1;
  cudaError_t state=cudaFree(0);
  if (state!=cudaSuccess){
    sprintf(message,"Unable to initialize a Cuda context: %s \n",cudaGetErrorString(state));
    *res=0;
    unset_dev_();
  }
}


// Get info from device  --------------
extern "C"
void  get_dev_info_(int* device,
		    char* name,
		    int* lenname,
		    int vers[2],
		    float* globalmem,
		    float* clockrate,
		    int* gflops,
		    int* constmem,
		    int* sharemem,
		    int* regist,
		    int* nprocs,
		    int* ncores
		    )
{
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, *device);
  strcpy(name,deviceProp.name);
  *lenname = strlen( name );
  vers[0] = deviceProp.major;
  vers[1] = deviceProp.minor;
  *globalmem = deviceProp.totalGlobalMem/1048576.;
  *clockrate = deviceProp.clockRate/1000000.;
  *nprocs = deviceProp.multiProcessorCount;
  *ncores = version_2_cores(deviceProp.major,deviceProp.minor);
  *gflops = int(deviceProp.multiProcessorCount*version_2_cores(deviceProp.major,deviceProp.minor)*(deviceProp.clockRate/1000000.));
  *constmem = deviceProp.totalConstMem;
  *sharemem =  deviceProp.sharedMemPerBlock;
  *regist = deviceProp.regsPerBlock;
}


// Get number of devices  --------------
extern "C"
void c_get_ndevice_(int* ndev)
{
  *ndev=0;
  int deviceCount;
  cudaGetDeviceCount(&deviceCount);
  for (int idev = 0; idev < deviceCount; ++idev)
    {
      cudaDeviceProp deviceProp;
      cudaGetDeviceProperties(&deviceProp, idev);
      //We check that no device is in "emu" mode
      if( deviceProp.major != 9999 ) {
#if defined HAVE_GPU_CUDA_DP
      //We check that double precision is available, c.c. >= 1.3 )
      if( (deviceProp.major>1)||(deviceProp.minor>2) )
#endif
	*ndev+=1;
      }
    }
}


// Get number of cores of device  --------------
//This function is present in cuda SDK: see ${CUDAROOT}/common/inc/helper_cuda_drvapi.h
//To be completed for new card versions
static
int version_2_cores(int major, int minor)
{
    // Defines for GPU Architecture types (using the SM version to determine the # of cores per SM
    typedef struct
    {
        int SM; // 0xMm (hexidecimal notation), M = SM Major version, and m = SM minor version
        int Cores;
    } sSMtoCores;
    sSMtoCores nGpuArchCoresPerSM[] =
    {
        { 0x10,  8 }, // Tesla Generation (SM 1.0) G80 class
        { 0x11,  8 }, // Tesla Generation (SM 1.1) G8x class
        { 0x12,  8 }, // Tesla Generation (SM 1.2) G9x class
        { 0x13,  8 }, // Tesla Generation (SM 1.3) GT200 class
        { 0x20, 32 }, // Fermi Generation (SM 2.0) GF100 class
        { 0x21, 48 }, // Fermi Generation (SM 2.1) GF10x class
        { 0x30, 192}, // Kepler Generation (SM 3.0) GK10x class
        { 0x32, 192}, // Kepler Generation (SM 3.2) GK10x class
        { 0x35, 192}, // Kepler Generation (SM 3.5) GK11x class
        { 0x50, 128}, // Maxwell Generation (SM 5.0) GM10x class
        { 0x60, 64 }, // Pascal Generation (SM 6.0) GP100 class
        { 0x61, 128}, // Pascal Generation (SM 6.1) GP10x class
        { 0x62, 128}, // Pascal Generation (SM 6.2) GP10x class
        { 0x70, 64 }, // Volta Generation (SM 7.0) GV100 class
        { 0x72, 64 }, // Volta Generation (SM 7.2) AGX class
        { 0x75, 64 }, // Turing Generation (SM 7.5) RTX class
        { 0x80, 64 }, // Ampere Generation (SM 8.0) A100 class
        { 0x86, 128}, // Ampere Generation (SM 8.6) RTX class
        { 0x87, 128}, // Ampere Generation (SM 8.7) AGX class
        { 0x89, 128}, // Ada Lovelace Generation (SM 8.9) RTX class
        { 0x90, 128}, // Hooper Generation (SM 9.0) H100 class
        {   -1, -1 }
    };
    int index = 0;
    while (nGpuArchCoresPerSM[index].SM != -1)
    {
        if (nGpuArchCoresPerSM[index].SM == ((major << 4) + minor))
        {
            return nGpuArchCoresPerSM[index].Cores;
        }
        index++;
    }

//  printf("MapSMtoCores for SM %d.%d is undefined.  Default to use %d Cores/SM\n", major, minor, nGpuArchCoresPerSM[7].Cores);
    return nGpuArchCoresPerSM[10].Cores;
}


/***************************************************************/
/*******                                                ********/
/*******      GPU MEMORY MANAGEMENT ROUTINES            ********/
/*******                                                ********/
/***************************************************************/

/*============================================================================*/
/* Print memory information (total amount and free available)                 */
/*============================================================================*/

extern "C" void check_gpu_mem_(const char* str)
{
  size_t free,total;
  cudaMemGetInfo(&free,&total);
  printf("[%s] *** GPU memory : Occupied => %4.2fMo   | Free =>  %4.2fMo   | Total =>  %4.2fMo ***\n",
         str, (total-free)*1e-6, free*1e-6, total*1e-6);
  fflush(stdout);
}

/*============================================================================*/
/* Allocate size byte in gpu memory and returns in gpu_ptr this location      */
/* INPUTS size= size in byte to allocate                                      */
/* OUTPUT gpu_ptr= C_PTR on gpu memory location that has been allocated       */
/*============================================================================*/

extern "C" void alloc_on_gpu_(void **gpu_ptr, const size_t* size)
{

  //check_gpu_mem_("alloc_on_gpu_");

  if (cudaMalloc(gpu_ptr,*size) != cudaSuccess)
  {
    fprintf(stderr, "ERROR: alloc_on_gpu failed allocating %ld bytes :%s\n", *size, cudaGetErrorString(cudaGetLastError()));
    fflush(stderr);
    abi_cabort();
  }
}

/*============================================================================*/
/* Free memory location pointed by gpu_ptr                                    */
/* OUTPUT gpu_ptr= C_PTR on gpu memory location that has been allocated       */
/* WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled   */
/*            the correct one is in xx_gpu_toolbox/dev_spec.cu                */
/*============================================================================*/

extern "C" void dealloc_on_gpu_(void **gpu_ptr)
{
  if(*gpu_ptr==NULL)
    return;

  if (cudaFree(*gpu_ptr) != cudaSuccess)
  {
    fprintf(stderr, "ERROR: dealloc_on_gpu failed :%s\n",cudaGetErrorString(cudaGetLastError()));
    fflush(stderr);
    abi_cabort();
  }
  *gpu_ptr=NULL;
}

/*============================================================================*/
/* Copy size byte from cpu pointer to gpu pointer.                            */
/* INPUTS                                                                     */
/*  size = size in byte to copy                                               */
/*  cpu_ptr = host memory location (LOC)                                      */
/* OUTPUT                                                                     */
/*  gpu_ptr = C_PTR : gpu memory location                                     */
/* WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled   */
/*            the correct one is in xx_gpu_toolbox/dev_spec.cu                */
/*============================================================================*/

extern "C" void copy_on_gpu_(void *cpu_ptr, void **gpu_ptr, const size_t* size)
{
  if (cudaMemcpy(*gpu_ptr, cpu_ptr, *size, cudaMemcpyHostToDevice) != cudaSuccess)
  {
    fprintf(stderr, "ERROR: copy_on_gpu failed : %s\n",cudaGetErrorString(cudaGetLastError()));
    fflush(stderr);
    abi_cabort();
  }
}

/*============================================================================*/
/* Copy size byte from gpu pointer to cpu pointer.                            */
/* INPUTS                                                                     */
/*  size = size in byte to copy                                               */
/*  gpu_ptr = C_PTR : gpu memory location                                     */
/* OUTPUT                                                                     */
/*  cpu_ptr = host memory location (LOC of an allocated array)                */
/*============================================================================*/

extern "C" void copy_from_gpu_(void *cpu_ptr, void **gpu_ptr, const size_t* size)
{
  if (cudaMemcpy(cpu_ptr, *gpu_ptr, *size, cudaMemcpyDeviceToHost) != cudaSuccess)
  {
    fprintf(stderr, "ERROR: copy_from_gpu failed : %s\n",cudaGetErrorString(cudaGetLastError()));
    fflush(stderr);
    abi_cabort();
  }
}

/*============================================================================*/
/* Copy size byte from gpu to gpu memory.                                     */
/* INPUTS                                                                     */
/*  size = size in byte to copy                                               */
/*  src_gpu_ptr                                                               */
/* OUTPUT                                                                     */
/*  dest_gpu_ptr = C_PTR on gpu memory location                               */
/* WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled   */
/*            the correct one is in xx_gpu_toolbox/dev_spec.cu                */
/*============================================================================*/

extern "C" void copy_gpu_to_gpu_cpp_(void **dest_gpu_ptr, void **src_gpu_ptr, const size_t* size)
{
  if (cudaMemcpy(*dest_gpu_ptr, *src_gpu_ptr, *size, cudaMemcpyDeviceToDevice) != cudaSuccess)
  {
    fprintf(stderr, "ERROR: copy_gpu_to_gpu failed (dest=%p, src=%p, size=%ld): %s\n",
            *dest_gpu_ptr, *src_gpu_ptr, *size, cudaGetErrorString(cudaGetLastError()));
    fflush(stderr);
    abi_cabort();
  }
}

/*============================================================================*/
/* Reset array (just wrapping cudaMemset)                                     */
/*                                                                            */
/* INPUTS                                                                     */
/*  gpu_ptr = C_PTR on gpu memory location                                    */
/*  value = integer used to initialize each bytes (should be in range [0,255])*/
/*  size = size in bytes of the region to be set                              */
/*                                                                            */
/* OUTPUT                                                                     */
/*  None                                                                      */
/*============================================================================*/

extern "C" void gpu_memset_cpp_(void **gpu_ptr, const int32_t* value, const size_t* size_in_bytes)
{
  if(cudaMemset(*gpu_ptr, *value, *size_in_bytes)!=cudaSuccess){
    fprintf(stderr, "ERROR: gpu_memset at address %p failed : %s\n",*gpu_ptr,cudaGetErrorString(cudaGetLastError()));
    fflush(stderr);
    abi_cabort();
  }
}

/*============================================================================*/
/* Kind of equivalent of fortran "allocated". Check if a gpu pointer          */
/* actually points to device allocated memory.                                */
/*                                                                            */
/* This is void function because I can't manage to bind it via iso_c_binding  */
/* as a fortran function; binding as a subroutine is ok though (?!)           */
/*                                                                            */
/* INPUTS                                                                     */
/*  gpu_ptr = C_PTR on gpu memory location                                    */
/*                                                                            */
/* OUTPUT                                                                     */
/*  boolean/logical (false = not allocated, true = allocated)                 */
/*============================================================================*/

extern "C" void gpu_allocated_impl_(void **gpu_ptr, bool* is_allocated)
{

  *is_allocated = false;

  cudaPointerAttributes attributes;

  CHECK_CUDA_ERROR(cudaPointerGetAttributes(&attributes, *gpu_ptr));

  if(attributes.devicePointer != NULL)
  {
    *is_allocated = true;
  }

} // gpu_allocated_impl_

/*============================================================================*/
/* Utility routine to print memory location of a cuda managed pointer.        */
/*                                                                            */
/* We check that the pointer has actually been allocated with                 */
/* cudaMallocManaged and then prints device and host addresses.               */
/*                                                                            */
/* INPUTS                                                                     */
/*  gpu_ptr = C_PTR on gpu memory location                                    */
/*                                                                            */
/* OUTPUT                                                                     */
/*  None.                                                                     */
/*============================================================================*/

extern "C" void gpu_managed_ptr_status_(void **gpu_ptr, const char* str)
{

  cudaPointerAttributes attributes;

  CHECK_CUDA_ERROR(cudaPointerGetAttributes(&attributes, *gpu_ptr));

  if(attributes.type == cudaMemoryTypeManaged)
  {
    printf("[%s] ptr %p is memory managed, host addr=%p, device addr=%p\n", str, *gpu_ptr,
           attributes.hostPointer,
           attributes.devicePointer);
    fflush(stdout);
  } else if(attributes.type == cudaMemoryTypeDevice) {
    printf("[%s] ptr %p is a device ptr.\n", str, *gpu_ptr);
    fflush(stdout);
  } else {
    printf("[%s] ptr %p is neither a memory managed pointer nor a device pointer.\n", str, *gpu_ptr);
    fflush(stdout);
  }

} // gpu_managed_ptr_status_
