/* dev_spec_hip.cpp*/

/*
 * Copyright (C) 2008-2024 ABINIT Group (MMancini,FDahm)
 * this file is distributed under the terms of the
 * gnu general public license, see ~abinit/COPYING
 * or http://www.gnu.org/copyleft/gpl.txt.
 *
 */

#include <stdio.h>
#include <abi_gpu_header_common.h>
#include <hip/hip_runtime_api.h>
#include "hip_api_error_check.h"

static int version_2_cores(int major, int minor);

/*=========================================================================*/
/*________________________ GPU_function called by HOST_____________________*/
/*=========================================================================*/
// display HIP device info
static void  prt_dev_info()
{
  int deviceCount;
  HIP_API_CHECK(hipGetDeviceCount(&deviceCount));
  for (int dev = 0; dev < deviceCount; ++dev)
    {
      hipDeviceProp_t deviceProp;
      HIP_API_CHECK(hipGetDeviceProperties(&deviceProp, dev));
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
  /* hip check errors */
  hipError_t hipError;
  hipError = hipGetLastError();
  if(hipError != hipSuccess)
    { fprintf(stderr, "HIP Runtime API Error reported : %s %d\n", hipGetErrorString(hipError),line);
      exit(EXIT_FAILURE);
    }
  return;
}


// Gives the number of GPU devices ---------
extern "C"
void get_gpu_ndev_(int* ndevice)
{
  int deviceCount;
  HIP_API_CHECK(hipGetDeviceCount(&deviceCount));
  *ndevice = deviceCount;

  return;
}

// Gives the max memory available for a GPU device ---------
extern "C"
void get_gpu_max_mem_(int* device, float* max_mem)
{
  hipDeviceProp_t deviceProp;
  HIP_API_CHECK(hipGetDeviceProperties(&deviceProp, *device));
  *max_mem = (float) deviceProp.totalGlobalMem;
   return;
}

// Gives currently free memory available for current GPU device ---------
extern "C"
void gpu_get_free_mem_cpp(size_t* free_mem)
{
   size_t max_mem;
   HIP_API_CHECK(hipMemGetInfo(&max_mem, free_mem));
   return;
}

// Set the device if it exists   -----------------
extern "C"
void set_dev_(int* gpudevice)
{
 if(*gpudevice >-1){
   hipError_t hipError;
   int deviceCount;
   HIP_API_CHECK(hipGetDeviceCount(&deviceCount));
   if(deviceCount>*gpudevice){
     HIP_API_CHECK(hipSetDevice(*gpudevice));
     hipError = hipGetLastError();
     if(hipError != hipSuccess){
       fprintf(stderr, "HIP Runtime API Error reported : %s\n", hipGetErrorString(hipError));
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
  HIP_API_CHECK(hipDeviceReset());
  return;
}

// Synchronize device (makes the CPU waits the GPU to finish all running kernels)
// this is required when using mamanged memory in order to reuse safely on CPU data
// that were processed / modified by the GPU
extern "C"
void gpu_device_synchronize_cpp()
{
  hipError_t hipError = hipDeviceSynchronize();
  if(hipError != hipSuccess)
    {
      fprintf(stderr, "HIP Runtime API Error reported : %s when trying to call hipDeviceSynchronize\n", hipGetErrorString(hipError));
      fflush(stderr);
    }
  return;
}

//
extern "C"
void gpu_get_device_cpp(int *deviceId)
{
  CHECK_HIP_ERROR( hipGetDevice(deviceId) );

  return;
}

//
extern "C"
void gpu_data_prefetch_async_cpp(const void* devPtr, size_t count, int deviceId)
{

  CHECK_HIP_ERROR( hipMemPrefetchAsync(devPtr, count, deviceId) );

  return;
}

//
extern "C"
void gpu_memory_advise_cpp(const void* devPtr, size_t count, hipMemoryAdvise advice, int deviceId)
{

  CHECK_HIP_ERROR( hipMemAdvise(devPtr, count, advice, deviceId) );

  return;
}

// Get context  -----------------------
extern "C"
void check_context_(int *res,char *message)
{
  *res=1;
  hipError_t state=hipFree(0);
  if (state!=hipSuccess){
    sprintf(message,"Unable to initialize a Cuda context: %s \n",hipGetErrorString(state));
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
  hipDeviceProp_t deviceProp;
  HIP_API_CHECK(hipGetDeviceProperties(&deviceProp, *device));
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
  HIP_API_CHECK(hipGetDeviceCount(&deviceCount));
  for (int idev = 0; idev < deviceCount; ++idev)
    {
      hipDeviceProp_t deviceProp;
      HIP_API_CHECK(hipGetDeviceProperties(&deviceProp, idev));
      //We check that no device is in "emu" mode
      if( deviceProp.major != 9999 ) {
#if defined HAVE_GPU_HIP_DP
      //We check that double precision is available, c.c. >= 1.3 )
      if( (deviceProp.major>1)||(deviceProp.minor>2) )
#endif
	*ndev+=1;
      }
    }
}


// Get number of cores of device  --------------
//This function is present in hip SDK: see ${HIPROOT}/common/inc/helper_hip_drvapi.h
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
        { 0x86, 128}, // Ampere Generation (SM 8.6)
        { 0x87, 128}, // Ampere Generation (SM 8.7)
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
  HIP_API_CHECK(hipMemGetInfo(&free,&total));
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

  if (hipMalloc(gpu_ptr,*size) != hipSuccess)
  {
    fprintf(stderr, "ERROR: alloc_on_gpu failed allocating %ld bytes :%s\n", *size, hipGetErrorString(hipGetLastError()));
    fflush(stderr);
    abi_cabort();
  }
}

/*============================================================================*/
/* Free memory location pointed by gpu_ptr                                    */
/* OUTPUT gpu_ptr= C_PTR on gpu memory location that has been allocated       */
/* WARNING! : this routine is a dummy one when HAVE_GPU_HIP is not enabled    */
/*            the correct one is in xx_gpu_toolbox/dev_spec.cu                */
/*============================================================================*/

extern "C" void dealloc_on_gpu_(void **gpu_ptr)
{
  if(*gpu_ptr==NULL)
    return;

  if (hipFree(*gpu_ptr) != hipSuccess)
  {
    fprintf(stderr, "ERROR: dealloc_on_gpu failed :%s\n",hipGetErrorString(hipGetLastError()));
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
/* WARNING! : this routine is a dummy one when HAVE_GPU_HIP is not enabled    */
/*            the correct one is in xx_gpu_toolbox/dev_spec.cu                */
/*============================================================================*/

extern "C" void copy_on_gpu_(void *cpu_ptr, void **gpu_ptr, const size_t* size)
{
  if (hipMemcpy(*gpu_ptr, cpu_ptr, *size, hipMemcpyHostToDevice) != hipSuccess)
  {
    fprintf(stderr, "ERROR: copy_on_gpu failed : %s\n",hipGetErrorString(hipGetLastError()));
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

extern "C" void copy_from_gpu_(void *cpu_ptr,void **gpu_ptr, const size_t* size)
{
  if (hipMemcpy(cpu_ptr, *gpu_ptr, *size, hipMemcpyDeviceToHost) != hipSuccess)
  {
    fprintf(stderr, "ERROR: copy_from_gpu failed : %s\n",hipGetErrorString(hipGetLastError()));
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
/* WARNING! : this routine is a dummy one when HAVE_GPU_HIP is not enabled    */
/*            the correct one is in xx_gpu_toolbox/dev_spec.cu                */
/*============================================================================*/

extern "C" void copy_gpu_to_gpu_cpp_(void **dest_gpu_ptr, void **src_gpu_ptr, const size_t* size)
{
  if (hipMemcpy(*dest_gpu_ptr, *src_gpu_ptr, *size, hipMemcpyDeviceToDevice) != hipSuccess)
  {
    fprintf(stderr, "ERROR: copy_gpu_to_gpu failed (dest=%p, src=%p, size=%ld): %s\n",
            *dest_gpu_ptr, *src_gpu_ptr, *size, hipGetErrorString(hipGetLastError()));
    fflush(stderr);
    abi_cabort();
  }
}

/*============================================================================*/
/* Reset array (just wrapping hipMemset)                                      */
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
  if(hipMemset(*gpu_ptr, *value, *size_in_bytes)!=hipSuccess){
    fprintf(stderr, "ERROR: gpu_memset at address %p failed : %s\n",*gpu_ptr,hipGetErrorString(hipGetLastError()));
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

  hipPointerAttribute_t attributes;

  CHECK_HIP_ERROR(hipPointerGetAttributes(&attributes, *gpu_ptr));

  if(attributes.devicePointer != NULL)
  {
    *is_allocated = true;
  }

} // gpu_allocated_impl_

/*============================================================================*/
/* Utility routine to print memory location of a hip managed pointer.         */
/*                                                                            */
/* We check that the pointer has actually been allocated with                 */
/* hipMallocManaged and then prints device and host addresses.                */
/*                                                                            */
/* INPUTS                                                                     */
/*  gpu_ptr = C_PTR on gpu memory location                                    */
/*                                                                            */
/* OUTPUT                                                                     */
/*  None.                                                                     */
/*============================================================================*/

extern "C" void gpu_managed_ptr_status_(void **gpu_ptr, const char* str)
{

  hipPointerAttribute_t attributes;

  CHECK_HIP_ERROR(hipPointerGetAttributes(&attributes, *gpu_ptr));

  if(attributes.type == hipMemoryTypeUnified)
  {
    printf("[%s] ptr %p is unified memory, host addr=%p, device addr=%p\n", str, *gpu_ptr,
           attributes.hostPointer,
           attributes.devicePointer);
    fflush(stdout);
  } else if(attributes.type == hipMemoryTypeDevice) {
    printf("[%s] ptr %p is a device ptr.\n", str, *gpu_ptr);
    fflush(stdout);
  } else {
    printf("[%s] ptr %p is neither a unified memory pointer nor a device pointer.\n", str, *gpu_ptr);
    fflush(stdout);
  }

} // gpu_managed_ptr_status_
