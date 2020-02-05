/********************************************
 *
 * finiteElementHeatPropagator.c, in the C programming language, written for MATLAB MEX function generation
 * C script for heat propagation based on Monte Carlo input
 *
 * Copyright 2017, 2018 by Anders K. Hansen, DTU Fotonik
 *
 * This file is part of MCmatlab.
 *
 * MCmatlab is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MCmatlab is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MCmatlab.  If not, see <https://www.gnu.org/licenses/>.
 *
 ** COMPILING ON WINDOWS
 * Can be compiled in MATLAB using the MinGW-w64 compiler (GCC) with
 * "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall -pedantic' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall -pedantic' -outdir helperfuncs\private .\src\finiteElementHeatPropagator.c ".\src\libut.lib""
 * ... or the Microsoft Visual C++ compiler (MSVC) with
 * "mex COMPFLAGS='/Zp8 /GR /EHs /nologo /MD /openmp /W4 /WX /wd4204 /wd4100' -outdir helperfuncs\private .\src\finiteElementHeatPropagator.c ".\src\libut.lib""
 * In my experience, GCC produces faster machine code than MSVC.
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Go to MATLAB's addon manager and tell it to install the "Support for MinGW-w64 compiler"
 * 2. Type "mex -setup" in the MATLAB command window and ensure that MATLAB has set the C compiler to MinGW64
 * 3. mex should now be able to compile the code using the above GCC command
 *
 * The source code in this file is written is such a way that it is
 * compilable by either C or C++ compilers, either with GCC, MSVC or
 * the Nvidia CUDA compiler called NVCC, which is based on MSVC. To
 * compile with CUDA GPU acceleration support, you must have MSVC
 * installed. As of January 2020, mexcuda does not work with MSVC 2019,
 * so I'd recommend MSVC 2017. You also need the Parallel Computing
 * Toolbox, which you will find in the MATLAB addon manager. To compile, run:
 * "copyfile ./src/finiteElementHeatPropagator.c ./src/finiteElementHeatPropagator_CUDA.cu; mexcuda -llibut COMPFLAGS='-use_fast_math -res-usage $COMPFLAGS' -outdir helperfuncs\private .\src\finiteElementHeatPropagator_CUDA.cu ".\src\nvml.lib""
 *
 ** COMPILING ON MAC
 * As of June 2017, the macOS compiler doesn't support libut (for ctrl+c 
 * breaking) or openmp (for multithreading).
 * This file can then be compiled with "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -std=c11 -Wall -pedantic' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -std=c11 -Wall -pedantic' -outdir helperfuncs/private ./src/finiteElementHeatPropagator.c"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Install XCode from the App Store
 * 2. Type "mex -setup" in the MATLAB command window
 ********************************************/
// printf("Reached line %d...\n",__LINE__);mexEvalString("drawnow;");mexEvalString("drawnow;");mexEvalString("drawnow;"); // For inserting into code for debugging purposes

#include <math.h>
#include "mex.h"
#ifdef _OPENMP
  #include "omp.h"
#endif
#ifndef __clang__
  #ifdef __cplusplus
  extern "C"
  #endif
  extern bool utIsInterruptPending(); // Allows catching ctrl+c while executing the mex function
#endif
#ifdef __GNUC__ // This is defined for GCC and CLANG but not for Microsoft Visual C++ compiler
  #include <time.h>
#else
  #include <windows.h>
#endif

#ifdef __NVCC__
  #include <nvml.h>
  #define TILE_DIM 32
#endif

#define R 8.3144598f // Gas constant [J/mol/K]
#define CELSIUSZERO 273.15f

struct debug {
  double             dbls[3];
  unsigned long long ulls[3];
};

struct parameters {
  long nx;
  long ny;
  long nz;
  long nM;
  float *Tfinal;
  float *T1;
  float *T2;
  float *Tyxz;
  long nt;
  float dt;
  unsigned char *M;
  unsigned char *Myxz;
  float *c;
  float *b;
  float *dTdt_abs;
  float *logA;
  float *E;
  bool lightsOn;
  long heatSinked;
  float *Omega;
  bool calcDamage;
  double *tempSensorCornerIdxs;
  long nSensors;
  double *tempSensorInterpWeights;
  double *sensorTemps;
  float *maxMediaTemps;
};

long long getMicroSeconds() {
  #ifdef __GNUC__
  struct timespec time; clock_gettime(CLOCK_MONOTONIC, &time);
  return time.tv_sec*1000000 + time.tv_nsec/1000;
  #else
  static LARGE_INTEGER freq; QueryPerformanceFrequency(&freq);
  LARGE_INTEGER pfmcntr; QueryPerformanceCounter(&pfmcntr);
  return pfmcntr.QuadPart*1000000/freq.QuadPart;
  #endif
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void atomicMaxWrapperFLT(float *ptr, float val) { // Atomic maximum wrapper for single-precision floats
  if(*ptr < val) {
    #ifdef __NVCC__ // If compiling for CUDA
    int old = *(int *)ptr, oldifnorace;

    do {
      oldifnorace = old;
      old = atomicCAS((int *)ptr, oldifnorace,__float_as_int(val));
    } while (oldifnorace != old && old < __float_as_int(val));
    #else
    #ifdef _OPENMP
    #pragma omp critical
    if(*ptr < val)
    #endif
      *ptr = val;
    #endif
  }
}

#ifdef __NVCC__ // If compiling for CUDA
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line) {
   if (code != cudaSuccess) {
      printf("GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      mexEvalString("drawnow;");
      while(true) {;}
   }
}

void createDeviceStructs(struct parameters *P, struct parameters **P_devptr,
                         struct debug *D, struct debug **D_devptr) {
  long L = P->nx*P->ny*P->nz;
  struct parameters P_tempvar = *P;

  gpuErrchk(cudaMalloc(&P_tempvar.T1,L*sizeof(float)));
  gpuErrchk(cudaMemcpy(P_tempvar.T1,P->T1,L*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.T2,L*sizeof(float)));
  gpuErrchk(cudaMalloc(&P_tempvar.Tyxz,L*sizeof(float)));
  gpuErrchk(cudaMalloc(&P_tempvar.M,L*sizeof(unsigned char)));
  gpuErrchk(cudaMemcpy(P_tempvar.M,P->M,L*sizeof(unsigned char),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.Myxz,L*sizeof(unsigned char)));
  unsigned char *Myxztemp = (unsigned char *)malloc(L*sizeof(unsigned char));
  for(long ix=0;ix<P->nx;ix++) {
    for(long iy=0;iy<P->ny;iy++) {
      for(long iz=0;iz<P->nz;iz++) {
        Myxztemp[iy + ix*P->ny + iz*P->ny*P->nx] = P->M[ix + iy*P->nx + iz*P->nx*P->ny];
      }
    }
  }
  gpuErrchk(cudaMemcpy(P_tempvar.Myxz,Myxztemp,L*sizeof(unsigned char),cudaMemcpyHostToDevice));
  free(Myxztemp);
  gpuErrchk(cudaMalloc(&P_tempvar.c,P->nM*P->nM*3*sizeof(float)));
  gpuErrchk(cudaMemcpy(P_tempvar.c,P->c,P->nM*P->nM*3*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.b,L*sizeof(float)));
  gpuErrchk(cudaMalloc(&P_tempvar.dTdt_abs,L*sizeof(float)));
  gpuErrchk(cudaMemcpy(P_tempvar.dTdt_abs,P->dTdt_abs,L*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.logA,P->nM*sizeof(double)));
  gpuErrchk(cudaMemcpy(P_tempvar.logA,P->logA,P->nM*sizeof(double),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.E,P->nM*sizeof(double)));
  gpuErrchk(cudaMemcpy(P_tempvar.E,P->E,P->nM*sizeof(double),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.Omega,(P->calcDamage? L: 1)*sizeof(float)));
  gpuErrchk(cudaMemcpy(P_tempvar.Omega,P->Omega,(P->calcDamage? L: 1)*sizeof(float),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.tempSensorCornerIdxs,P->nSensors*sizeof(double)));
  gpuErrchk(cudaMemcpy(P_tempvar.tempSensorCornerIdxs,P->tempSensorCornerIdxs,P->nSensors*sizeof(double),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.tempSensorInterpWeights,P->nSensors*3*sizeof(double)));
  gpuErrchk(cudaMemcpy(P_tempvar.tempSensorInterpWeights,P->tempSensorInterpWeights,P->nSensors*3*sizeof(double),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&P_tempvar.sensorTemps,P->nSensors*(P->nt+1)*sizeof(double)));
  gpuErrchk(cudaMalloc(&P_tempvar.maxMediaTemps,P->nM*sizeof(float)));
  gpuErrchk(cudaMemcpy(P_tempvar.maxMediaTemps,P->maxMediaTemps,P->nM*sizeof(float),cudaMemcpyHostToDevice));

  gpuErrchk(cudaMalloc(P_devptr, sizeof(struct parameters)));
  gpuErrchk(cudaMemcpy(*P_devptr,&P_tempvar,sizeof(struct parameters),cudaMemcpyHostToDevice));

  // Allocate and copy debug struct
  struct debug D_tempvar = *D;
  gpuErrchk(cudaMalloc(D_devptr, sizeof(struct debug)));
  gpuErrchk(cudaMemcpy(*D_devptr,&D_tempvar,sizeof(struct debug),cudaMemcpyHostToDevice));
}

void retrieveAndFreeDeviceStructs(struct parameters *P, struct parameters *P_dev,
                                  struct debug *D, struct debug *D_dev) {
  long L = P->nx*P->ny*P->nz;
  struct parameters P_temp; gpuErrchk(cudaMemcpy(&P_temp,P_dev,sizeof(struct parameters),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(P->Tfinal,P_temp.T1,L*sizeof(float),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(P->Omega,P_temp.Omega,(P->calcDamage? L: 1)*sizeof(float),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(P->sensorTemps,P_temp.sensorTemps,P->nSensors*(P->nt+1)*sizeof(double),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaMemcpy(P->maxMediaTemps,P_temp.maxMediaTemps,P->nM*sizeof(float),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaFree(P_temp.T1));
  gpuErrchk(cudaFree(P_temp.T2));
  gpuErrchk(cudaFree(P_temp.Tyxz));
  gpuErrchk(cudaFree(P_temp.M));
  gpuErrchk(cudaFree(P_temp.Myxz));
  gpuErrchk(cudaFree(P_temp.c));
  gpuErrchk(cudaFree(P_temp.b));
  gpuErrchk(cudaFree(P_temp.dTdt_abs));
  gpuErrchk(cudaFree(P_temp.logA));
  gpuErrchk(cudaFree(P_temp.E));
  gpuErrchk(cudaFree(P_temp.Omega));
  gpuErrchk(cudaFree(P_temp.tempSensorCornerIdxs));
  gpuErrchk(cudaFree(P_temp.tempSensorInterpWeights));
  gpuErrchk(cudaFree(P_temp.sensorTemps));
  gpuErrchk(cudaFree(P_temp.maxMediaTemps));
  gpuErrchk(cudaFree(P_dev));

  gpuErrchk(cudaMemcpy(D, D_dev, sizeof(struct debug),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaFree(D_dev));
}
#endif

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void substep1a(struct parameters *P_global) {
  // Explicit part of substep 1 out of 3
  #ifdef __NVCC__
  __shared__ struct parameters P_var;
  struct parameters *P = &P_var;
  if(!threadIdx.x) P_var = *P_global; // Only let one thread per block do the copying. 
  __syncthreads(); // All threads in the block wait for the copy to have finished
  __shared__ float tile[TILE_DIM][TILE_DIM+1]; // +1 is to avoid memory bank conflicts
  
  long xTiles = (P->nx + TILE_DIM - 1)/TILE_DIM;
  long yTiles = (P->ny + TILE_DIM - 1)/TILE_DIM;
  for(long tileNum=blockIdx.x; tileNum<xTiles*yTiles*P->nz; tileNum += gridDim.x) {
    long tilexoffset = TILE_DIM*(tileNum%xTiles);
    long tileyoffset = TILE_DIM*(tileNum/xTiles%yTiles);
    long iz          =           tileNum/xTiles/yTiles;
    long ix = tilexoffset + threadIdx.x;
    long iy = tileyoffset + threadIdx.y;

    if(ix<P->nx && iy<P->ny) {
      long i = ix + iy*P->nx + iz*P->nx*P->ny;
      float dTdt = 0;
      if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) { // If not constant-temperature boundary voxel
        dTdt = P->lightsOn? P->dTdt_abs[i]: 0; // This is where the heat is deposited
        if(ix != 0      ) dTdt += (P->T1[i-1]           - P->T1[i])*P->c[P->M[i] + P->M[i-1          ]*P->nM                ]/2; // Factor ½ is because the first substep splits the x heat flow over an explicit and an implicit part
        if(ix != P->nx-1) dTdt += (P->T1[i+1]           - P->T1[i])*P->c[P->M[i] + P->M[i+1          ]*P->nM                ]/2; // Factor ½ is because the first substep splits the x heat flow over an explicit and an implicit part
        if(iy != 0      ) dTdt += (P->T1[i-P->nx]       - P->T1[i])*P->c[P->M[i] + P->M[i-P->nx      ]*P->nM +   P->nM*P->nM];
        if(iy != P->ny-1) dTdt += (P->T1[i+P->nx]       - P->T1[i])*P->c[P->M[i] + P->M[i+P->nx      ]*P->nM +   P->nM*P->nM];
        if(iz != 0      ) dTdt += (P->T1[i-P->nx*P->ny] - P->T1[i])*P->c[P->M[i] + P->M[i-P->nx*P->ny]*P->nM + 2*P->nM*P->nM];
        if(iz != P->nz-1) dTdt += (P->T1[i+P->nx*P->ny] - P->T1[i])*P->c[P->M[i] + P->M[i+P->nx*P->ny]*P->nM + 2*P->nM*P->nM];
      }
      tile[threadIdx.y][threadIdx.x] = P->T1[i] + P->dt*dTdt;
    }
    __syncthreads();
    // Save permuted (transposed) xyz -> yxz
    ix = tilexoffset + threadIdx.y;
    iy = tileyoffset + threadIdx.x;
    if(ix<P->nx && iy<P->ny) P->Tyxz[iz*P->ny*P->nx + ix*P->ny + iy] = tile[threadIdx.x][threadIdx.y];
    __syncthreads();
  }
  #else
  long ix,iy,iz;
  struct parameters *P = P_global;
  #ifdef _OPENMP
  #pragma omp for schedule(dynamic)
  #endif
  for(iz=0; iz<P->nz; iz++) {
    for(iy=0; iy<P->ny; iy++) {
      for(ix=0; ix<P->nx; ix++) {
        long i = ix + iy*P->nx + iz*P->nx*P->ny;
        float dTdt = 0;
        if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) { // If not constant-temperature boundary voxel
          dTdt = P->lightsOn? P->dTdt_abs[i]: 0; // This is where the heat is deposited
          if(ix != 0      ) dTdt += (P->T1[i-1]           - P->T1[i])*P->c[P->M[i] + P->M[i-1          ]*P->nM                ]/2; // Factor ½ is because the first substep splits the x heat flow over an explicit and an implicit part
          if(ix != P->nx-1) dTdt += (P->T1[i+1]           - P->T1[i])*P->c[P->M[i] + P->M[i+1          ]*P->nM                ]/2; // Factor ½ is because the first substep splits the x heat flow over an explicit and an implicit part
          if(iy != 0      ) dTdt += (P->T1[i-P->nx]       - P->T1[i])*P->c[P->M[i] + P->M[i-P->nx      ]*P->nM +   P->nM*P->nM];
          if(iy != P->ny-1) dTdt += (P->T1[i+P->nx]       - P->T1[i])*P->c[P->M[i] + P->M[i+P->nx      ]*P->nM +   P->nM*P->nM];
          if(iz != 0      ) dTdt += (P->T1[i-P->nx*P->ny] - P->T1[i])*P->c[P->M[i] + P->M[i-P->nx*P->ny]*P->nM + 2*P->nM*P->nM];
          if(iz != P->nz-1) dTdt += (P->T1[i+P->nx*P->ny] - P->T1[i])*P->c[P->M[i] + P->M[i+P->nx*P->ny]*P->nM + 2*P->nM*P->nM];
        }
        P->T2[i] = P->T1[i] + P->dt*dTdt;
      }
    }
  }
  #endif
}

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void substep1b(struct parameters *P_global) {
  // Implicit part of substep 1 out of 3
  #ifdef __NVCC__
  long threadNum = threadIdx.x + threadIdx.y*blockDim.x + blockIdx.x*blockDim.x*blockDim.y;
  __shared__ struct parameters P_var;
  struct parameters *P = &P_var;
  if(!threadIdx.x) P_var = *P_global; // Only let one thread per block do the copying.
  __syncthreads(); // All threads in the block wait for the copy to have finished
  for(long iyz=threadNum;iyz<P->ny*P->nz;iyz+=gridDim.x*TILE_DIM*TILE_DIM){
    long iy = iyz%P->ny;
    long iz = iyz/P->ny;
    for(long ix=0; ix<P->nx; ix++) {
      long i = iy + ix*P->ny + iz*P->ny*P->nx;
    // Thomson algorithm, sweeps up from 0 to nx-1 and then down from nx-1 to 0:
    // Algorithm is taken from https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
      P->b[i] = 1;
      if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) {
        if(ix < P->nx-1) P->b[i] += P->dt/2*P->c[P->Myxz[i] + P->nM*P->Myxz[i+P->ny]];
        if(ix > 0) {
          P->b[i] +=  P->dt/2*P->c[P->Myxz[i] + P->nM*P->Myxz[i-P->ny]];
          float w  = -P->dt/2*P->c[P->Myxz[i] + P->nM*P->Myxz[i-P->ny]]/P->b[i - P->ny];
          if(ix > 1 || !P->heatSinked) P->b[i] += w*P->dt/2*P->c[P->Myxz[i-P->ny] + P->nM*P->Myxz[i]];
          P->Tyxz[i] -= w*P->Tyxz[i-P->ny];
        }
      }
    }

    for(long ix=P->nx-1; ix>=0; ix--) {
      if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) {
        long i = iy + ix*P->ny + iz*P->ny*P->nx;
        P->Tyxz[i] = (P->Tyxz[i] + (ix == P->nx-1? 0: P->dt/2*P->c[P->Myxz[i] + P->nM*P->Myxz[i+P->ny]]*P->Tyxz[i+P->ny]))/P->b[i];
      }
    }
  }
  #else
  long ix,iy,iz;
  struct parameters *P = P_global;
  #ifdef _OPENMP
  long threadNum = omp_get_thread_num();
  #pragma omp for schedule(dynamic)
  #else
  long threadNum = 0;
  #endif
  for(iz=0; iz<P->nz; iz++) {
    for(iy=0; iy<P->ny; iy++) {
      for(ix=0; ix<P->nx; ix++) {
        long ib = ix + threadNum*P->nx;
        long i = ix + iy*P->nx + iz*P->nx*P->ny;
      // Thomson algorithm, sweeps up from 0 to nx-1 and then down from nx-1 to 0:
      // Algorithm is taken from https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
        P->b[ib] = 1;
        if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) {
          if(ix < P->nx-1) P->b[ib] += P->dt/2*P->c[P->M[i] + P->nM*P->M[i+1]];
          if(ix > 0) {
            P->b[ib] +=  P->dt/2*P->c[P->M[i] + P->nM*P->M[i-1]];
            float w   = -P->dt/2*P->c[P->M[i] + P->nM*P->M[i-1]]/P->b[ib - 1];
            if(ix > 1 || !P->heatSinked) P->b[ib] += w*P->dt/2*P->c[P->M[i-1] + P->nM*P->M[i]];
            P->T2[i] -= w*P->T2[i-1];
          }
        }
      }

      for(ix=P->nx-1; ix>=0; ix--) {
        if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) {
          long i = ix + iy*P->nx + iz*P->nx*P->ny;
          long ib = ix + threadNum*P->nx;
          P->T2[i] = (P->T2[i] + (ix == P->nx-1? 0: P->dt/2*P->c[P->M[i] + P->nM*P->M[i+1]]*P->T2[i+1]))/P->b[ib];
        }
      }
    }
  }
  #endif
}

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void substep2a(struct parameters *P_global) {
  // Explicit part of substep 2 out of 3
  #ifdef __NVCC__
  __shared__ struct parameters P_var;
  struct parameters *P = &P_var;
  if(!threadIdx.x) P_var = *P_global; // Only let one thread per block do the copying. 
  __syncthreads(); // All threads in the block wait for the copy to have finished
  __shared__ float tile[TILE_DIM][TILE_DIM+1]; // +1 is to avoid memory bank conflicts

  long xTiles = (P->nx + TILE_DIM - 1)/TILE_DIM;
  long yTiles = (P->ny + TILE_DIM - 1)/TILE_DIM;
  for(long tileNum=blockIdx.x; tileNum<xTiles*yTiles*P->nz; tileNum += gridDim.x) {
    long tilexoffset = TILE_DIM*(tileNum%xTiles);
    long tileyoffset = TILE_DIM*(tileNum/xTiles%yTiles);
    long iz          =           tileNum/xTiles/yTiles;
    long ix = tilexoffset + threadIdx.y;
    long iy = tileyoffset + threadIdx.x;

    __syncthreads(); // All threads in the block wait for any previous tile usage to have completed
    if(ix<P->nx && iy<P->ny) tile[threadIdx.x][threadIdx.y] = P->Tyxz[iz*P->ny*P->nx + ix*P->ny + iy]; // load yx data and store in shared memory in tile, which is xy
    __syncthreads(); // All threads in the block wait for the copy to have finished

    ix = tilexoffset + threadIdx.x;
    iy = tileyoffset + threadIdx.y;
    if(ix<P->nx && iy<P->ny) {
      long i = ix + iy*P->nx + iz*P->nx*P->ny;
      float dTdt = 0;
      if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) { // If not constant-temperature boundary voxel
        if(iy != 0      ) dTdt -= (P->T1[i-P->nx] - P->T1[i])*P->c[P->M[i] + P->M[i-P->nx]*P->nM + P->nM*P->nM]/2;
        if(iy != P->ny-1) dTdt -= (P->T1[i+P->nx] - P->T1[i])*P->c[P->M[i] + P->M[i+P->nx]*P->nM + P->nM*P->nM]/2;
      }
      P->T2[i] = tile[threadIdx.y][threadIdx.x] + P->dt*dTdt;
    }
  }
  #else
  long ix,iy,iz;
  struct parameters *P = P_global;
  #ifdef _OPENMP
  #pragma omp for schedule(dynamic)
  #endif
  for(iz=0; iz<P->nz; iz++) {
    for(iy=0; iy<P->ny; iy++) {
      for(ix=0; ix<P->nx; ix++) {
        long i = ix + iy*P->nx + iz*P->nx*P->ny;
        if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) { // If not constant-temperature boundary voxel
          float dTdt = 0;
          if(iy != 0      ) dTdt -= (P->T1[i-P->nx] - P->T1[i])*P->c[P->M[i] + P->M[i-P->nx]*P->nM + P->nM*P->nM]/2;
          if(iy != P->ny-1) dTdt -= (P->T1[i+P->nx] - P->T1[i])*P->c[P->M[i] + P->M[i+P->nx]*P->nM + P->nM*P->nM]/2;
          P->T2[i] = P->T2[i] + P->dt*dTdt;
        }
      }
    }
  }
  #endif
}

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void substep2b(struct parameters *P_global) {
  // Implicit part of substep 2 out of 3
  #ifdef __NVCC__
  long threadNum = threadIdx.x + threadIdx.y*blockDim.x + blockIdx.x*blockDim.x*blockDim.y;
  __shared__ struct parameters P_var;
  struct parameters *P = &P_var;
  if(!threadIdx.x) P_var = *P_global; // Only let one thread per block do the copying. 
  __syncthreads(); // All threads in the block wait for the copy to have finished
  for(long ixz=threadNum;ixz<P->nx*P->nz;ixz+=gridDim.x*TILE_DIM*TILE_DIM) {
    long ix = ixz%P->nx;
    long iz = ixz/P->nx;
    for(long iy=0; iy<P->ny; iy++) {
      long i = ix + iy*P->nx + iz*P->nx*P->ny;
      P->b[i] = 1;
      if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) {
        if(iy < P->ny-1) P->b[i] += P->dt/2*P->c[P->M[i] + P->nM*P->M[i+P->nx] + P->nM*P->nM];
        if(iy > 0) {
          P->b[i] +=  P->dt/2*P->c[P->M[i] + P->nM*P->M[i-P->nx] + P->nM*P->nM];
          float w  = -P->dt/2*P->c[P->M[i] + P->nM*P->M[i-P->nx] + P->nM*P->nM]/P->b[i - P->nx];
          if(iy > 1 || !P->heatSinked) P->b[i] += w*P->dt/2*P->c[P->M[i-P->nx] + P->nM*P->M[i] + P->nM*P->nM];
          P->T2[i] -= w*P->T2[i-P->nx];
        }
      }
    }

    for(long iy=P->ny-1; iy>=0; iy--) {
      if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) {
        long i = ix + iy*P->nx + iz*P->nx*P->ny;
        P->T2[i] = (P->T2[i] + (iy == P->ny-1? 0: P->dt/2*P->c[P->M[i] + P->nM*P->M[i+P->nx] + P->nM*P->nM]*P->T2[i+P->nx]))/P->b[i];
      }
    }
  }
  #else
  long ix,iy,iz;
  struct parameters *P = P_global;
  #ifdef _OPENMP
  long threadNum = omp_get_thread_num();
  #pragma omp for schedule(dynamic)
  #else
  long threadNum = 0;
  #endif
  for(iz=0; iz<P->nz; iz++) {
    for(ix=0; ix<P->nx; ix++) {
      for(iy=0; iy<P->ny; iy++) {
        long ib = iy + threadNum*P->ny;
        P->b[ib] = 1;
        if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) {
          long i = ix + iy*P->nx + iz*P->nx*P->ny;
          if(iy < P->ny-1) P->b[ib] += P->dt/2*P->c[P->M[i] + P->nM*P->M[i+P->nx] + P->nM*P->nM];
          if(iy > 0) {
            P->b[ib] +=  P->dt/2*P->c[P->M[i] + P->nM*P->M[i-P->nx] + P->nM*P->nM];
            float w = -P->dt/2*P->c[P->M[i] + P->nM*P->M[i-P->nx] + P->nM*P->nM]/P->b[ib - 1];
            if(iy > 1 || !P->heatSinked) P->b[ib] += w*P->dt/2*P->c[P->M[i-P->nx] + P->nM*P->M[i] + P->nM*P->nM];
            P->T2[i] -= w*P->T2[i-P->nx];
          }
        }
      }

      for(iy=P->ny-1; iy>=0; iy--) {
        if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) {
          long ib = iy + threadNum*P->ny;
          long i = ix + iy*P->nx + iz*P->nx*P->ny;
          P->T2[i] = (P->T2[i] + (iy == P->ny-1? 0: P->dt/2*P->c[P->M[i] + P->nM*P->M[i+P->nx] + P->nM*P->nM]*P->T2[i+P->nx]))/P->b[ib];
        }
      }
    }
  }
  #endif
}

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void substep3a(struct parameters *P_global) {
  // Explicit part of substep 3 out of 3
  #ifdef __NVCC__
  long threadNum = threadIdx.x + threadIdx.y*blockDim.x + blockIdx.x*blockDim.x*blockDim.y;
  __shared__ struct parameters P_var;
  struct parameters *P = &P_var;
  if(!threadIdx.x) P_var = *P_global; // Only let one thread per block do the copying. 
  __syncthreads(); // All threads in the block wait for the copy to have finished
  for(long i=threadNum; i<P->nx*P->ny*P->nz; i+=gridDim.x*TILE_DIM*TILE_DIM) {
    long ix = i%P->nx;
    long iy = i/P->nx%P->ny;
    long iz = i/P->nx/P->ny;
    if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) { // If not constant-temperature boundary voxel
      float dTdt = 0;
      if(iz != 0      ) dTdt -= (P->T1[i-P->nx*P->ny] - P->T1[i])*P->c[P->M[i] + P->M[i-P->nx*P->ny]*P->nM + 2*P->nM*P->nM]/2;
      if(iz != P->nz-1) dTdt -= (P->T1[i+P->nx*P->ny] - P->T1[i])*P->c[P->M[i] + P->M[i+P->nx*P->ny]*P->nM + 2*P->nM*P->nM]/2;
      P->T2[i] = P->T2[i] + P->dt*dTdt;
    }
  }
  #else
  long ix,iy,iz;
  struct parameters *P = P_global;
  #ifdef _OPENMP
  #pragma omp for schedule(dynamic)
  #endif
  for(iz=0; iz<P->nz; iz++) {
    for(iy=0; iy<P->ny; iy++) {
      for(ix=0; ix<P->nx; ix++) {
        long i = ix + iy*P->nx + iz*P->nx*P->ny;
        if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) { // If not constant-temperature boundary voxel
          float dTdt = 0;
          if(iz != 0      ) dTdt -= (P->T1[i-P->nx*P->ny] - P->T1[i])*P->c[P->M[i] + P->M[i-P->nx*P->ny]*P->nM + 2*P->nM*P->nM]/2;
          if(iz != P->nz-1) dTdt -= (P->T1[i+P->nx*P->ny] - P->T1[i])*P->c[P->M[i] + P->M[i+P->nx*P->ny]*P->nM + 2*P->nM*P->nM]/2;
          P->T2[i] = P->T2[i] + P->dt*dTdt;
        }
      }
    }
  }
  #endif
}

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void substep3b(struct parameters *P_global) {
  // Implicit part of substep 3 out of 3
  #ifdef __NVCC__
  long threadNum = threadIdx.x + threadIdx.y*blockDim.x + blockIdx.x*blockDim.x*blockDim.y;
  __shared__ struct parameters P_var;
  struct parameters *P = &P_var;
  if(!threadIdx.x) P_var = *P_global; // Only let one thread per block do the copying. 
  __syncthreads(); // All threads in the block wait for the copy to have finished
  for(long ixy=threadNum;ixy<P->nx*P->ny;ixy+=gridDim.x*TILE_DIM*TILE_DIM) {
    long ix = ixy%P->nx;
    long iy = ixy/P->nx;
    for(long iz=0; iz<P->nz; iz++) {
      long i = ix + iy*P->nx + iz*P->nx*P->ny;
      P->b[i] = 1;
      if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) {
        if(iz < P->nz-1) P->b[i] += P->dt/2*P->c[P->M[i] + P->nM*P->M[i+P->nx*P->ny] + 2*P->nM*P->nM];
        if(iz > 0) {
          P->b[i]  +=  P->dt/2*P->c[P->M[i] + P->nM*P->M[i-P->nx*P->ny] + 2*P->nM*P->nM];
          float w   = -P->dt/2*P->c[P->M[i] + P->nM*P->M[i-P->nx*P->ny] + 2*P->nM*P->nM]/P->b[i - P->nx*P->ny];
          if(iz > 1 || !P->heatSinked) P->b[i] += w*P->dt/2*P->c[P->M[i-P->nx*P->ny] + P->nM*P->M[i] + 2*P->nM*P->nM];
          P->T2[i] -= w*P->T2[i-P->nx*P->ny];
        }
      }
    }
    for(long iz=P->nz-1; iz>=0; iz--) {
      long i = ix + iy*P->nx + iz*P->nx*P->ny;
      if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) {
        P->T2[i] = (P->T2[i] + (iz == P->nz-1? 0: P->dt/2*P->c[P->M[i] + P->nM*P->M[i+P->nx*P->ny] + 2*P->nM*P->nM]*P->T2[i+P->nx*P->ny]))/P->b[i];
      }
      if(P->calcDamage && P->E[P->M[i]]) { // Arrhenius damage integral evaluation
        P->Omega[i] += P->dt*(expf(P->logA[P->M[i]]-P->E[P->M[i]]/(R*((P->T2[i] + P->T1[i])/2 + CELSIUSZERO))));
      }
      atomicMaxWrapperFLT(&P->maxMediaTemps[P->M[i]], P->T2[i]); // If the temperature in this voxel is a new record high for this medium, write the value into the maxMediaTemps array
    }
  }
  #else
  long ix,iy,iz;
  struct parameters *P = P_global;
  #ifdef _OPENMP
  long threadNum = omp_get_thread_num();
  #pragma omp for schedule(dynamic)
  #else
  long threadNum = 0;
  #endif
  for(iy=0; iy<P->ny; iy++) {
    for(iz=0; iz<P->nz; iz++) {
      for(ix=0; ix<P->nx; ix++) { // The x sweep is put inside the z sweep on purpose, since it improves performance a lot in this step. The b coefficients used here for the Thomson algorithm are stored in a 2D (x,z) matrix.
        long ib = ix + iz*P->nx + threadNum*P->nx*P->nz;
        P->b[ib] = 1;
        if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) {
          long i = ix + iy*P->nx + iz*P->nx*P->ny;
          if(iz < P->nz-1) P->b[ib] += P->dt/2*P->c[P->M[i] + P->nM*P->M[i+P->nx*P->ny] + 2*P->nM*P->nM];
          if(iz > 0) {
            P->b[ib]  +=  P->dt/2*P->c[P->M[i] + P->nM*P->M[i-P->nx*P->ny] + 2*P->nM*P->nM];
            float w    = -P->dt/2*P->c[P->M[i] + P->nM*P->M[i-P->nx*P->ny] + 2*P->nM*P->nM]/P->b[ib - P->nx];
            if(iz > 1 || !P->heatSinked) P->b[ib] += w*P->dt/2*P->c[P->M[i-P->nx*P->ny] + P->nM*P->M[i] + 2*P->nM*P->nM];
            P->T2[i] -= w*P->T2[i-P->nx*P->ny];
          }
        }
      }
    }
    for(iz=P->nz-1; iz>=0; iz--) {
      for(ix=0; ix<P->nx; ix++) {
        long ib = ix + iz*P->nx + threadNum*P->nx*P->nz;
        long i = ix + iy*P->nx + iz*P->nx*P->ny;
        if(!(P->heatSinked && (ix == 0 || ix == P->nx-1 || iy == 0 || iy == P->ny-1 || iz == 0 || iz == P->nz-1))) {
          P->T2[i] = (P->T2[i] + (iz == P->nz-1? 0: P->dt/2*P->c[P->M[i] + P->nM*P->M[i+P->nx*P->ny] + 2*P->nM*P->nM]*P->T2[i+P->nx*P->ny]))/P->b[ib];
        }
        if(P->calcDamage && P->E[P->M[i]]) { // Arrhenius damage integral evaluation
          P->Omega[i] += P->dt*(expf(P->logA[P->M[i]]-P->E[P->M[i]]/(R*((P->T2[i] + P->T1[i])/2 + CELSIUSZERO))));
        }
        atomicMaxWrapperFLT(&P->maxMediaTemps[P->M[i]], P->T2[i]); // If the temperature in this voxel is a new record high for this medium, write the value into the maxMediaTemps array
      }
    }
  }
  #endif
}

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void recordTempSensors(struct parameters *P,long n) {
  for(long j=0; j<P->nSensors; j++) { // Interpolate to get the new temperatures on the temperature sensors
    long idx = (long)P->tempSensorCornerIdxs[j];
    float wx = (float)P->tempSensorInterpWeights[j              ];
    float wy = (float)P->tempSensorInterpWeights[j+  P->nSensors];
    float wz = (float)P->tempSensorInterpWeights[j+2*P->nSensors];
    P->sensorTemps[j+(n+1)*P->nSensors] = (1-wx)*(1-wy)*(1-wz)*P->T1[idx        ] + (1-wx)*(1-wy)*wz*P->T1[idx        +P->nx*P->ny] +
                                          (1-wx)*   wy *(1-wz)*P->T1[idx  +P->nx] + (1-wx)*   wy *wz*P->T1[idx  +P->nx+P->nx*P->ny] +
                                             wx *(1-wy)*(1-wz)*P->T1[idx+1      ] +    wx *(1-wy)*wz*P->T1[idx+1      +P->nx*P->ny] +
                                             wx *   wy *(1-wz)*P->T1[idx+1+P->nx] +    wx *   wy *wz*P->T1[idx+1+P->nx+P->nx*P->ny];
  }
}

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void swapTempPointers(struct parameters *P, long n) {
  #ifdef __NVCC__
  float *temp = P->T1;
  P->T1 = P->T2;
  P->T2 = temp;
  #else
  if(n>0) { // Swap T1 and T2
    float *temp = P->T1;
    P->T1 = P->T2;
    P->T2 = temp;
  } else if(P->nt%2) {
    P->T1 = P->T2;
    P->T2 = (float *)malloc(P->nx*P->ny*P->nz*sizeof(float));
  } else {
    P->T1 = P->T2;
    P->T2 = P->Tfinal;
  }
  #endif
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
  bool ctrlc_caught = false;           // Has a ctrl+c been passed from MATLAB?
  
  struct parameters P_var;
  struct parameters *P = &P_var;

  mwSize const *dimPtr = mxGetDimensions(prhs[0]);
  plhs[0] = mxCreateNumericArray(3,dimPtr,mxSINGLE_CLASS,mxREAL);
  P->nx = (long)dimPtr[0];
  P->ny = (long)dimPtr[1];
  P->nz = (long)dimPtr[2];
  P->nM = (long)mxGetM(mxGetField(prhs[2],0,"dTdtperdeltaT")); // Number of media in the simulation.
  
  P->nt               = (long)*mxGetPr(mxGetField(prhs[2],0,"steps")); // nt is an integer (number of time steps to perform)
  P->T1               = (float *)mxGetData(prhs[0]); // T1 is an nx*ny*nz array of floats (singles in MATLAB)
  P->Tfinal           = (float *)mxGetData(plhs[0]);
  #ifndef __NVCC__
  P->T2               = (float *)(P->nt%2? P->Tfinal: malloc(P->nx*P->ny*P->nz*sizeof(float))); // T2 is an nx*ny*nz array of floats
  #endif
  P->dt               = (float)*mxGetPr(mxGetField(prhs[2],0,"dt")); // dt in MATLAB is a double that is here converted to a float (duration of time step)
  P->M                = (unsigned char *)mxGetData(mxGetField(prhs[2],0,"M")); // M is a nx*ny*nz array of uint8 (unsigned char) containing values from 0..nM-1
  P->c                = (float *)mxGetData(mxGetField(prhs[2],0,"dTdtperdeltaT")); // c is an nM*nM*3 array of floats
  P->dTdt_abs         = (float *)mxGetData(mxGetField(prhs[2],0,"dTdt_abs")); // dTdt_abs is an nx*ny*nz array of floats
  double *Amatlab     = mxGetPr(mxGetField(prhs[2],0,"A")); // A is a nM array of doubles
  P->logA             = (float *)malloc(P->nM*sizeof(float));
  double *Ematlab     = mxGetPr(mxGetField(prhs[2],0,"E")); // E is a nM array of doubles
  P->E                = (float *)malloc(P->nM*sizeof(float));
  for(long j=0;j<P->nM;j++) {
    P->logA[j] = (float)log(Amatlab[j]);
    P->E[j] = (float)Ematlab[j];
  }
  P->lightsOn         = mxIsLogicalScalarTrue(mxGetField(prhs[2],0,"lightsOn"));
  P->heatSinked       = (long)*mxGetPr(mxGetField(prhs[2],0,"heatBoundaryType")); // 0: Insulating boundaries, 1: Constant-temperature boundaries
  
  float *Omega_in     = (float *)mxGetData(prhs[1]); // Omega is an nx*ny*nz array of floats if we are supposed to calculate damage, a single NaN element otherwise
  P->calcDamage     = !mxIsNaN(Omega_in[0]); // If the Omega input is just one NaN element then we shouldn't bother with thermal damage calculation
  plhs[1] = P->calcDamage? mxCreateNumericArray(3,dimPtr,mxSINGLE_CLASS,mxREAL): mxCreateNumericMatrix(1,1,mxSINGLE_CLASS,mxREAL); // Omega is the same dimensions as Omega_in
  P->Omega = (float *)mxGetData(plhs[1]);
  if(!P->calcDamage) P->Omega[0] = NAN;
  
  P->tempSensorCornerIdxs = mxGetPr(mxGetField(prhs[2],0,"tempSensorCornerIdxs"));
  P->nSensors = (long)mxGetM(mxGetField(prhs[2],0,"tempSensorCornerIdxs"));
  P->tempSensorInterpWeights = mxGetPr(mxGetField(prhs[2],0,"tempSensorInterpWeights"));
  
  plhs[2] = P->nSensors? mxCreateDoubleMatrix(P->nSensors,P->nt+1,mxREAL): mxCreateDoubleMatrix(0,0,mxREAL);
  P->sensorTemps = mxGetPr(plhs[2]);
  
  plhs[3] = mxCreateNumericMatrix(1,P->nM,mxSINGLE_CLASS,mxREAL);
  P->maxMediaTemps = (float *)mxGetData(plhs[3]);
  for(long j=0;j<P->nM;j++) P->maxMediaTemps[j] = -INFINITY;
  
  if(P->calcDamage) for(long i=0; i<P->nx*P->ny*P->nz; i++) P->Omega[i] = Omega_in[i]; // Initialize omega array
  
  #ifdef __NVCC__
  int temp, nBlocks; gpuErrchk(cudaOccupancyMaxPotentialBlockSize(&nBlocks,&temp,&swapTempPointers,0,0));
  dim3 blockDims(TILE_DIM,TILE_DIM,1);

  P->b = (float *)malloc((P->nx*P->ny*P->nz)*sizeof(float));
  struct parameters *P_dev;
  struct debug D_var = {{0.0,0.0,0.0},{0,0,0}};
  struct debug *D = &D_var;
  struct debug *D_dev;
  createDeviceStructs(P,&P_dev,D,&D_dev);
  recordTempSensors<<<1,1>>>(P_dev,-1);
  #else
  recordTempSensors(P,-1);
  #ifdef _OPENMP
  bool useAllCPUs       = mxIsLogicalScalarTrue(mxGetField(prhs[2],0,"useAllCPUs"));
  long numThreads = useAllCPUs || omp_get_num_procs() == 1? omp_get_num_procs(): omp_get_num_procs()-1;
  P->b = malloc(numThreads*(P->nx*P->nz > P->ny? P->nx*P->nz: P->ny)*sizeof(float));
  #pragma omp parallel num_threads(numThreads)
  #else
  P->b = malloc((P->nx*P->nz > P->ny? P->nx*P->nz: P->ny)*sizeof(float));
  #endif // _OPENMP
  #endif // __NVCC__
  {
    for(long n=0; n<P->nt; n++) {
      if(ctrlc_caught) break;
      
      // In the following, the Douglas-Gunn Alternating Direction Implicit (DG-ADI) method is used to propagate the heat. This is based on equations (3.23a) and (3.23b) in me690-lctr-nts.pdf.
      #ifdef __NVCC__
      substep1a<<<nBlocks, blockDims>>>(P_dev); // xyz -> yxz
      substep1b<<<nBlocks, blockDims>>>(P_dev); // yxz -> yxz
      substep2a<<<nBlocks, blockDims>>>(P_dev); // yxz -> xyz
      substep2b<<<nBlocks, blockDims>>>(P_dev); // xyz -> xyz
      substep3a<<<nBlocks, blockDims>>>(P_dev); // xyz -> xyz
      substep3b<<<nBlocks, blockDims>>>(P_dev); // xyz -> xyz
      gpuErrchk(cudaDeviceSynchronize()); // Wait until all kernels have finished
      #else
      substep1a(P);
      substep1b(P);
      substep2a(P);
      substep2b(P);
      substep3a(P);
      substep3b(P);
      #endif

      #ifdef _OPENMP
      #pragma omp master
      #endif
      {
        #ifndef __clang__
        if(utIsInterruptPending()) {
          ctrlc_caught = true;
          printf("\nCtrl+C detected, stopping.\n");
        }
        #endif

        #ifdef __NVCC__
        swapTempPointers<<<1,1>>>(P_dev,n);
        recordTempSensors<<<1,1>>>(P_dev,n);
        #else
        swapTempPointers(P,n);
        recordTempSensors(P,n);
        #endif
      }
      #ifdef _OPENMP
      #pragma omp barrier
      #endif
    }
  }
  #ifdef __NVCC__
  gpuErrchk(cudaDeviceSynchronize()); // Wait until all kernels have finished
  retrieveAndFreeDeviceStructs(P,P_dev,D,D_dev);
//   printf("\nDebug: %.18e %.18e %.18e %llu %llu %llu\n          ",D->dbls[0],D->dbls[1],D->dbls[2],D->ulls[0],D->ulls[1],D->ulls[2]);
  #else
  if(P->nt > 1) free(P->T2);
  #endif

  free(P->b);
  free(P->logA);
  free(P->E);
  return;
}
