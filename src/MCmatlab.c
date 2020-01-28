/********************************************
 *
 * MCmatlab.c, in the C programming language, written for MATLAB MEX function generation
 * C script for Monte Carlo Simulation of Photon Transport in 3D
 * 
 * Copyright 2017, 2018 by Anders K. Hansen, DTU Fotonik
 * inspired by mcxyz.c by Steven Jacques, Ting Li, Scott Prahl at the Oregon Health & Science University
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
 * (Set the MATLAB current folder to the one with all the example files)
 * Can be compiled in MATLAB using the MinGW-w64 compiler (GCC) with
 * "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' -outdir helperfuncs\private .\src\MCmatlab.c ".\src\libut.lib""
 * ... or the Microsoft Visual C++ compiler (MSVC) with
 * "mex COMPFLAGS='/Zp8 /GR /EHs /nologo /MD /openmp /W4 /WX /wd4204 /wd4100' -outdir helperfuncs\private .\src\MCmatlab.c ".\src\libut.lib""
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
 * "copyfile ./src/MCmatlab.c ./src/MCmatlab_CUDA.cu; mexcuda -llibut COMPFLAGS='-use_fast_math -res-usage $COMPFLAGS' -outdir helperfuncs\private .\src\MCmatlab_CUDA.cu ".\src\nvml.lib""
 * 
 ** COMPILING ON MAC
 * As of June 2017, the macOS compiler doesn't support libut (for ctrl+c 
 * breaking) or openmp (for multithreading).
 * Compile in MATLAB with "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -std=c11 -Wall' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -std=c11 -Wall' -outdir helperfuncs/private ./src/MCmatlab.c"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Install XCode from the App Store
 * 2. Type "mex -setup" in the MATLAB command window
 ********************************************/
// printf("Debug 1...\n");mexEvalString("drawnow;");mexEvalString("drawnow;");mexEvalString("drawnow;"); // For inserting into code for debugging purposes

#include "mex.h"
#include <math.h>

#ifdef __GNUC__ // This is defined for GCC and CLANG but not for Microsoft Visual C++ compiler
  #include <time.h>
  #define min(a,b) ({__typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b? _b: _a;})
  #define max(a,b) ({__typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b? _a: _b;})
#else
  #include <float.h>
  #include <windows.h>
#endif
#include "lambert.c" // For calculating the Lambert W function, originally part of the GNU Scientific Library, created by K. Briggs, G. Jungman and B. Gough and slightly modified by A. Hansen for easier MCmatlab integration

#ifndef __clang__ // If compiling with CLANG (on Mac)
  #ifdef __cplusplus
  extern "C"
  #endif
  extern bool utIsInterruptPending(void); // Allows catching ctrl+c while executing the mex function
#endif

#define USEFLOATSFORGPU 1 // Comment out this line if you want the Monte Carlo routine to work internally with double-precision floating point numbers
#if defined(USEFLOATSFORGPU) && defined(__NVCC__)
  typedef float FLOATORDBL;
  #define FLOATORDBLEPS FLT_EPSILON
  #define SIN(x)     sinf(x)
  #define ASIN(x)    asinf(x)
  #define COS(x)     cosf(x)
  #define ACOS(x)    acosf(x)
  #define TAN(x)     tanf(x)
  #define ATAN(x)    atanf(x)
  #define ATAN2(x,y) atan2f(x,y)
  #define EXP(x)     expf(x)
  #define EXPM1(x)   expm1f(x)
  #define LOG(x)     logf(x)
  #define SQRT(x)    sqrtf(x)
  #define FLOOR(x)   floorf(x)
  #define FABS(x)    fabsf(x)
  #define SQR(x)     sqrf(x)
#else
  typedef double FLOATORDBL;
  #define FLOATORDBLEPS DBL_EPSILON
  #define SIN(x)     sin(x)
  #define ASIN(x)    asin(x)
  #define COS(x)     cos(x)
  #define ACOS(x)    acos(x)
  #define TAN(x)     tan(x)
  #define ATAN(x)    atan(x)
  #define ATAN2(x,y) atan2(x,y)
  #define EXP(x)     exp(x)
  #define EXPM1(x)   expm1(x)
  #define LOG(x)     log(x)
  #define SQRT(x)    sqrt(x)
  #define FLOOR(x)   floor(x)
  #define FABS(x)    fabs(x)
  #define SQR(x)     sqr(x)
#endif

#ifdef __NVCC__ // If compiling for CUDA
  #ifdef USEFLOATSFORGPU
    #define RandomNum   curand_uniform(&P->PRNGstate) // Calls for a random float in (0,1]
  #else
    #define RandomNum   curand_uniform_double(&P->PRNGstate) // Calls for a random double in (0,1]
  #endif
  #define DEVICE 0
  #define GPUHEAPMEMORYLIMIT 500000000 // Bytes to allow the GPU to allocate, mostly for tracking photons prior to storage in NFRdet. Must be less than the available memory.
  #define KERNELTIME 50000 // Microseconds to spend in each kernel call (watchdog timer is said to termine kernel function that take more than a few seconds)
  #include <nvml.h>
  #include <curand_kernel.h>
  typedef curandStateXORWOW_t PRNG_t;
//   typedef curandStateMRG32k3a_t PRNG_t;
//   typedef curandStatePhilox4_32_10_t PRNG_t;
  #define THREADNUM (threadIdx.x + blockDim.x*blockIdx.x)
#else
  #define DSFMT_MEXP 19937 // Mersenne exponent for dSFMT
  #include "dSFMT-src-2.2.3/dSFMT.c" // Double precision SIMD oriented Fast Mersenne Twister(dSFMT)
  #define RandomNum   dsfmt_genrand_open_close(&P->PRNGstate) // Calls for a random number in (0,1]
  typedef dsfmt_t PRNG_t;
  #ifdef _OPENMP
    #include <omp.h>
    #define THREADNUM omp_get_thread_num()
  #else
    #define THREADNUM 0
  #endif
#endif

#define PI          ACOS(-1.0f)
#define C           (FLOATORDBL)29979245800 // speed of light in vacuum in cm/s
#define THRESHOLD   (FLOATORDBL)0.01    // used in roulette
#define CHANCE      (FLOATORDBL)0.1      // used in roulette
#define SIGN(x)     ((x)>=0? 1:-1)
#define INITIALPATHSSIZE 2000
#define INITIALRECORDSIZE 1000
#define KILLRANGE   5 // Must be odd integer
    /* KILLRANGE determines the region that photons are allowed to stay 
     * alive in in multiples of the cuboid size (if outside, the probability
     * of returning to the region of interest is judged as too low). When
     * launching an infinite plane wave without boundaries, photons will be
     * launched in this whole extended region. */

#include "MCmatlablib.c"

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void threadInitAndLoop(struct beam *B_global, struct geometry *G_global,
          struct lightCollector *LC_global, struct paths *Pa, struct outputs *O_global, long nM,
          long long simulationTimeStart, long long microSecondsOrGPUCycles, unsigned long long nPhotonsRequested,
          bool *ctrlc_caughtPtr, bool silentMode, struct debug *D) {
  struct photon P_var;
  struct photon *P = &P_var;

  P->recordSize    = O_global->NFRdet? INITIALRECORDSIZE: 0; // If we're supposed to calculate NFRdet, we start the record at a size of 1000 elements - it will be dynamically expanded later if needed
  P->j_record      = O_global->NFRdet? (long *)malloc(P->recordSize*sizeof(long)): NULL;
  P->weight_record = O_global->NFRdet? (FLOATORDBL *)malloc(P->recordSize*sizeof(FLOATORDBL)): NULL;

  #ifdef __NVCC__ // If compiling for CUDA
  // Copy structs from global device memory to shared device memory, which is orders of magnitude faster since it is on-chip
  __shared__ struct beam B_var;
  struct beam *B = &B_var;
  __shared__ struct geometry G_var;
  struct geometry *G = &G_var;
  __shared__ struct lightCollector LC_var;
  struct lightCollector *LC = &LC_var;
  __shared__ struct outputs O_var;
  struct outputs *O = &O_var;
  extern __shared__ FLOATORDBL smallArrays[]; // Dynamically allocated shared memory with size implicitly specified by the third cuda kernel launch parameter
  if(!threadIdx.x) { // Only let one thread per block do the copying
    // Copy contents of structs and smallArrays (not deep copies, pointers still point to global device memory)
    B_var = *B_global;
    G_var = *G_global;
    LC_var = *LC_global;
    O_var = *O_global;
    memcpy(smallArrays,G->muav,(nM*3 + G->n[2] + B->L_NF1 + B->L_FF1 + B->L_NF2 + B->L_FF2)*sizeof(FLOATORDBL));
    // Set all array pointers to the correct new locations in the shared memory version of smallArrays
    G->muav = smallArrays;
    G->musv = G->muav + nM;
    G->gv = G->musv + nM;
    G->RIv = G->gv + nM;
    B->NFdist1 = G->RIv + G->n[2];
    B->FFdist1 = B->NFdist1 + B->L_NF1;
    B->NFdist2 = B->FFdist1 + B->L_FF1;
    B->FFdist2 = B->NFdist2 + B->L_NF2;
  }
  __syncthreads(); // All threads in the block wait for the copy to have finished
  
  // Initialize the PRNG and timing for the major loop
  curand_init(simulationTimeStart, THREADNUM, 0, &P->PRNGstate);
  long long tEnd = clock64() + microSecondsOrGPUCycles;

  // Launch major loop
  while(clock64() < tEnd && O->nPhotons + THREADNUM < nPhotonsRequested) { // "+ THREADNUM" ensures that we avoid race conditions that might launch more than nPhotonsRequested photons
  #else
  struct beam *B = B_global;
  struct geometry *G = G_global;
  struct lightCollector *LC= LC_global;
  struct outputs *O = O_global;
  // Check for failed memory allocations and initialize the PRNG
  if(P->recordSize && !P->j_record) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");
  if(P->recordSize && !P->weight_record) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");
  dsfmt_init_gen_rand(&P->PRNGstate,(unsigned long)simulationTimeStart + THREADNUM); // Seed the photon's random number generator
  int pctProgress = 0;      // Simulation progress in percent
  // Launch major loop
  while(pctProgress < 100 && O->nPhotons + THREADNUM < nPhotonsRequested && !*ctrlc_caughtPtr) { // "+ THREADNUM" ensures that we avoid race conditions that might launch more than nPhotonsRequested photons
  #endif
    atomicAddWrapperULL(&O_global->nPhotons,1); // We have to store the photon number in the global memory so it's visible to all blocks
    launchPhoton(P,B,G,Pa);

    while(P->alive) {
      while(P->stepLeft>0) {
        if(!P->sameVoxel) {
          checkEscape(P,G,LC,O);
          if(!P->alive) break;
          getNewVoxelProperties(P,G); // If photon has just entered a new voxel or has just been launched
        }
        propagatePhoton(P,G,O,Pa,D);
      }
      checkRoulette(P);
      scatterPhoton(P,G);
    }

    #ifndef __NVCC__
      // Check progress
      int pctTimeProgress = (int)(100.0*(getMicroSeconds() - simulationTimeStart)/microSecondsOrGPUCycles);
      int pctPhotonsProgress = (int)(100.0*O->nPhotons/nPhotonsRequested);

      #ifdef _OPENMP
      #pragma omp master
      #endif
      {
        // Check whether ctrl+c has been pressed
        #ifndef __clang__
        if(utIsInterruptPending()) {
          *ctrlc_caughtPtr = true;
          printf("\nCtrl+C detected, stopping.");
        }
        #endif

        // Print out message about progress.
        int newPctProgress = max(pctTimeProgress,pctPhotonsProgress);
        if(newPctProgress != pctProgress && !silentMode && !*ctrlc_caughtPtr) {
          printf("\b\b\b\b\b\b\b\b\b%3.i%% done", newPctProgress<100? newPctProgress: 100);
          mexEvalString("drawnow;");
        }
      }
      pctProgress = max(pctTimeProgress,pctPhotonsProgress);
    #endif
  }

  free(P->j_record); // Will do nothing if P->j_record == NULL
  free(P->weight_record); // Will do nothing if P->weight_record == NULL
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
  struct debug D_var = {{0.0,0.0,0.0},{0,0,0}};
  struct debug *D = &D_var;
  
  long idx;                  // General-purpose non-thread-specific index variable
  bool simFluorescence = *mxGetPr(prhs[1]) == 2;
  mxArray *MatlabMC = mxGetField(prhs[0],0,simFluorescence? "FMC": "MC");
  
  bool silentMode = mxIsLogicalScalarTrue(mxGetField(MatlabMC,0,"silentMode"));
  bool calcNFR    = mxIsLogicalScalarTrue(mxGetField(MatlabMC,0,"calcNFR")); // Are we supposed to calculate the NFR matrix?
  bool calcNFRdet = mxIsLogicalScalarTrue(mxGetField(MatlabMC,0,"calcNFRdet")); // Are we supposed to calculate the NFRdet matrix?
  
  // To know if we are simulating fluorescence, we check if a "sourceDistribution" field exists. If so, we will use it later in the beam definition.
  mxArray *MatlabBeam = mxGetField(MatlabMC,0,"beam");
  double *S_PDF       = (double *)mxGetData(mxGetField(MatlabBeam,0,"sourceDistribution")); // Power emitted by the individual voxels per unit volume. Can be percieved as an unnormalized probability density function of the 3D source distribution
  
  // Variables for timekeeping and number of photons
  bool            simulationTimed = mxIsNaN(*mxGetPr(mxGetField(MatlabMC,0,"nPhotonsRequested")));
  bool            ctrlc_caught = false;
  double          simulationTimeRequested = simulationTimed? *mxGetPr(mxGetField(MatlabMC,0,"simulationTimeRequested")): INFINITY;
  unsigned long long nPhotonsRequested = simulationTimed? ULLONG_MAX: (unsigned long long)*mxGetPr(mxGetField(MatlabMC,0,"nPhotonsRequested"));
  
  // Find out how much total memory to allocate for the small arrays (the array that for GPUs would be stored in GPU shared memory)
  mxArray *mediaProperties = mxGetField(MatlabMC,0,"mediaProperties");
  int     nM = (int)mxGetN(mediaProperties);
  mwSize const *dimPtr = mxGetDimensions(mxGetField(MatlabMC,0,"M"));
  long beamType = S_PDF? -1: (int)*mxGetPr(mxGetField(MatlabBeam,0,"beamType"));
  mxArray *MatlabBeamNF = S_PDF? 0: mxGetField(MatlabBeam,0,"NF");
  mxArray *MatlabBeamFF = S_PDF? 0: mxGetField(MatlabBeam,0,"FF");
  long L_NF1 = beamType >= 4? (long)mxGetN(mxGetField(MatlabBeamNF,0,beamType == 4? "radialDistr": "XDistr")): 0;
  long L_FF1 = beamType >= 4? (long)mxGetN(mxGetField(MatlabBeamFF,0,beamType == 4? "radialDistr": "XDistr")): 0;
  long L_NF2 = beamType == 5? (long)mxGetN(mxGetField(MatlabBeamNF,0,"YDistr")): 0;
  long L_FF2 = beamType == 5? (long)mxGetN(mxGetField(MatlabBeamFF,0,"YDistr")): 0;
  size_t size_smallArrays = (nM*3 + dimPtr[2] + L_NF1 + L_FF1 + L_NF2 + L_FF2)*sizeof(FLOATORDBL);
  FLOATORDBL *smallArrays = (FLOATORDBL *)malloc(size_smallArrays);

  // Geometry struct definition
  long L = (long)(dimPtr[0]*dimPtr[1]*dimPtr[2]); // Total number of voxels in cuboid
  mxArray *MatlabG         = mxGetField(prhs[0],0,"G");
  struct geometry G_var;
  struct geometry *G = &G_var;
  G->d[0] = (FLOATORDBL)*mxGetPr(mxGetField(MatlabG,0,"dx"));
  G->d[1] = (FLOATORDBL)*mxGetPr(mxGetField(MatlabG,0,"dy"));
  G->d[2] = (FLOATORDBL)*mxGetPr(mxGetField(MatlabG,0,"dz"));
  G->n[0] = (long)dimPtr[0];
  G->n[1] = (long)dimPtr[1];
  G->n[2] = (long)dimPtr[2];
  G->farFieldRes = (long)*mxGetPr(mxGetField(MatlabMC,0,"farFieldRes"));
  G->boundaryType = (int)*mxGetPr(mxGetField(MatlabMC,0,"boundaryType")); // boundaryType
  G->muav = smallArrays;
  G->musv = G->muav + nM;
  G->gv   = G->musv + nM;
  G->M = (unsigned char *)malloc(L*sizeof(unsigned char)); // M
  G->RIv  = G->gv + nM;

  // Fill the geometry arrays
  for(idx=0;idx<nM;idx++) {
    G->muav[idx] = (FLOATORDBL)*mxGetPr(mxGetField(mediaProperties,idx,"mua"));
    G->musv[idx] = (FLOATORDBL)*mxGetPr(mxGetField(mediaProperties,idx,"mus"));
    G->gv[idx]   = (FLOATORDBL)*mxGetPr(mxGetField(mediaProperties,idx,"g"));
  }
  unsigned char *M_matlab = (unsigned char *)mxGetData(mxGetField(MatlabMC,0,"M"));
  for(idx=0;idx<L;idx++) G->M[idx] = M_matlab[idx] - 1; // Convert from MATLAB 1-based indexing to C 0-based indexing
  double *RI_matlab = (double *)mxGetData(mxGetField(MatlabMC,0,"RI"));
  for(idx=0;idx<dimPtr[2];idx++) G->RIv[idx] = (FLOATORDBL)RI_matlab[idx];
  
  // Beam struct definition
  FLOATORDBL power        = 0;
  FLOATORDBL *S           = NULL; // Cumulative distribution function
  if(S_PDF) {
    S = (FLOATORDBL *)malloc((L+1)*sizeof(FLOATORDBL));
    if(!S) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");
    S[0] = 0;
    for(idx=1;idx<(L+1);idx++) S[idx] = S[idx-1] + (FLOATORDBL)S_PDF[idx-1];
    power = S[L]*G->d[0]*G->d[1]*G->d[2];
    for(idx=1;idx<(L+1);idx++) S[idx] /= S[L];
  }

  FLOATORDBL *NFdist1 = NULL, *NFdist2 = NULL, *FFdist1 = NULL, *FFdist2 = NULL;
  if(beamType >= 4) {
    double *MatlabNFdist1 = mxGetPr(mxGetField(MatlabBeamNF,0,beamType == 4? "radialDistr": "XDistr"));
    NFdist1 = G->RIv + dimPtr[2];
    if(L_NF1 == 1) {
      *NFdist1 = -1 - (FLOATORDBL)MatlabNFdist1[0];
    } else {
      NFdist1[0] = 0;
      for(idx=1;idx<L_NF1;idx++) NFdist1[idx] = NFdist1[idx-1] + (beamType == 4? idx-1: 1)*(FLOATORDBL)MatlabNFdist1[idx-1] + (beamType == 4? idx: 1)*(FLOATORDBL)MatlabNFdist1[idx];
      for(idx=1;idx<L_NF1;idx++) NFdist1[idx] /= NFdist1[L_NF1-1];
    }
    
    FFdist1 = NFdist1 + L_NF1;
    double *MatlabFFdist1 = mxGetPr(mxGetField(MatlabBeamFF,0,beamType == 4? "radialDistr": "XDistr"));
    if(L_FF1 == 1) {
      *FFdist1 = -1 - (FLOATORDBL)MatlabFFdist1[0];
    } else {
      FFdist1[0] = 0;
      for(idx=1;idx<L_FF1;idx++) FFdist1[idx] = FFdist1[idx-1] + (beamType == 4? idx-1: 1)*(FLOATORDBL)MatlabFFdist1[idx-1] + (beamType == 4? idx: 1)*(FLOATORDBL)MatlabFFdist1[idx];
      for(idx=1;idx<L_FF1;idx++) FFdist1[idx] /= FFdist1[L_FF1-1];
    }
    
    if(beamType == 5) {
      NFdist2 = FFdist1 + L_FF1;
      double *MatlabNFdist2 = mxGetPr(mxGetField(MatlabBeamNF,0,"YDistr"));
      if(L_NF2 == 1) {
        *NFdist2 = -1 - (FLOATORDBL)MatlabNFdist2[0];
      } else {
        NFdist2[0] = 0;
        for(idx=1;idx<L_NF2;idx++) NFdist2[idx] = NFdist2[idx-1] + (FLOATORDBL)MatlabNFdist2[idx-1] + (FLOATORDBL)MatlabNFdist2[idx];
        for(idx=1;idx<L_NF2;idx++) NFdist2[idx] /= NFdist2[L_NF2-1];
      }
    
      FFdist2 = NFdist2 + L_NF2;
      double *MatlabFFdist2 = mxGetPr(mxGetField(MatlabBeamFF,0,"YDistr"));
      if(L_FF2 == 1) {
        *FFdist2 = -1 - (FLOATORDBL)MatlabFFdist2[0];
      } else {
        FFdist2[0] = 0;
        for(idx=1;idx<L_FF2;idx++) FFdist2[idx] = FFdist2[idx-1] + (FLOATORDBL)MatlabFFdist2[idx-1] + (FLOATORDBL)MatlabFFdist2[idx];
        for(idx=1;idx<L_FF2;idx++) FFdist2[idx] /= FFdist2[L_FF2-1];
      }
    }
  }

  FLOATORDBL tb = (FLOATORDBL)(S? 0: *mxGetPr(mxGetField(MatlabBeam,0,"theta")));
  FLOATORDBL pb = (FLOATORDBL)(S? 0: *mxGetPr(mxGetField(MatlabBeam,0,"phi")));
  FLOATORDBL psib = (FLOATORDBL)(S? 0: *mxGetPr(mxGetField(MatlabBeam,0,"psi")));

  FLOATORDBL u[3] = {SIN(tb)*COS(pb),SIN(tb)*SIN(pb),COS(tb)}; // Temporary array
  
  FLOATORDBL w[3] = {1,0,0};
  FLOATORDBL v[3] = {0,0,1};
  if(u[2]!=1) unitcrossprod(u,v,w);
  axisrotate(w,u,psib,v);
  unitcrossprod(u,v,w);
  
  struct beam B_var = {
    S? 0: beamType,
    S? 0: NFdist1,
    S? 0: L_NF1,
    S? 0: (FLOATORDBL)(beamType == 5? *mxGetPr(mxGetField(MatlabBeamNF,0,"XWidth")): *mxGetPr(mxGetField(MatlabBeamNF,0,"radialWidth"))),
    S? 0: NFdist2,
    S? 0: L_NF2,
    S? 0: (FLOATORDBL)(beamType == 5? *mxGetPr(mxGetField(MatlabBeamNF,0,"YWidth")): 0),
    S? 0: FFdist1,
    S? 0: L_FF1,
    S? 0: (FLOATORDBL)(beamType == 5? *mxGetPr(mxGetField(MatlabBeamFF,0,"XWidth")): *mxGetPr(mxGetField(MatlabBeamFF,0,"radialWidth"))),
    S? 0: FFdist2,
    S? 0: L_FF2,
    S? 0: (FLOATORDBL)(beamType == 5? *mxGetPr(mxGetField(MatlabBeamFF,0,"YWidth")): 0),
    S,
    power,
    {(FLOATORDBL)(S? 0: *mxGetPr(mxGetField(MatlabBeam,0,"xFocus"))),
     (FLOATORDBL)(S? 0: *mxGetPr(mxGetField(MatlabBeam,0,"yFocus"))),
     (FLOATORDBL)(S? 0: *mxGetPr(mxGetField(MatlabBeam,0,"zFocus")))},
    {u[0],u[1],u[2]},
    {v[0],v[1],v[2]}, // normal vector to beam center axis
    {w[0],w[1],w[2]}
  };
  struct beam *B = &B_var;

  // Light Collector struct definition
  mxArray *MatlabLC      = mxGetField(MatlabMC,0,"LC");
  bool useLightCollector = mxIsLogicalScalarTrue(mxGetField(MatlabMC,0,"useLightCollector"));
  
  FLOATORDBL theta     = (FLOATORDBL)*mxGetPr(mxGetField(MatlabLC,0,"theta"));
  FLOATORDBL phi       = (FLOATORDBL)*mxGetPr(mxGetField(MatlabLC,0,"phi"));
  FLOATORDBL f         = (FLOATORDBL)*mxGetPr(mxGetField(MatlabLC,0,"f"));
  long   nTimeBins = S? 0: (long)*mxGetPr(mxGetField(MatlabLC,0,"nTimeBins"));

  struct lightCollector LC_var = {
    {(FLOATORDBL)*mxGetPr(mxGetField(MatlabLC,0,"x")) - (isfinite(f)? f*SIN(theta)*COS(phi):0), // r field, xyz coordinates of center of light collector
     (FLOATORDBL)*mxGetPr(mxGetField(MatlabLC,0,"y")) - (isfinite(f)? f*SIN(theta)*SIN(phi):0),
     (FLOATORDBL)*mxGetPr(mxGetField(MatlabLC,0,"z")) - (isfinite(f)? f*COS(theta)         :0)},
    theta,
    phi,
    f,
    (FLOATORDBL)*mxGetPr(mxGetField(MatlabLC,0,"diam")),
    (FLOATORDBL)*mxGetPr(mxGetField(MatlabLC,0,isfinite(f)?"fieldSize":"NA")),
    {(long)*mxGetPr(mxGetField(MatlabLC,0,"res")), nTimeBins? nTimeBins+2: 1},
    (FLOATORDBL)(S? 0: *mxGetPr(mxGetField(MatlabLC,0,"tStart"))),
    (FLOATORDBL)(S? 0: *mxGetPr(mxGetField(MatlabLC,0,"tEnd")))
  };
  struct lightCollector *LC = &LC_var;

  // Paths definitions (example photon trajectories)
  struct paths Pa_var;
  struct paths *Pa = &Pa_var;
  Pa->nExamplePaths = (long)*mxGetPr(mxGetField(MatlabMC,0,"nExamplePaths")); // How many photon path examples are we supposed to store?
  Pa->nMasterPhotonsLaunched = 0;
  Pa->pathsSize = INITIALPATHSSIZE*Pa->nExamplePaths; // Will be dynamically increased later if needed
  Pa->pathsElems = 0;
  Pa->data = Pa->nExamplePaths? (FLOATORDBL *)malloc(4*Pa->pathsSize*sizeof(FLOATORDBL)): NULL;
  if(Pa->pathsSize && !Pa->data) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");

  // Prepare output MATLAB struct
  plhs[0] = mxDuplicateArray(prhs[0]);
  
  mxArray *MCout = mxGetField(plhs[0],0,S? "FMC": "MC");
  mxArray *LCout = mxGetField(MCout,0,"LC");
  mxDestroyArray(mxGetField(MCout,0,"nPhotons"));
  mxDestroyArray(mxGetField(MCout,0,"nThreads"));
  mxDestroyArray(mxGetField(MCout,0,"simulationTime"));
  
  if(useLightCollector) mxDestroyArray(mxGetField(LCout,0,"image"));
  if(calcNFR) mxDestroyArray(mxGetField(MCout,0,"NFR"));
  if(calcNFRdet) mxDestroyArray(mxGetField(MCout,0,"NFRdet"));
  if(G->farFieldRes) mxDestroyArray(mxGetField(MCout,0,"farField"));
  if(G->boundaryType == 1) {
    mxDestroyArray(mxGetField(MCout,0,"NI_xpos"));
    mxDestroyArray(mxGetField(MCout,0,"NI_xneg"));
    mxDestroyArray(mxGetField(MCout,0,"NI_ypos"));
    mxDestroyArray(mxGetField(MCout,0,"NI_yneg"));
    mxDestroyArray(mxGetField(MCout,0,"NI_zpos"));
    mxDestroyArray(mxGetField(MCout,0,"NI_zneg"));
  } else if(G->boundaryType == 2) mxDestroyArray(mxGetField(MCout,0,"NI_zneg"));
  
  
  if(calcNFR)           mxSetField(MCout,0,"NFR",mxCreateNumericArray(3,dimPtr,mxDOUBLE_CLASS,mxREAL));
  if(calcNFRdet)        mxSetField(MCout,0,"NFRdet",mxCreateNumericArray(3,dimPtr,mxDOUBLE_CLASS,mxREAL));
  if(useLightCollector) {
    mwSize LCsize[3] = {(mwSize)LC->res[0],(mwSize)LC->res[0],(mwSize)LC->res[1]};
    mxSetField(LCout,0,"image", mxCreateNumericArray(3,LCsize,mxDOUBLE_CLASS,mxREAL));
  }
  if(G->farFieldRes)    mxSetField(MCout,0,"farField", mxCreateDoubleMatrix(G->farFieldRes,G->farFieldRes,mxREAL));
  if(G->boundaryType == 1) {
    mxSetField(MCout,0,"NI_xpos", mxCreateDoubleMatrix(G->n[1],G->n[2],mxREAL));
    mxSetField(MCout,0,"NI_xneg", mxCreateDoubleMatrix(G->n[1],G->n[2],mxREAL));
    mxSetField(MCout,0,"NI_ypos", mxCreateDoubleMatrix(G->n[0],G->n[2],mxREAL));
    mxSetField(MCout,0,"NI_yneg", mxCreateDoubleMatrix(G->n[0],G->n[2],mxREAL));
    mxSetField(MCout,0,"NI_zpos", mxCreateDoubleMatrix(G->n[0],G->n[1],mxREAL));
    mxSetField(MCout,0,"NI_zneg", mxCreateDoubleMatrix(G->n[0],G->n[1],mxREAL));
  } else if(G->boundaryType == 2) mxSetField(MCout,0,"NI_zneg", mxCreateDoubleMatrix(KILLRANGE*G->n[0],KILLRANGE*G->n[1],mxREAL));
  mxSetField(MCout,0,"simulationTime",mxCreateDoubleMatrix(1,1,mxREAL));
  mxSetField(MCout,0,"nPhotons",mxCreateDoubleMatrix(1,1,mxREAL));
  mxSetField(MCout,0,"nThreads",mxCreateDoubleMatrix(1,1,mxREAL));
  double *simTimePtr = mxGetPr(mxGetField(MCout,0,"simulationTime"));
  double *nPhotonsPtr = mxGetPr(mxGetField(MCout,0,"nPhotons"));
  double *nThreadsPtr = mxGetPr(mxGetField(MCout,0,"nThreads"));
  
  struct outputs O_var = {
    0, // nPhotons
  	calcNFR? mxGetPr(mxGetField(MCout,0,"NFR")): NULL,
  	calcNFRdet? mxGetPr(mxGetField(MCout,0,"NFRdet")): NULL,
    useLightCollector? mxGetPr(mxGetField(LCout,0,"image")): NULL,
    G->farFieldRes? mxGetPr(mxGetField(MCout,0,"farField")): NULL,
    G->boundaryType == 1? mxGetPr(mxGetField(MCout,0,"NI_xpos")): NULL,
    G->boundaryType == 1? mxGetPr(mxGetField(MCout,0,"NI_xneg")): NULL,
    G->boundaryType == 1? mxGetPr(mxGetField(MCout,0,"NI_ypos")): NULL,
    G->boundaryType == 1? mxGetPr(mxGetField(MCout,0,"NI_yneg")): NULL,
    G->boundaryType == 1? mxGetPr(mxGetField(MCout,0,"NI_zpos")): NULL,
    G->boundaryType != 0? mxGetPr(mxGetField(MCout,0,"NI_zneg")): NULL
  };
  struct outputs *O = &O_var;

  if(!silentMode) {
    // Display progress indicator
    printf(simFluorescence?"-----------Fluorescence Monte Carlo Simulation-----------\n":
                           "-----------------Monte Carlo Simulation------------------\n");
    if(simulationTimed) printf("Simulation duration = %0.3f min\n",simulationTimeRequested);
    else                printf("Requested # of photons = %0.2e\n",(double)nPhotonsRequested);
  }


  // ============================ MAJOR CYCLE ========================
  #ifdef __NVCC__ // If compiling for CUDA
  int threadsPerBlock, blocks; gpuErrchk(cudaOccupancyMaxPotentialBlockSize(&blocks,&threadsPerBlock,&threadInitAndLoop,size_smallArrays,0));
  *nThreadsPtr = threadsPerBlock*blocks;

  struct cudaDeviceProp CDP; gpuErrchk(cudaGetDeviceProperties(&CDP,DEVICE));
//   if(!silentMode) printf("Using %s with CUDA compute capability %d.%d,\n    launching %d blocks, each with %d threads,\n    using %d/%d bytes of shared memory per block\n",CDP.name,CDP.major,CDP.minor,blocks,threadsPerBlock,384+size_smallArrays,CDP.sharedMemPerBlock);
  if(!silentMode) printf("Using %s with CUDA compute capability %d.%d\n",CDP.name,CDP.major,CDP.minor);

  unsigned long long heapSizeLimit; gpuErrchk(cudaDeviceGetLimit(&heapSizeLimit,cudaLimitMallocHeapSize));
  if(calcNFRdet && heapSizeLimit < GPUHEAPMEMORYLIMIT) {
    gpuErrchk(cudaDeviceReset());
    gpuErrchk(cudaDeviceSetLimit(cudaLimitMallocHeapSize,GPUHEAPMEMORYLIMIT));
  }

  nvmlInit();
  nvmlDevice_t nvmldevice; int returnvar = nvmlDeviceGetHandleByIndex(DEVICE,&nvmldevice);
  unsigned int clock; nvmlDeviceGetClockInfo(nvmldevice,NVML_CLOCK_GRAPHICS,&clock);
//   unsigned int fanspeed; nvmlDeviceGetFanSpeed(nvmldevice,&fanspeed);
//   nvmlMemory_t memory; nvmlDeviceGetMemoryInfo(nvmldevice,&memory);
//   unsigned int T; nvmlDeviceGetTemperature(nvmldevice,NVML_TEMPERATURE_GPU,&T);
//   nvmlUtilization_t utilization; nvmlDeviceGetUtilizationRates(nvmldevice,&utilization);

  struct geometry *G_dev;
  struct beam *B_dev;
  struct lightCollector *LC_dev;
  struct paths *Pa_dev;
  struct outputs *O_dev;
  struct debug *D_dev;
  createDeviceStructs(G,&G_dev,B,&B_dev,LC,&LC_dev,Pa,&Pa_dev,O,&O_dev,nM,L,D,&D_dev);

  int pctProgress = 0;
  int callnumber = 0;
  long long timeLeft = simulationTimed? (long long)(simulationTimeRequested*60000000): LLONG_MAX; // microseconds
  if(!silentMode) {
//     printf("Calculating...   0%% done [GPU      MHz, memory used    /    GB]");
    printf("Calculating...   0%% done");
    mexEvalString("drawnow;");
  }
  long long simulationTimeStart = getMicroSeconds();
  long long prevtime = simulationTimeStart;
  do {
    // Run kernel
    threadInitAndLoop<<<blocks, threadsPerBlock, size_smallArrays>>>(B_dev,G_dev,LC_dev,Pa_dev,O_dev,nM,prevtime,clock*min((long long)KERNELTIME,timeLeft),nPhotonsRequested,NULL,false,D_dev);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());
    // Progress indicator
    long long newtime = getMicroSeconds();
    if(simulationTimed) {
      timeLeft = (long long)(simulationTimeRequested*60000000) - newtime + simulationTimeStart; // In microseconds
      pctProgress = (int)(100.0*(1.0 - timeLeft/(simulationTimeRequested*60000000.0)));
    } else {
      gpuErrchk(cudaMemcpy(O, O_dev, sizeof(unsigned long long),cudaMemcpyDeviceToHost)); // Copy just nPhotons, the first 8 bytes of O
      pctProgress = (int)(100.0*O->nPhotons/nPhotonsRequested);
    }
    nvmlDeviceGetClockInfo(nvmldevice,NVML_CLOCK_GRAPHICS,&clock);
    if(!(callnumber++%5) && !silentMode) {
//       nvmlDeviceGetFanSpeed(nvmldevice,&fanspeed);
//       nvmlDeviceGetMemoryInfo(nvmldevice,&memory);
//       nvmlDeviceGetTemperature(nvmldevice,NVML_TEMPERATURE_GPU,&T);
//       nvmlDeviceGetUtilizationRates(nvmldevice,&utilization);
//       printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%3d%% done [GPU %4u MHz, memory used %3.1f/%3.1f GB]",pctProgress,clock,(double)memory.used/1024/1024/1024,(double)memory.total/1024/1024/1024);
      printf("\b\b\b\b\b\b\b\b\b%3d%% done",pctProgress);
      mexEvalString("drawnow;");
    }
    prevtime = newtime;
    // Resize Pa->data if necessary
    if(Pa->nExamplePaths) {
      struct paths Pa_temp;
      gpuErrchk(cudaMemcpy(&Pa_temp, Pa_dev, sizeof(struct paths),cudaMemcpyDeviceToHost));
      if(Pa_temp.pathsElems == -1) {
        Pa->pathsSize *= 10; // Increase the record's size by a factor of 10
        Pa->data = (FLOATORDBL *)realloc(Pa->data,4*Pa->pathsSize*sizeof(FLOATORDBL));
        if(!Pa->data) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");
        gpuErrchk(cudaFree(Pa_temp.data));
        Pa_temp = *Pa;
        gpuErrchk(cudaMalloc(&Pa_temp.data, 4*Pa->pathsSize*sizeof(FLOATORDBL)));
        gpuErrchk(cudaMemcpy(Pa_dev,&Pa_temp,sizeof(struct paths),cudaMemcpyHostToDevice));
      }
    }
    // Check for ctrl+c
    if(utIsInterruptPending()) {
      ctrlc_caught = true;
      printf("\nCtrl+C detected, stopping.");
    }
  } while(timeLeft > 0 && O->nPhotons < nPhotonsRequested && !ctrlc_caught);
//   if(!silentMode) printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b100%% done [GPU %4u MHz, memory used %3.1f/%3.1f GB]",clock,(double)memory.used/1024/1024/1024,(double)memory.total/1024/1024/1024);
  if(!silentMode) printf("\b\b\b\b\b\b\b\b\b%3d%% done",pctProgress);

  nvmlShutdown();
  retrieveAndFreeDeviceStructs(G,G_dev,B,B_dev,LC,LC_dev,Pa,Pa_dev,O,O_dev,L,D,D_dev);

  #else

  if(!silentMode) {
    printf("Calculating...   0%% done");
    mexEvalString("drawnow;");
  }
  long long simulationTimeStart = getMicroSeconds();
  #ifdef _OPENMP
  bool useAllCPUs = mxIsLogicalScalarTrue(mxGetField(MatlabMC,0,"useAllCPUs"));
  *nThreadsPtr = useAllCPUs? omp_get_num_procs(): max(omp_get_num_procs()-1,1);
  #pragma omp parallel num_threads((long)*nThreadsPtr)
  #else
  *nThreadsPtr = 1;
  #endif
  {
    threadInitAndLoop(B,G,LC,Pa,O,nM,simulationTimeStart,(long long)(simulationTimeRequested*60000000),nPhotonsRequested,&ctrlc_caught,silentMode,D);
  }
  #endif

  *nPhotonsPtr = (double)O->nPhotons;
  *simTimePtr = (getMicroSeconds() - simulationTimeStart)/60000000.0; // In minutes
  if(!silentMode) {
    if(simulationTimed) printf("\nSimulated %0.2e photons at a rate of %0.2e photons per minute\n",*nPhotonsPtr, *nPhotonsPtr/(*simTimePtr));
    else printf("\nSimulated for %0.2e minutes at a rate of %0.2e photons per minute\n",*simTimePtr, *nPhotonsPtr/(*simTimePtr));
    mexEvalString("drawnow;");
  }
  normalizeDeposition(B,G,LC,O); // Convert data to relative fluence rate
  free(G->M);
  free(B->S);
  free(smallArrays);
  
  if(Pa->nExamplePaths) {
    mxDestroyArray(mxGetField(MCout,0,"examplePaths"));
    mxSetField(MCout,0,"examplePaths",mxCreateNumericMatrix(4,Pa->pathsElems,mxDOUBLE_CLASS,mxREAL));
    for(idx=0;idx<4*Pa->pathsElems;idx++) mxGetPr(mxGetField(MCout,0,"examplePaths"))[idx] = Pa->data[idx];
    free(Pa->data);
  }
  
//   printf("\nDebug: %.18e %.18e %.18e %llu %llu %llu\n",D->dbls[0],D->dbls[1],D->dbls[2],D->ulls[0],D->ulls[1],D->ulls[2]);
}
