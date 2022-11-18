/********************************************
 *
 * MCmatlab.c, in the C programming language, written for MATLAB MEX function generation
 * C script for Monte Carlo Simulation of Photon Transport in 3D
 * 
 * Copyright 2017, 2018 by Anders K. Hansen
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
 * "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' -outdir +MCmatlab\@model\private .\+MCmatlab\src\MCmatlab.c .\+MCmatlab\src\libut.lib"
 * ... or the Microsoft Visual C++ compiler (MSVC) with
 * "mex COMPFLAGS='/Zp8 /GR /EHs /nologo /MD /openmp /W4 /WX /wd4204 /wd4100' -outdir +MCmatlab\@model\private .\+MCmatlab\src\MCmatlab.c .\+MCmatlab\src\libut.lib"
 * In my experience, GCC produces faster machine code than MSVC.
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Go to MATLAB's addon manager and tell it to install the "Support for MinGW-w64 compiler"
 * 2. Type "mex -setup" in the MATLAB command window and ensure that MATLAB has set the C compiler to MinGW64
 * 3. mex should now be able to compile the code using the above GCC command
 *
 * For C++ compiling:
 * "copyfile ./+MCmatlab/src/MCmatlab.c ./+MCmatlab/src/MCmatlab.cpp; mex COMPFLAGS='/Zp8 /GR /EHs /nologo /MD /openmp /W4 /WX /wd4204 /wd4100' -outdir +MCmatlab\@model\private .\+MCmatlab\src\MCmatlab.cpp .\+MCmatlab\src\libut.lib"
 * 
 * The source code in this file is written is such a way that it is
 * compilable by either C or C++ compilers, either with GCC, MSVC or
 * the Nvidia CUDA compiler called NVCC, which is based on MSVC. To
 * compile with CUDA GPU acceleration support, you must have MSVC
 * installed. As of January 2020, mexcuda does not work with MSVC 2019,
 * so I'd recommend MSVC 2017. You also need the Parallel Computing
 * Toolbox, which you will find in the MATLAB addon manager. To compile, run:
 * "copyfile ./+MCmatlab/src/MCmatlab.c ./+MCmatlab/src/MCmatlab_CUDA.cu; mexcuda -llibut COMPFLAGS='-use_fast_math -res-usage $COMPFLAGS' -outdir +MCmatlab\@model\private .\+MCmatlab\src\MCmatlab_CUDA.cu"
 * 
 ** COMPILING ON MAC
 * As of November 2022, the macOS compiler doesn't support libut (for ctrl+c 
 * breaking) or openmp (for multithreading).
 * Compile in MATLAB with "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -std=c11 -Wall' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -std=c11 -Wall' -outdir +MCmatlab/@model/private ./+MCmatlab/src/MCmatlab.c"
 * You can enable openmp and multithreading if you're willing to install a custom version of llvm: Check "enable-openmp-on-macos.txt"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Install XCode from the App Store and, subsequently, the Apple Command Line Tools: fire up a terminal and type "xcode-select --install"
 * 2. Type "mex -setup" in the MATLAB command window
 *
 ** Compiling on Linux
 * "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall' -outdir +MCmatlab/@model/private ./+MCmatlab/src/MCmatlab.c ./+MCmatlab/src/libut.so"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Use a package manager like apt to install GCC (on Ubuntu, part of the build-essential package)
 * 2. Type "mex -setup" in the MATLAB command window
 ********************************************/
// printf("Reached line %d...\n",__LINE__);mexEvalString("drawnow; pause(.005);");mexEvalString("drawnow; pause(.005);");mexEvalString("drawnow; pause(.005);"); // For inserting into code for debugging purposes

#include "mex.h"
#include "print.h"
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
  #define GPUHEAPMEMORYLIMIT 500000000 // Bytes to allow the GPU to allocate, mostly for tracking photons prior to storage in NFRdet. Must be less than the available memory.
  #define KERNELTIME 50000 // Microseconds to spend in each kernel call (watchdog timer is said to termine kernel function that take more than a few seconds)
  #include <curand_kernel.h>
  typedef curandStateXORWOW_t PRNG_t;
//   typedef curandStateMRG32k3a_t PRNG_t;
//   typedef curandStatePhilox4_32_10_t PRNG_t;
  #define THREADNUM (threadIdx.x + blockDim.x*blockIdx.x)
  #define ISFINITE(x) (isfinite(x))
  #define ISNAN(x) (isnan(x))
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
  #define ISFINITE(x) (mxIsFinite(x))
  #define ISNAN(x) (mxIsNaN(x))
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
#define CDFSIZE     201 // Each custom phase function CDF contains this many elements (has to be one plus the value set in getOpticalMediaProperties.m)

#include "MCmatlablib.c"

#ifdef __NVCC__ // If compiling for CUDA
__global__
#endif
void threadInitAndLoop(struct source *B_global, struct geometry *G_global,
          struct lightCollector *LC_global, struct paths *Pa, struct outputs *O_global, struct depositionCriteria *DC_global, long nM, size_t size_smallArrays,
          long long simulationTimeStart, long long microSecondsOrGPUCycles, unsigned long long nPhotonsRequested,
          int iL, int nL,
          bool *abortingPtr, bool silentMode, struct debug *D) {
  struct photon P_var;
  struct photon *P = &P_var;

  P->recordSize    = O_global->NFRdet? INITIALRECORDSIZE: 0; // If we're supposed to calculate NFRdet, we start the record at a size of 1000 elements - it will be dynamically expanded later if needed
  P->j_record      = O_global->NFRdet? (long *)malloc(P->recordSize*sizeof(long)): NULL;
  P->weight_record = O_global->NFRdet? (FLOATORDBL *)malloc(P->recordSize*sizeof(FLOATORDBL)): NULL;

  #ifdef __NVCC__ // If compiling for CUDA
  // Copy structs from global device memory to shared device memory, which is orders of magnitude faster since it is on-chip
  __shared__ struct source B_var;
  struct source *B = &B_var;
  __shared__ struct geometry G_var;
  struct geometry *G = &G_var;
  __shared__ struct lightCollector LC_var;
  struct lightCollector *LC = &LC_var;
  __shared__ struct outputs O_var;
  struct outputs *O = &O_var;
  __shared__ struct depositionCriteria DC_var;
  struct depositionCriteria *DC = &DC_var;
  extern __shared__ FLOATORDBL smallArrays[]; // Dynamically allocated shared memory with size implicitly specified by the third cuda kernel launch parameter
  if(!threadIdx.x) { // Only let one thread per block do the copying
    // Copy contents of structs and smallArrays (not deep copies, pointers still point to global device memory)
    B_var = *B_global;
    G_var = *G_global;
    LC_var = *LC_global;
    O_var = *O_global;
    DC_var = *DC_global;
    memcpy(smallArrays,G->muav,size_smallArrays);
    // Set all array pointers to the correct new locations in the shared memory version of smallArrays
    long CDFarraySize = B->FPIDdist1? B->FPIDdist1 - G->CDFs: (FLOATORDBL *)G->CDFidxv - G->CDFs;
    G->muav = smallArrays;
    G->musv = G->muav + nM;
    G->gv = G->musv + nM;
    G->RIv = G->gv + nM;
    G->CDFs = G->RIv + nM;
    B->FPIDdist1 = G->CDFs + CDFarraySize;
    B->AIDdist1 = B->FPIDdist1 + B->L_FPID1;
    B->FPIDdist2 = B->AIDdist1 + B->L_AID1;
    B->AIDdist2 = B->FPIDdist2 + B->L_FPID2;
    G->CDFidxv = (unsigned char *)(B->AIDdist2 + B->L_AID2);
  }
  __syncthreads(); // All threads in the block wait for the copy to have finished
  
  // Initialize the PRNG and timing for the major loop
  curand_init(simulationTimeStart, THREADNUM, 0, &P->PRNGstate);
  long long tEnd = clock64() + microSecondsOrGPUCycles;

  // Launch major loop
  while(clock64() < tEnd && O->nPhotons + THREADNUM < nPhotonsRequested) { // "+ THREADNUM" ensures that we avoid race conditions that might launch more than nPhotonsRequested photons
  #else
  struct source *B = B_global;
  struct geometry *G = G_global;
  struct lightCollector *LC= LC_global;
  struct outputs *O = O_global;
  struct depositionCriteria *DC = DC_global;
  // Check for failed memory allocations and initialize the PRNG
  if(P->recordSize && !P->j_record) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");
  if(P->recordSize && !P->weight_record) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");
  dsfmt_init_gen_rand(&P->PRNGstate,(unsigned long)simulationTimeStart + THREADNUM); // Seed the photon's random number generator
  int pctProgressThisWavelength = 0;      // Simulation progress in percent
  int pctProgress = 0;
  // Launch major loop
  while(pctProgressThisWavelength < 100 && O->nPhotons + THREADNUM < nPhotonsRequested && !*abortingPtr) { // "+ THREADNUM" ensures that we avoid race conditions that might launch more than nPhotonsRequested photons
  #endif
    launchPhoton(P,B,G,Pa,DC,abortingPtr,D);
    if(P->alive) getNewVoxelProperties(P,G,D);
    if(P->alive) atomicAddWrapperULL(&O_global->nPhotons,1); // We have to store the photon number in the global memory so it's visible to all blocks

    while(P->alive) { // keep doing scattering events
      while(P->alive && P->stepLeft>0) { // keep propagating
        propagatePhoton(P,G,O,DC,Pa,D);
        if(!P->sameVoxel) {
          checkEscape(P,Pa,G,LC,O,DC); // photon may die here
          if(P->alive) getNewVoxelProperties(P,G,D);
        }
      }
      if(P->alive) checkRoulette(P);
      if(P->alive) scatterPhoton(P,G,Pa,DC,D);
    }

    #ifndef __NVCC__
      // Check progress
      int pctTimeProgressThisWavelength = (int)(100.0*(getMicroSeconds() - simulationTimeStart)/microSecondsOrGPUCycles);
      int pctPhotonsProgressThisWavelength = (int)(100.0*O->nPhotons/nPhotonsRequested);
      int pctTimeProgress = (int)(100.0*(iL + (double)(getMicroSeconds() - simulationTimeStart)/microSecondsOrGPUCycles)/nL);
      int pctPhotonsProgress = (int)(100.0*(iL + (double)O->nPhotons/nPhotonsRequested)/nL);

      #ifdef _OPENMP
      #pragma omp master
      #endif
      {
        // Check whether ctrl+c has been pressed
        #ifndef __clang__
        if(utIsInterruptPending()) {
          *abortingPtr = true;
          printf("\nCtrl+C detected, stopping.");
        }
        #endif
        if(O->nPhotons) { // If launches did not fail
          // Print out message about progress.
          int newPctProgress = max(pctTimeProgress,pctPhotonsProgress);
          if(newPctProgress != pctProgress && !silentMode && !*abortingPtr) {
            mexPrintf("\b\b\b\b\b\b\b\b\b%3.i%% done", newPctProgress<100? newPctProgress: 100);
            mexEvalString("drawnow; pause(.005);");
          }
        }
      }
      pctProgressThisWavelength = max(pctTimeProgressThisWavelength,pctPhotonsProgressThisWavelength);
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
  mxArray *MatlabMC = mxGetPropertyShared(prhs[0],0,simFluorescence? "FMC": "MC");
  
  bool silentMode = mxIsLogicalScalarTrue(mxGetPropertyShared(MatlabMC,0,"silentMode"));
  bool calcNFR    = mxIsLogicalScalarTrue(mxGetPropertyShared(MatlabMC,0,"calcNFR")); // Are we supposed to calculate the NFR matrix?
  bool calcNFRdet = mxIsLogicalScalarTrue(mxGetPropertyShared(MatlabMC,0,"calcNFRdet")); // Are we supposed to calculate the NFRdet matrix?

  // To know if we are simulating fluorescence, we check if a "sourceDistribution" field exists. If so, we will use it later in the beam definition.
  mxArray *MatlabLS = mxGetProperty(MatlabMC,0,"LS");
  double *S_PDF       = simFluorescence? (double *)mxGetData(mxGetProperty(MatlabMC,0,"sourceDistribution")): NULL; // Power emitted by the individual voxels per unit volume. Can be percieved as an unnormalized probability density function of the 3D source distribution

  // Variables for timekeeping and number of photons
  bool            simulationTimed = mxIsNaN(*mxGetPr(mxGetPropertyShared(MatlabMC,0,"nPhotonsRequested")));
  bool            aborting = false;
  double          simulationTimeRequested = simulationTimed? *mxGetPr(mxGetPropertyShared(MatlabMC,0,"simulationTimeRequested")): INFINITY;
  unsigned long long nPhotonsRequested = simulationTimed? ULLONG_MAX: (unsigned long long)*mxGetPr(mxGetPropertyShared(MatlabMC,0,"nPhotonsRequested"));
  double          nThreads = 1; // Will be updated later with the correct number

  if(!silentMode) {
    // Display progress indicator
    printf(simFluorescence?"-----------Fluorescence Monte Carlo Simulation-----------\n":
                           "-----------------Monte Carlo Simulation------------------\n");
    if(simulationTimed) printf("Simulation duration = %0.3f min\n",simulationTimeRequested);
    else                printf("Requested # of photons = %0.2e\n",(double)nPhotonsRequested);
  }

  // Find out how much total memory to allocate for the small arrays (the array that for GPUs would be stored in GPU shared memory)
  mxArray *mediaProperties = mxGetPropertyShared(MatlabMC,0,"mediaProperties");
  double *muaMATLAB    = mxGetPr(mxGetField(mediaProperties,0,"mua"));
  double *musMATLAB    = mxGetPr(mxGetField(mediaProperties,0,"mus"));
  double *gMATLAB      = mxGetPr(mxGetField(mediaProperties,0,"g"));
  double *nMATLAB      = mxGetPr(mxGetField(mediaProperties,0,"n"));
  double *CDFidxMATLAB = mxGetPr(mxGetField(mediaProperties,0,"CDFidx"));
  int nL = (int)mxGetN(mxGetField(mediaProperties,0,"mua")); // Number of wavelengths (Lambdas)
  int nM = (int)mxGetM(mxGetField(mediaProperties,0,"mua")); // Number of media

  double nPhotonsCumulative = 0;
  double nPhotonsDetectedCumulative = 0;
  double simulationTimeCumulative = 0;
  mwSize const *dimPtr = mxGetDimensions(mxGetPropertyShared(MatlabMC,0,"M"));
  int sourceType = S_PDF? -1: (int)*mxGetPr(mxGetPropertyShared(MatlabLS,0,"sourceType"));
  FLOATORDBL emitterLength = S_PDF? 0: (FLOATORDBL)*mxGetPr(mxGetPropertyShared(MatlabLS,0,"emitterLength"));
  mxArray *MatlabSourceFPID = S_PDF? 0: mxGetPropertyShared(MatlabLS,0,"FPID");
  mxArray *MatlabSourceAID = S_PDF? 0: mxGetPropertyShared(MatlabLS,0,"AID");
  long L_FPID1 = sourceType >= 4? (long)mxGetN(mxGetPropertyShared(MatlabSourceFPID,0,sourceType == 4? "radialDistr": "XDistr")): 0;
  long L_AID1 = sourceType >= 4? (long)mxGetN(mxGetPropertyShared(MatlabSourceAID,0,sourceType == 4? "radialDistr": "XDistr")): 0;
  long L_FPID2 = sourceType == 5? (long)mxGetN(mxGetPropertyShared(MatlabSourceFPID,0,"YDistr")): 0;
  long L_AID2 = sourceType == 5? (long)mxGetN(mxGetPropertyShared(MatlabSourceAID,0,"YDistr")): 0;

  // Geometry struct definition
  long L = (long)(dimPtr[0]*dimPtr[1]*dimPtr[2]); // Total number of voxels in cuboid
  mxArray *MatlabG         = mxGetPropertyShared(prhs[0],0,"G");
  struct geometry G_var;
  struct geometry *G = &G_var;
  G->d[0] = (FLOATORDBL)*mxGetPr(mxGetPropertyShared(MatlabG,0,"dx"));
  G->d[1] = (FLOATORDBL)*mxGetPr(mxGetPropertyShared(MatlabG,0,"dy"));
  G->d[2] = (FLOATORDBL)*mxGetPr(mxGetPropertyShared(MatlabG,0,"dz"));
  G->n[0] = (long)dimPtr[0];
  G->n[1] = (long)dimPtr[1];
  G->n[2] = (long)dimPtr[2];
  G->farFieldRes = (long)*mxGetPr(mxGetPropertyShared(MatlabMC,0,"farFieldRes"));
  G->boundaryType = (int)*mxGetPr(mxGetPropertyShared(MatlabMC,0,"boundaryType"));
  G->M = (unsigned char *)malloc(L*sizeof(unsigned char)); // M
  G->interfaceNormals = (float *)mxGetData(mxGetPropertyShared(MatlabMC,0,"interfaceNormals"));
  unsigned char *M_matlab = (unsigned char *)mxGetData(mxGetPropertyShared(MatlabMC,0,"M"));
  for(idx=0;idx<L;idx++) G->M[idx] = M_matlab[idx] - 1; // Convert from MATLAB 1-based indexing to C 0-based indexing

  mxArray *MatlabCDFs = mxGetPropertyShared(MatlabMC,0,"CDFs");
  long CDFarraySize = (long)mxGetNumberOfElements(MatlabCDFs);
  size_t size_smallArrays = (nM*4 + CDFarraySize + L_FPID1 + L_AID1 + L_FPID2 + L_AID2)*sizeof(FLOATORDBL) + nM*sizeof(unsigned char);
  char *smallArrays = (char *)malloc(size_smallArrays); // Because smallArrays contain different data types, we just use pointer to char (1 byte) here, and make it the correct types in the derived pointers
  G->muav = (FLOATORDBL *)smallArrays;
  G->musv = G->muav + nM;
  G->gv   = G->musv + nM;
  G->RIv  = G->gv + nM;
  G->CDFs = G->RIv + nM;
  for(idx=0;idx<CDFarraySize;idx++) G->CDFs[idx] = (FLOATORDBL)mxGetPr(MatlabCDFs)[idx];

  // Fill depositionCriteria struct
  mxArray *MatlabDC = mxGetProperty(MatlabMC,0,"depositionCriteria");
  struct depositionCriteria DC_var = {
    infCast(*mxGetPr(mxGetPropertyShared(MatlabDC,0,"minScatterings"))),
    infCast(*mxGetPr(mxGetPropertyShared(MatlabDC,0,"maxScatterings"))),
    infCast(*mxGetPr(mxGetPropertyShared(MatlabDC,0,"minRefractions"))),
    infCast(*mxGetPr(mxGetPropertyShared(MatlabDC,0,"maxRefractions"))),
    infCast(*mxGetPr(mxGetPropertyShared(MatlabDC,0,"minReflections"))),
    infCast(*mxGetPr(mxGetPropertyShared(MatlabDC,0,"maxReflections"))),
    infCast(*mxGetPr(mxGetPropertyShared(MatlabDC,0,"minInterfaceTransitions"))),
    infCast(*mxGetPr(mxGetPropertyShared(MatlabDC,0,"maxInterfaceTransitions")))
  };
  struct depositionCriteria *DC = &DC_var;

  // Light Collector struct definition
  mxArray *MatlabLC      = mxGetPropertyShared(MatlabMC,0,"LC");
  bool useLightCollector = mxIsLogicalScalarTrue(mxGetPropertyShared(MatlabMC,0,"useLightCollector"));
  
  FLOATORDBL theta     = (FLOATORDBL)*mxGetPr(mxGetPropertyShared(MatlabLC,0,"theta"));
  FLOATORDBL phi       = (FLOATORDBL)*mxGetPr(mxGetPropertyShared(MatlabLC,0,"phi"));
  FLOATORDBL f         = (FLOATORDBL)*mxGetPr(mxGetPropertyShared(MatlabLC,0,"f"));
  long   nTimeBins = S_PDF? 0: (long)*mxGetPr(mxGetPropertyShared(MatlabLC,0,"nTimeBins"));

  struct lightCollector LC_var = {
    {(FLOATORDBL)*mxGetPr(mxGetPropertyShared(MatlabLC,0,"x")) - (mxIsFinite(f)? f*SIN(theta)*COS(phi):0), // r field, xyz coordinates of center of light collector
     (FLOATORDBL)*mxGetPr(mxGetPropertyShared(MatlabLC,0,"y")) - (mxIsFinite(f)? f*SIN(theta)*SIN(phi):0),
     (FLOATORDBL)*mxGetPr(mxGetPropertyShared(MatlabLC,0,"z")) - (mxIsFinite(f)? f*COS(theta)         :0)},
    theta,
    phi,
    f,
    (FLOATORDBL)*mxGetPr(mxGetPropertyShared(MatlabLC,0,"diam")),
    (FLOATORDBL)*mxGetPr(mxGetPropertyShared(MatlabLC,0,mxIsFinite(f)?"fieldSize":"NA")),
    {(long)*mxGetPr(mxGetPropertyShared(MatlabLC,0,"res")), nTimeBins? nTimeBins+2: 1},
    (FLOATORDBL)(S_PDF? 0: *mxGetPr(mxGetPropertyShared(MatlabLC,0,"tStart"))),
    (FLOATORDBL)(S_PDF? 0: *mxGetPr(mxGetPropertyShared(MatlabLC,0,"tEnd")))
  };
  struct lightCollector *LC = &LC_var;

  // Paths definitions (example photon trajectories)
  struct paths Pa_var;
  struct paths *Pa = &Pa_var;
  Pa->nExamplePaths = (long)*mxGetPr(mxGetPropertyShared(MatlabMC,0,"nExamplePaths")); // How many photon path examples are we supposed to store?
  Pa->nExamplePhotonPathsStarted = 0;
  Pa->pathsSize = INITIALPATHSSIZE*Pa->nExamplePaths; // Will be dynamically increased later if needed
  Pa->pathsElems = 0;
  Pa->data = Pa->nExamplePaths? (FLOATORDBL *)malloc(4*Pa->pathsSize*sizeof(FLOATORDBL)): NULL;
  if(Pa->pathsSize && !Pa->data) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");

  // Prepare output MATLAB arrays and the temporary struct to store output data in
  plhs[0] = mxDuplicateArray(prhs[0]);
  
  mxArray *MCout = mxGetPropertyShared(plhs[0],0,S_PDF? "FMC": "MC");
  mxArray *LCout = mxGetPropertyShared(MCout,0,"LC");
  
  mwSize outDimPtr[4] = {dimPtr[0], dimPtr[1], dimPtr[2], (mwSize)nL};
  if(calcNFR)           mxSetPropertyShared(MCout,0,"NFR",mxCreateNumericArray(4,outDimPtr,mxSINGLE_CLASS,mxREAL));
  if(calcNFRdet)        mxSetPropertyShared(MCout,0,"NFRdet",mxCreateNumericArray(4,outDimPtr,mxSINGLE_CLASS,mxREAL));
  if(useLightCollector) {
    mwSize LCsize[4] = {(mwSize)LC->res[0],(mwSize)LC->res[0],(mwSize)LC->res[1],(mwSize)nL};
    mxSetPropertyShared(LCout,0,"image", mxCreateNumericArray(4,LCsize,mxSINGLE_CLASS,mxREAL));
  }
  if(G->farFieldRes)    {
    mwSize FFsize[3] = {(mwSize)G->farFieldRes,(mwSize)G->farFieldRes,(mwSize)nL};
    mxSetPropertyShared(MCout,0,"farField", mxCreateNumericArray(3,FFsize,mxSINGLE_CLASS,mxREAL));
  }
  if(G->boundaryType == 1) {
    mwSize NI_xSize[3] = {(mwSize)G->n[1],(mwSize)G->n[2],(mwSize)nL};
    mwSize NI_ySize[3] = {(mwSize)G->n[0],(mwSize)G->n[2],(mwSize)nL};
    mwSize NI_zSize[3] = {(mwSize)G->n[0],(mwSize)G->n[1],(mwSize)nL};
    mxSetPropertyShared(MCout,0,"NI_xpos", mxCreateNumericArray(3,NI_xSize,mxSINGLE_CLASS,mxREAL));
    mxSetPropertyShared(MCout,0,"NI_xneg", mxCreateNumericArray(3,NI_xSize,mxSINGLE_CLASS,mxREAL));
    mxSetPropertyShared(MCout,0,"NI_ypos", mxCreateNumericArray(3,NI_ySize,mxSINGLE_CLASS,mxREAL));
    mxSetPropertyShared(MCout,0,"NI_yneg", mxCreateNumericArray(3,NI_ySize,mxSINGLE_CLASS,mxREAL));
    mxSetPropertyShared(MCout,0,"NI_zpos", mxCreateNumericArray(3,NI_zSize,mxSINGLE_CLASS,mxREAL));
    mxSetPropertyShared(MCout,0,"NI_zneg", mxCreateNumericArray(3,NI_zSize,mxSINGLE_CLASS,mxREAL));
  } else if(G->boundaryType == 2) {
    mwSize NI_zSize[3] = {(mwSize)(KILLRANGE*G->n[0]),(mwSize)(KILLRANGE*G->n[1]),(mwSize)nL};
    mxSetPropertyShared(MCout,0,"NI_zneg", mxCreateNumericArray(3,NI_zSize,mxSINGLE_CLASS,mxREAL));
  } else if(G->boundaryType == 3) {
    mwSize NI_zSize[3] = {(mwSize)G->n[0],(mwSize)G->n[1],(mwSize)nL};
    mxSetPropertyShared(MCout,0,"NI_zpos", mxCreateNumericArray(3,NI_zSize,mxSINGLE_CLASS,mxREAL));
    mxSetPropertyShared(MCout,0,"NI_zneg", mxCreateNumericArray(3,NI_zSize,mxSINGLE_CLASS,mxREAL));
  }

  struct MATLABoutputs O_MATLAB_var = {
    calcNFR? (float *)mxGetPr(mxGetPropertyShared(MCout,0,"NFR")): NULL,
    calcNFRdet? (float *)mxGetPr(mxGetPropertyShared(MCout,0,"NFRdet")): NULL,
    useLightCollector? (float *)mxGetPr(mxGetPropertyShared(LCout,0,"image")): NULL,
    G->farFieldRes? (float *)mxGetPr(mxGetPropertyShared(MCout,0,"farField")): NULL,
    G->boundaryType == 1? (float *)mxGetPr(mxGetPropertyShared(MCout,0,"NI_xpos")): NULL,
    G->boundaryType == 1? (float *)mxGetPr(mxGetPropertyShared(MCout,0,"NI_xneg")): NULL,
    G->boundaryType == 1? (float *)mxGetPr(mxGetPropertyShared(MCout,0,"NI_ypos")): NULL,
    G->boundaryType == 1? (float *)mxGetPr(mxGetPropertyShared(MCout,0,"NI_yneg")): NULL,
    G->boundaryType == 1 || G->boundaryType == 3?
                          (float *)mxGetPr(mxGetPropertyShared(MCout,0,"NI_zpos")): NULL,
    G->boundaryType != 0? (float *)mxGetPr(mxGetPropertyShared(MCout,0,"NI_zneg")): NULL
  };
  struct MATLABoutputs *O_MATLAB = &O_MATLAB_var;
  
  struct outputs O_var = {
    0, // nPhotons
    0, // nPhotonsDetected
    calcNFR? (FLOATORDBL *)calloc(G->n[0]*G->n[1]*G->n[2],sizeof(FLOATORDBL)): NULL,
    calcNFRdet? (FLOATORDBL *)calloc(G->n[0]*G->n[1]*G->n[2],sizeof(FLOATORDBL)): NULL,
    useLightCollector? (FLOATORDBL *)calloc(LC->res[0]*LC->res[0]*LC->res[1],sizeof(FLOATORDBL)): NULL,
    G->farFieldRes? (FLOATORDBL *)calloc(G->farFieldRes*G->farFieldRes,sizeof(FLOATORDBL)): NULL,
    G->boundaryType == 1? (FLOATORDBL *)calloc(G->n[1]*G->n[2],sizeof(FLOATORDBL)): NULL,
    G->boundaryType == 1? (FLOATORDBL *)calloc(G->n[1]*G->n[2],sizeof(FLOATORDBL)): NULL,
    G->boundaryType == 1? (FLOATORDBL *)calloc(G->n[0]*G->n[2],sizeof(FLOATORDBL)): NULL,
    G->boundaryType == 1? (FLOATORDBL *)calloc(G->n[0]*G->n[2],sizeof(FLOATORDBL)): NULL,
    G->boundaryType == 1 || G->boundaryType == 3?
                          (FLOATORDBL *)calloc(G->n[0]*G->n[1],sizeof(FLOATORDBL)): NULL,
    G->boundaryType != 0? (FLOATORDBL *)calloc(G->n[0]*G->n[1]*(G->boundaryType == 2? KILLRANGE*KILLRANGE: 1),sizeof(FLOATORDBL)): NULL
  };
  struct outputs *O = &O_var;

  // Beam struct definition
  FLOATORDBL power        = 0;
  FLOATORDBL *S           = NULL; // Cumulative distribution function
  if(S_PDF) {
    S = (FLOATORDBL *)malloc((L+1)*sizeof(FLOATORDBL));
    if(!S) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");
    S[0] = 0;
  }

  FLOATORDBL *FPIDdist1 = G->CDFs + CDFarraySize;
  FLOATORDBL *AIDdist1 = FPIDdist1 + L_FPID1;
  FLOATORDBL *FPIDdist2 = AIDdist1 + L_AID1;
  FLOATORDBL *AIDdist2 = FPIDdist2 + L_FPID2;
  if(sourceType >= 4) {
    double *MatlabFPIDdist1 = mxGetPr(mxGetPropertyShared(MatlabSourceFPID,0,sourceType == 4? "radialDistr": "XDistr"));
    if(L_FPID1 == 1) {
      *FPIDdist1 = -1 - (FLOATORDBL)MatlabFPIDdist1[0];
    } else {
      FPIDdist1[0] = 0;
      for(idx=1;idx<L_FPID1;idx++) FPIDdist1[idx] = FPIDdist1[idx-1] + (sourceType == 4? idx-1: 1)*(FLOATORDBL)MatlabFPIDdist1[idx-1] + (sourceType == 4? idx: 1)*(FLOATORDBL)MatlabFPIDdist1[idx];
      for(idx=1;idx<L_FPID1;idx++) FPIDdist1[idx] /= FPIDdist1[L_FPID1-1];
    }
    
    double *MatlabAIDdist1 = mxGetPr(mxGetPropertyShared(MatlabSourceAID,0,sourceType == 4? "radialDistr": "XDistr"));
    if(L_AID1 == 1) {
      *AIDdist1 = -1 - (FLOATORDBL)MatlabAIDdist1[0];
    } else {
      AIDdist1[0] = 0;
      for(idx=1;idx<L_AID1;idx++) AIDdist1[idx] = AIDdist1[idx-1] + (sourceType == 4? idx-1: 1)*(FLOATORDBL)MatlabAIDdist1[idx-1] + (sourceType == 4? idx: 1)*(FLOATORDBL)MatlabAIDdist1[idx];
      for(idx=1;idx<L_AID1;idx++) AIDdist1[idx] /= AIDdist1[L_AID1-1];
    }
    
    if(sourceType == 5) {
      double *MatlabFPIDdist2 = mxGetPr(mxGetPropertyShared(MatlabSourceFPID,0,"YDistr"));
      if(L_FPID2 == 1) {
        *FPIDdist2 = -1 - (FLOATORDBL)MatlabFPIDdist2[0];
      } else {
        FPIDdist2[0] = 0;
        for(idx=1;idx<L_FPID2;idx++) FPIDdist2[idx] = FPIDdist2[idx-1] + (FLOATORDBL)MatlabFPIDdist2[idx-1] + (FLOATORDBL)MatlabFPIDdist2[idx];
        for(idx=1;idx<L_FPID2;idx++) FPIDdist2[idx] /= FPIDdist2[L_FPID2-1];
      }
    
      double *MatlabAIDdist2 = mxGetPr(mxGetPropertyShared(MatlabSourceAID,0,"YDistr"));
      if(L_AID2 == 1) {
        *AIDdist2 = -1 - (FLOATORDBL)MatlabAIDdist2[0];
      } else {
        AIDdist2[0] = 0;
        for(idx=1;idx<L_AID2;idx++) AIDdist2[idx] = AIDdist2[idx-1] + (FLOATORDBL)MatlabAIDdist2[idx-1] + (FLOATORDBL)MatlabAIDdist2[idx];
        for(idx=1;idx<L_AID2;idx++) AIDdist2[idx] /= AIDdist2[L_AID2-1];
      }
    }
  }
  G->CDFidxv = (unsigned char *)(AIDdist2 + L_AID2); // Array that describes which CDF array is the one that applies to a particular medium

  FLOATORDBL tb = (FLOATORDBL)(S? 0: *mxGetPr(mxGetPropertyShared(MatlabLS,0,"theta")));
  FLOATORDBL pb = (FLOATORDBL)(S? 0: *mxGetPr(mxGetPropertyShared(MatlabLS,0,"phi")));
  FLOATORDBL psib = (FLOATORDBL)(S? 0: *mxGetPr(mxGetPropertyShared(MatlabLS,0,"psi")));

  FLOATORDBL u[3] = {SIN(tb)*COS(pb),SIN(tb)*SIN(pb),COS(tb)}; // Temporary array
  
  FLOATORDBL w[3] = {1,0,0};
  FLOATORDBL v[3] = {0,0,1};
  if(u[2]!=1) unitcrossprod(u,v,w);
  axisrotate(w,u,psib,v);
  unitcrossprod(u,v,w);
  
  struct source B_var = {
    S? 0: sourceType,
    S? 0: emitterLength,
    S? 0: FPIDdist1,
    S? 0: L_FPID1,
    S? 0: (FLOATORDBL)(sourceType == 5? *mxGetPr(mxGetPropertyShared(MatlabSourceFPID,0,"XWidth")): *mxGetPr(mxGetPropertyShared(MatlabSourceFPID,0,"radialWidth"))),
    S? 0: FPIDdist2,
    S? 0: L_FPID2,
    S? 0: (FLOATORDBL)(sourceType == 5? *mxGetPr(mxGetPropertyShared(MatlabSourceFPID,0,"YWidth")): 0),
    S? 0: AIDdist1,
    S? 0: L_AID1,
    S? 0: (FLOATORDBL)(sourceType == 5? *mxGetPr(mxGetPropertyShared(MatlabSourceAID,0,"XWidth")): *mxGetPr(mxGetPropertyShared(MatlabSourceAID,0,"radialWidth"))),
    S? 0: AIDdist2,
    S? 0: L_AID2,
    S? 0: (FLOATORDBL)(sourceType == 5? *mxGetPr(mxGetPropertyShared(MatlabSourceAID,0,"YWidth")): 0),
    S,
    power,
    {(FLOATORDBL)(S? 0: *mxGetPr(mxGetPropertyShared(MatlabLS,0,"xFocus"))),
     (FLOATORDBL)(S? 0: *mxGetPr(mxGetPropertyShared(MatlabLS,0,"yFocus"))),
     (FLOATORDBL)(S? 0: *mxGetPr(mxGetPropertyShared(MatlabLS,0,"zFocus")))},
    {u[0],u[1],u[2]},
    {v[0],v[1],v[2]}, // normal vector to beam center axis
    {w[0],w[1],w[2]}
  };
  struct source *B = &B_var;

  for(int iL = 0; iL < nL && !aborting; iL++) {
    for(idx=0;idx<nM;idx++) {
      G->muav[idx]    = (FLOATORDBL)      muaMATLAB[idx + iL*nM];
      G->musv[idx]    = (FLOATORDBL)      musMATLAB[idx + iL*nM];
      G->gv[idx]      = (FLOATORDBL)        gMATLAB[idx + iL*nM];
      G->RIv[idx]     = (FLOATORDBL)        nMATLAB[idx + iL*nM];
      G->CDFidxv[idx] = (unsigned char)CDFidxMATLAB[idx + iL*nM];
    }

    if(S) {
      for(idx=1;idx<(L+1);idx++) S[idx] = S[idx-1] + (FLOATORDBL)S_PDF[iL*L + idx-1];
      B->power        = (FLOATORDBL)(S[L]*G->d[0]*G->d[1]*G->d[2]);
      for(idx=1;idx<(L+1);idx++) S[idx] /= S[L];
    } else {
      B->power        = (FLOATORDBL)mxGetPr(mxGetPropertyShared(MatlabMC,0,"spectrum"))[iL];
    }

    // ============================ MAJOR CYCLE ========================
    unsigned long long nPhotonsRequested_ThisWavelength = simulationTimed? ULLONG_MAX: (unsigned long long)(iL == nL - 1? nPhotonsRequested - nPhotonsCumulative: nPhotonsRequested/nL);
    double simulationTimeRequested_ThisWavelength = simulationTimeRequested/nL;
    #ifdef __NVCC__ // If compiling for CUDA
    long GPUdevice = *(long *)mxGetPr(mxGetPropertyShared(MatlabMC,0,"GPUdevice"));
    gpuErrchk(cudaSetDevice(GPUdevice));
  
    int threadsPerBlock, blocks; gpuErrchk(cudaOccupancyMaxPotentialBlockSize(&blocks,&threadsPerBlock,&threadInitAndLoop,size_smallArrays,0));
    nThreads = threadsPerBlock*blocks;
  
    struct cudaDeviceProp CDP; gpuErrchk(cudaGetDeviceProperties(&CDP,GPUdevice));
  //   if(!silentMode) printf("Using %s with CUDA compute capability %d.%d,\n    launching %d blocks, each with %d threads,\n    using %d/%d bytes of shared memory per block\n",CDP.name,CDP.major,CDP.minor,blocks,threadsPerBlock,384+size_smallArrays,CDP.sharedMemPerBlock);
    if(!silentMode && iL == 0) printf("Using device %d: %s with CUDA compute capability %d.%d\n",GPUdevice,CDP.name,CDP.major,CDP.minor);
  
    size_t heapSizeLimit; gpuErrchk(cudaDeviceGetLimit(&heapSizeLimit,cudaLimitMallocHeapSize));
    if(calcNFRdet && heapSizeLimit < GPUHEAPMEMORYLIMIT) {
      gpuErrchk(cudaDeviceReset());
      gpuErrchk(cudaDeviceSetLimit(cudaLimitMallocHeapSize,GPUHEAPMEMORYLIMIT));
    }
  
    int clock; gpuErrchk(cudaDeviceGetAttribute(&clock,cudaDevAttrClockRate,GPUdevice));
  
    struct geometry *G_dev;
    struct source *B_dev;
    struct lightCollector *LC_dev;
    struct paths *Pa_dev;
    struct outputs *O_dev;
    struct depositionCriteria *DC_dev;
    struct debug *D_dev;
    createDeviceStructs(G,&G_dev,B,&B_dev,LC,&LC_dev,Pa,&Pa_dev,O,&O_dev,DC,&DC_dev,nM,L,size_smallArrays,D,&D_dev);
  
    int pctProgress = 0;
    long long timeLeft = simulationTimed? (long long)(simulationTimeRequested_ThisWavelength*60000000): LLONG_MAX; // microseconds
  
    if(!silentMode && iL == 0) {
      printf("Calculating...   0%% done");
      mexEvalString("drawnow; pause(.005);");
    }
    long long simulationTimeStart = getMicroSeconds();
    long long prevtime = simulationTimeStart;
    do {
      // Run kernel
      threadInitAndLoop<<<blocks, threadsPerBlock, size_smallArrays>>>(B_dev,G_dev,LC_dev,Pa_dev,O_dev,DC_dev,nM,size_smallArrays,prevtime,clock/1000*min((long long)KERNELTIME,timeLeft),nPhotonsRequested_ThisWavelength,0,0,NULL,false,D_dev);
      gpuErrchk(cudaPeekAtLastError());
      gpuErrchk(cudaDeviceSynchronize());
      // Progress indicator
      long long newtime = getMicroSeconds();
      gpuErrchk(cudaMemcpy(O, O_dev, sizeof(unsigned long long),cudaMemcpyDeviceToHost)); // Copy just nPhotons, the first 8 bytes of O
      if(simulationTimed) {
        timeLeft = (long long)(simulationTimeRequested_ThisWavelength*60000000) - newtime + simulationTimeStart; // In microseconds
        pctProgress = (int)(100.0*(iL + 1.0 - timeLeft/(simulationTimeRequested_ThisWavelength*60000000.0))/nL);
      } else {
        pctProgress = (int)(100.0*(iL + O->nPhotons/nPhotonsRequested_ThisWavelength)/nL);
      }
  
      // Check for ctrl+c
      if(utIsInterruptPending()) {
        aborting = true;
        printf("\nCtrl+C detected, stopping.");
      }
  
      if(!silentMode) {
        printf("\b\b\b\b\b\b\b\b\b%3d%% done",pctProgress);
        mexEvalString("drawnow; pause(.005);");
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
    } while(timeLeft > 0 && O->nPhotons < nPhotonsRequested_ThisWavelength && !aborting && O->nPhotons); // The O->nPhotons is there to stop looping if launches failed
    if(!silentMode && O->nPhotons) printf("\b\b\b\b\b\b\b\b\b%3d%% done",pctProgress);
  
    retrieveAndFreeDeviceStructs(G,G_dev,B,B_dev,LC,LC_dev,Pa,Pa_dev,O,O_dev,DC_dev,L,D,D_dev);
  
    #else
  
    if(!silentMode && iL == 0) {
      printf("Calculating...   0%% done");
      mexEvalString("drawnow; pause(.005);");
    }
    long long simulationTimeStart = getMicroSeconds();
    #ifdef _OPENMP
    bool useAllCPUs = mxIsLogicalScalarTrue(mxGetPropertyShared(MatlabMC,0,"useAllCPUs"));
    nThreads = useAllCPUs? omp_get_num_procs(): max(omp_get_num_procs()-1,1);
    #pragma omp parallel num_threads((long)nThreads)
    #else
    nThreads = 1;
    #endif
    {
      threadInitAndLoop(B,G,LC,Pa,O,DC,nM,0,simulationTimeStart,(long long)(simulationTimeRequested_ThisWavelength*60000000),nPhotonsRequested_ThisWavelength,iL,nL,&aborting,silentMode,D);
    }
    #endif
  
    double nPhotons = (double)O->nPhotons;
    nPhotonsCumulative += nPhotons;
    double nPhotonsDetected = (double)O->nPhotonsDetected;
    nPhotonsDetectedCumulative += nPhotonsDetected;
    double simTime = (getMicroSeconds() - simulationTimeStart)/60000000.0; // In minutes
    simulationTimeCumulative += simTime;
    if(!silentMode) {
      if(!nPhotons) printf("\nERROR: All photons launch outside simulation cuboid. Check your model definition.\n");
      else if(iL == nL - 1) {
        printf("\b\b\b\b\b\b\b\b\b100%% done");
        printf("\nSimulated %0.2e photons over %0.2e minutes at a rate of %0.2e photons per minute\n",nPhotonsCumulative, simulationTimeCumulative, nPhotonsCumulative/simulationTimeCumulative);
      }
      mexEvalString("drawnow; pause(.005);");
    }
    normalizeDepositionAndResetO(B,G,LC,O,O_MATLAB,iL,B->power); // Convert data to relative fluence rate and save in O_MATLAB
  }

  mxArray *output = mxCreateDoubleMatrix(1,1,mxREAL);
  *mxGetPr(output) = nPhotonsCumulative;
  mxSetProperty(MCout,0,"nPhotons",output);
  *mxGetPr(output) = nPhotonsDetectedCumulative;
  mxSetProperty(MCout,0,"nPhotonsDetected",output);
  *mxGetPr(output) = nThreads;
  mxSetProperty(MCout,0,"nThreads",output);
  *mxGetPr(output) = simulationTimeCumulative;
  mxSetProperty(MCout,0,"simulationTime",output);
  mxDestroyArray(output);
  if(LC->res[0]*LC->res[1]*nL == 1) {
    output = mxCreateNumericMatrix(1,1,mxSINGLE_CLASS,mxREAL);
    *(float *)mxGetPr(output) = *O_MATLAB->image;
    mxSetProperty(LCout,0,"image",output);
    mxDestroyArray(output);
  }

  if(Pa->nExamplePaths) {
    mxSetPropertyShared(MCout,0,"examplePaths",mxCreateNumericMatrix(4,Pa->pathsElems,mxDOUBLE_CLASS,mxREAL));
    for(idx=0;idx<4*Pa->pathsElems;idx++) mxGetPr(mxGetPropertyShared(MCout,0,"examplePaths"))[idx] = Pa->data[idx];
    free(Pa->data);
  }

  free(B->S);
  free(smallArrays);
  free(G->M);
  free(O->NFR);
  free(O->NFRdet);
  free(O->image);
  free(O->FF);
  free(O->NI_xpos);
  free(O->NI_xneg);
  free(O->NI_ypos);
  free(O->NI_yneg);
  free(O->NI_zpos);
  free(O->NI_zneg);
//   printf("\nDebug: %.18e %.18e %.18e %llu %llu %llu\n",D->dbls[0],D->dbls[1],D->dbls[2],D->ulls[0],D->ulls[1],D->ulls[2]);
}
