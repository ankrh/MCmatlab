/********************************************
 *
 * getDownsampledParamVals.c, in the C programming language, written for MATLAB MEX function generation
 * C script for Monte Carlo Simulation of Photon Transport in 3D
 * 
 * Copyright 2017, 2018 by Anders K. Hansen, DTU Fotonik
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
 * Can be compiled in MATLAB with "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall -pedantic' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall -pedantic' -outdir helperfuncs\private .\src\getDownsampledParamVals.c"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Go to MATLAB's addon manager and tell it to install the "Support for MinGW-w64 compiler"
 * 2. Type "mex -setup" in the MATLAB command window and ensure that MATLAB has set the C compiler to MinGW64
 * 3. mex should now be able to compile the code using the above command
 *
 ** COMPILING ON MAC
 * As of June 2017, the macOS compiler doesn't support libut (for ctrl+c 
 * breaking) or openmp (for multithreading).
 * Compile in MATLAB with "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -std=c11 -Wall -pedantic' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -std=c11 -Wall -pedantic' -outdir helperfuncs/private ./src/getDownsampledParamVals.c"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Install XCode from the App Store
 * 2. Type "mex -setup" in the MATLAB command window
 ********************************************/

#include "mex.h"
#include <math.h>
#ifdef _WIN32 // This is defined on both win32 and win64 systems. We use this preprocessor condition to avoid loading openmp or libut on, e.g., Mac
#include <omp.h>
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
  long N1 = mxGetN(prhs[0]); // Number of parameter sets provided as input
  long N2 = *mxGetPr(prhs[1]); // Requested number of downsampled parameter sets
  float *paramVals = mxGetData(prhs[0]);
  long i; // General-purpose index
  plhs[0] = mxCreateNumericMatrix(3,N2,mxSINGLE_CLASS,mxREAL);
  float *downsampledParamVals = mxGetData(plhs[0]);
  plhs[1] = mxCreateNumericMatrix(1,N1,mxUINT8_CLASS,mxREAL);
  unsigned char *binidx = mxGetData(plhs[1]);
  
  float max_mua = paramVals[0];
  float min_mua = paramVals[0];
  float sum_mua = paramVals[0];
  float max_mus = paramVals[1];
  float min_mus = paramVals[1];
  float sum_mus = paramVals[1];
  float max_g   = paramVals[2];
  float min_g   = paramVals[2];
  float sum_g   = paramVals[2];
  for(i=1;i<N1;i++) { // i is the current column
    max_mua  = paramVals[  3*i] > max_mua? paramVals[  3*i]: max_mua;
    min_mua  = paramVals[  3*i] < min_mua? paramVals[  3*i]: min_mua;
    sum_mua += paramVals[  3*i];
    max_mus  = paramVals[1+3*i] > max_mus? paramVals[1+3*i]: max_mus;
    min_mus  = paramVals[1+3*i] < min_mus? paramVals[1+3*i]: min_mus;
    sum_mus += paramVals[1+3*i];
    max_g    = paramVals[2+3*i] > max_g?   paramVals[2+3*i]: max_g;
    min_g    = paramVals[2+3*i] < min_g?   paramVals[2+3*i]: min_g;
    sum_g   += paramVals[2+3*i];
  }
  bool mua_varying = max_mua != min_mua;
  bool mus_varying = max_mus != min_mus;
  bool g_varying   = max_g   != min_g;

  long mus_idx = mua_varying; // May be 0 or 1
  long g_idx   = mua_varying + mus_varying; // May be 0, 1 or 2
  
  int rows = mua_varying + mus_varying + g_varying; // How many rows we need in normParamVals
  float *normParamVals = malloc(rows*N1*sizeof(float));
  for(i=0;i<N1;i++) {
    if(mua_varying) normParamVals[          rows*i] = paramVals[  3*i]/(max_mua - min_mua);
    if(mus_varying) normParamVals[mus_idx + rows*i] = paramVals[1+3*i]/(max_mus - min_mus);
    if(g_varying  ) normParamVals[g_idx   + rows*i] = paramVals[2+3*i]/(max_g   - min_g  );
  }
  
  long maxDistIdx = 0;
  float *sqrDists = calloc(N1,sizeof(float));
  for(i=0;i<N1;i++) {
    if(mua_varying) sqrDists[i] += pow(normParamVals[          rows*i] - sum_mua/(max_mua - min_mua)/N1,2);
    if(mus_varying) sqrDists[i] += pow(normParamVals[mus_idx + rows*i] - sum_mus/(max_mus - min_mus)/N1,2);
    if(g_varying  ) sqrDists[i] += pow(normParamVals[g_idx   + rows*i] - sum_g  /(max_g   - min_g  )/N1,2);
    
    if(sqrDists[i] > sqrDists[maxDistIdx]) maxDistIdx = i;
  }
  
  float *normDownsampledParamVals = malloc(rows*N2*sizeof(float));
  // Write first downsampled point
  if(mua_varying) normDownsampledParamVals[0      ] = normParamVals[          rows*maxDistIdx];
  if(mus_varying) normDownsampledParamVals[mus_idx] = normParamVals[mus_idx + rows*maxDistIdx];
  if(g_varying)   normDownsampledParamVals[g_idx  ] = normParamVals[g_idx   + rows*maxDistIdx];
  
  // Find squared distances to the first downsampled point
  maxDistIdx = 0;
  for(i=0;i<N1;i++) {
    sqrDists[i] = 0;
    if(mua_varying) sqrDists[i] += pow(normParamVals[          rows*i] - normDownsampledParamVals[0      ],2);
    if(mus_varying) sqrDists[i] += pow(normParamVals[mus_idx + rows*i] - normDownsampledParamVals[mus_idx],2);
    if(g_varying)   sqrDists[i] += pow(normParamVals[g_idx   + rows*i] - normDownsampledParamVals[g_idx  ],2);
    
    binidx[i] = 1; // Using MATLAB indexing, starting from 1
    if(sqrDists[i] > sqrDists[maxDistIdx]) maxDistIdx = i;
  }
  
  // Loop over remaining points
  #ifdef _WIN32
  long nThreads = omp_get_num_procs() == 1? omp_get_num_procs(): omp_get_num_procs()-1;
  #else
  long nThreads = 1;
  #endif
  long *maxDistIdxs = calloc(nThreads,sizeof(long));
  
  #ifdef _WIN32
  #pragma omp parallel num_threads(nThreads)
  #endif
  {
    long j,k;
    float newSqrDist;
    
    for(j=1;j<N2;j++) {
      #ifdef _WIN32
      #pragma omp master
      #endif
      {
      // Write next downsampled point
      if(mua_varying) normDownsampledParamVals[          rows*j] = normParamVals[          rows*maxDistIdx];
      if(mus_varying) normDownsampledParamVals[mus_idx + rows*j] = normParamVals[mus_idx + rows*maxDistIdx];
      if(g_varying)   normDownsampledParamVals[g_idx   + rows*j] = normParamVals[g_idx   + rows*maxDistIdx];
      }
      #ifdef _WIN32
      #pragma omp barrier
      #endif
      
      // Find the distances all the points have to this new downsampled point
      #ifdef _WIN32
      #pragma omp for schedule(auto)
      #endif
      for(k=0;k<N1;k++) {
        newSqrDist = 0;
        if(mua_varying) newSqrDist += pow(normParamVals[          rows*k] - normDownsampledParamVals[          rows*j],2);
        if(mus_varying) newSqrDist += pow(normParamVals[mus_idx + rows*k] - normDownsampledParamVals[mus_idx + rows*j],2);
        if(g_varying)   newSqrDist += pow(normParamVals[g_idx   + rows*k] - normDownsampledParamVals[g_idx   + rows*j],2);

        if(newSqrDist < sqrDists[k]) {
          sqrDists[k] = newSqrDist;
          binidx[k] = j+1; // Using MATLAB indexing. The k'th input point is closest to the j'th downsampled point
        }

        #ifdef _WIN32
        if(sqrDists[k] > sqrDists[maxDistIdxs[omp_get_thread_num()]]) maxDistIdxs[omp_get_thread_num()] = k;
        #else
        if(sqrDists[k] > sqrDists[maxDistIdxs[0]]) maxDistIdxs[0] = k;
        #endif
      }

      #ifdef _WIN32
      #pragma omp master
      #endif
      {
      maxDistIdx = 0;
      for(i=0;i<nThreads;i++) if(sqrDists[maxDistIdxs[i]] > sqrDists[maxDistIdx]) maxDistIdx = maxDistIdxs[i];
      }
    #ifdef _WIN32
    #pragma omp barrier
    #endif
    }
  }
  
  // Rescale normDownsampledParamVals to get downsampledParamVals
  for(i=0;i<N2;i++) {
    downsampledParamVals[    3*i] = mua_varying? normDownsampledParamVals[          rows*i]*(max_mua - min_mua): paramVals[0];
    downsampledParamVals[1 + 3*i] = mus_varying? normDownsampledParamVals[mus_idx + rows*i]*(max_mus - min_mus): paramVals[1];
    downsampledParamVals[2 + 3*i] = g_varying?   normDownsampledParamVals[g_idx   + rows*i]*(max_g   - min_g  ): paramVals[2];
  }
  
  free(maxDistIdxs);
  free(normDownsampledParamVals);
  free(sqrDists);
  free(normParamVals);
}
