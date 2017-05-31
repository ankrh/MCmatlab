/*=================================================================
 *
 * C script for heat propagation based on Monte Carlo input
 * arguments: (nt,[[[Temp]]],[[[Tissue]]],[[[dQperdeltaT]]],[[[dQ_abs]]],[HC])
 * 
 * Can be compiled using "mex CFLAGS='$CFLAGS -fopenmp' LDFLAGS='$LDFLAGS -fopenmp' COPTIMFLAGS='$COPTIMFLAGS -fopenmp -Ofast' LDOPTIMFLAGS='$LDOPTIMFLAGS -fopenmp -Ofast' DEFINES='$DEFINES -fopenmp' propagator.c"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Go to MATLAB's addon manager and tell it to install the "Support for MinGW-w64 compiler"
 * 2. Type "mex -setup" in the MATLAB command line and ensure that MATLAB has set the C compiler to MinGW64
 * 3. Copy the files "libgomp.a" and "libgomp.spec" to the folder with a path similar to "C:\ProgramData\MATLAB\SupportPackages\R2017a\MW_MinGW_4_9\lib\gcc\x86_64-w64-mingw32\4.9.2"
 * 4. mex should now be able to compile the code using the above command but in order to run, it needs to have the file "libgomp_64-1.dll" copied to the same folder as the mex file.
 *=================================================================*/
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include "omp.h"

void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, mxArray const *prhs[] ) {
    
    int nx,ny,nz,nT;
    /* Check for proper number of arguments */
    
    if (nrhs != 6) {
        mexErrMsgIdAndTxt( "MATLAB:heatSim:invalidNumInputs",
                "Six input arguments required (nt,[[[Temp]]],[[[Tissue]]],[[[dQperdeltaT]]],[[[dQ_abs]]],[HC])");
    } else if (nlhs > 1) {
        mexErrMsgIdAndTxt( "MATLAB:heatSim:maxlhs",
                "Too many output arguments.");
    }
    
    mwSize const *dimPtr = mxGetDimensions(prhs[1]);
    nx = dimPtr[0];
    ny = dimPtr[1];
    nz = dimPtr[2];
    nT = mxGetN(prhs[5]); // Number of tissues in simulation
    
    plhs[0] = mxCreateNumericArray(3,dimPtr,mxDOUBLE_CLASS,mxREAL);
    
    double *outputTemp    = mxGetPr(plhs[0]); // newTemp is an nx*ny*nz array of doubles
    int nt                = (int)*mxGetPr(prhs[0]); // nt is an integer (number of time steps to perform)
    double *Temp          = mxGetPr(prhs[1]); // Temp is an nx*ny*nz array of doubles
    unsigned char *T      = (unsigned char *)mxGetData(prhs[2]); // T is a nx*ny*nz array of uint8 (unsigned char) containing values from 0..nT-1
    double *dQperdeltaT   = mxGetPr(prhs[3]); // dQperdeltaT is an nT*nT*3 array of doubles
    double *dQ_abs        = mxGetPr(prhs[4]); // dQ_abs is an nx*ny*nz array of doubles
    double *HC            = mxGetPr(prhs[5]); // HC is an nT array of doubles
    
    double *tempTemp = (double *)malloc(nx*ny*nz*sizeof(double)); // Temporary temperature matrix
    
    double **srcTemp = &Temp; // Pointer to the pointer to whichever temperature matrix is to be read from
    double **tgtTemp; // Pointer to the pointer to whichever temperature matrix is to be written to
    
    if(nt%2) { // If nt is odd then we can start by writing to the output temperature matrix (pointed to by plhs[0]), otherwise we start by writing to the temporary temperature matrix
        tgtTemp = &outputTemp;
    } else {
        tgtTemp = &tempTemp;
    }
    
//#pragma omp parallel num_threads(omp_get_num_procs())
#pragma omp parallel num_threads(8)
    {
        int ix,iy,iz;
        int n,i;
        double dQ;
        
        //printf("Thread %u iterating from %u to %u\r\n",omp_get_thread_num(),nx*ny*nz*omp_get_thread_num()/omp_get_num_threads(),nx*ny*nz*(omp_get_thread_num()+1)/omp_get_num_threads()-1);
        for(n=0;n<nt;n++) {
            for(i = nx*ny*nz*omp_get_thread_num()/omp_get_num_threads();i<nx*ny*nz*(omp_get_thread_num()+1)/omp_get_num_threads();i++) {
                ix = i%nx;    // ix = (i/(1    ))%(nx);
                iy = i/nx%ny; // iy = (i/(nx   ))%(ny);
                iz = i/nx/ny; // iz = (i/(nx*ny))%(nz);
                
                dQ = dQ_abs[i];
                if(ix != 0   ) {dQ += ((*srcTemp)[i-1]     - (*srcTemp)[i])*dQperdeltaT[T[i] + T[i-1    ]*nT          ];};
                if(ix != nx-1) {dQ += ((*srcTemp)[i+1]     - (*srcTemp)[i])*dQperdeltaT[T[i] + T[i+1    ]*nT          ];};
                if(iy != 0   ) {dQ += ((*srcTemp)[i-nx]    - (*srcTemp)[i])*dQperdeltaT[T[i] + T[i-nx   ]*nT + nT*nT  ];};
                if(iy != ny-1) {dQ += ((*srcTemp)[i+nx]    - (*srcTemp)[i])*dQperdeltaT[T[i] + T[i+nx   ]*nT + nT*nT  ];};
                if(iz != 0   ) {dQ += ((*srcTemp)[i-nx*ny] - (*srcTemp)[i])*dQperdeltaT[T[i] + T[i-nx*ny]*nT + 2*nT*nT];};
                if(iz != nz-1) {dQ += ((*srcTemp)[i+nx*ny] - (*srcTemp)[i])*dQperdeltaT[T[i] + T[i+nx*ny]*nT + 2*nT*nT];};
                (*tgtTemp)[i] = dQ/HC[T[i]] + (*srcTemp)[i];
            }
            #pragma omp barrier
            #pragma omp master
            {
                if(tgtTemp == &outputTemp) {
                    tgtTemp = &tempTemp;
                    srcTemp = &outputTemp;
                } else {
                    tgtTemp = &outputTemp;
                    srcTemp = &tempTemp;
                }
            }
            #pragma omp barrier
        }
    }
    
    
    free(tempTemp);
    return;
}
