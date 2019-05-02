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
 * Can be compiled using "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall -pedantic' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall -pedantic' -outdir helperfuncs\private .\src\finiteElementHeatPropagator.c ".\src\libut.lib""
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Go to MATLAB's addon manager and tell it to install the "Support for MinGW-w64 compiler"
 * 2. Type "mex -setup" in the MATLAB command window and ensure that MATLAB has set the C compiler to MinGW64
 * 3. mex should now be able to compile the code using the above command
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

#include <math.h>
#include "mex.h"
#include "expint.c" // For calculating the special function E1 (scaled by exp), adapted from the GNU Scientific Library
#ifdef _WIN32 // This is defined on both win32 and win64 systems. We use this preprocessor condition to avoid loading openmp or libut on, e.g., Mac
#include "omp.h"
extern bool utIsInterruptPending(); // Allows catching ctrl+c while executing the mex function
#endif

#define R 8.3144598 // Gas constant [J/mol/K]
#define CELSIUSZERO 273.15
#define T1 (*srcTemp)
#define T2 (*tgtTemp)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    bool ctrlc_caught = false;           // Has a ctrl+c been passed from MATLAB?
    long nx,ny,nz,nM;
    
    mwSize const *dimPtr = mxGetDimensions(prhs[0]);
    plhs[0] = mxCreateNumericArray(3,dimPtr,mxSINGLE_CLASS,mxREAL);
    nx = dimPtr[0];
    ny = dimPtr[1];
    nz = dimPtr[2];
    nM = mxGetM(mxGetField(prhs[2],0,"dTdtperdeltaT")); // Number of media in the simulation.
    
    float *Temp           = mxGetData(prhs[0]); // Temp is an nx*ny*nz array of floats (singles in MATLAB)
    float *outputTemp     = mxGetData(plhs[0]); // outputTemp is an nx*ny*nz array of floats
    long nt               = *mxGetPr(mxGetField(prhs[2],0,"steps")); // nt is an integer (number of time steps to perform)
    float dt              = *mxGetPr(mxGetField(prhs[2],0,"dt")); // dt is a float (duration of time step)
    unsigned char *M      = mxGetData(mxGetField(prhs[2],0,"M")); // M is a nx*ny*nz array of uint8 (unsigned char) containing values from 0..nM-1
    float *c              = mxGetData(mxGetField(prhs[2],0,"dTdtperdeltaT")); // c is an nM*nM*3 array of floats
    float *dTdt_abs       = mxGetData(mxGetField(prhs[2],0,"dTdt_abs")); // dTdt_abs is an nx*ny*nz array of floats
    double *A             = mxGetPr(mxGetField(prhs[2],0,"A")); // A is a nM array of doubles
    double *E             = mxGetPr(mxGetField(prhs[2],0,"E")); // E is a nM array of doubles
    bool lightsOn         = mxIsLogicalScalarTrue(mxGetField(prhs[2],0,"lightsOn"));
    long heatSinked       = *mxGetPr(mxGetField(prhs[2],0,"heatBoundaryType")); // 0: Insulating boundaries, 1: Constant-temperature boundaries
    
    float *Omega        = mxGetData(prhs[1]); // Omega is an nx*ny*nz array of floats if we are supposed to calculate damage, a single NaN element otherwise
    bool calcDamage     = !mxIsNaN(Omega[0]); // If the Omega input is just one NaN element then we shouldn't bother with thermal damage calculation
    plhs[1] = calcDamage? mxCreateNumericArray(3,dimPtr,mxSINGLE_CLASS,mxREAL): mxCreateNumericArray(2,(mwSize[]){1, 1},mxSINGLE_CLASS,mxREAL); // outputOmega is the same dimensions as Omega
    float *outputOmega = mxGetData(plhs[1]);
    if(!calcDamage) outputOmega[0] = NAN;
    
    double *tempSensorCornerIdxs = mxGetPr(mxGetField(prhs[2],0,"tempSensorCornerIdxs"));
    long nSensors = mxGetM(mxGetField(prhs[2],0,"tempSensorCornerIdxs"));
    double *tempSensorInterpWeights = mxGetPr(mxGetField(prhs[2],0,"tempSensorInterpWeights"));
    
    plhs[2] = nSensors? mxCreateDoubleMatrix(nSensors,nt+1,mxREAL): mxCreateDoubleMatrix(0,0,mxREAL);
    double *sensorTemps = mxGetPr(plhs[2]);
    
    plhs[3] = mxCreateNumericArray(2,(mwSize[]){nM, 1},mxSINGLE_CLASS,mxREAL);
    float *maxMediaTemps = mxGetData(plhs[3]);
    
    float *tempTemp = malloc(nx*ny*nz*sizeof(float)); // Temporary temperature matrix
    
    float **srcTemp = &Temp; // Pointer to the pointer to whichever temperature matrix is to be read from
    float **tgtTemp; // Pointer to the pointer to whichever temperature matrix is to be written to
    tgtTemp = (nt%2)? &outputTemp: &tempTemp; // If nt is odd then we can start by writing to the output temperature matrix (pointed to by plhs[0]), otherwise we start by writing to the temporary temperature matrix
    
    for(long j=0; j<nSensors; j++) { // Interpolate to get the starting temperatures on the temperature sensors
        long idx = tempSensorCornerIdxs[j];
        double wx = tempSensorInterpWeights[j           ];
        double wy = tempSensorInterpWeights[j+  nSensors];
        double wz = tempSensorInterpWeights[j+2*nSensors];
        sensorTemps[j] = (1-wx)*(1-wy)*(1-wz)*Temp[idx     ] + (1-wx)*(1-wy)*wz*Temp[idx     +nx*ny] +
                         (1-wx)*   wy *(1-wz)*Temp[idx  +nx] + (1-wx)*   wy *wz*Temp[idx  +nx+nx*ny] +
                            wx *(1-wy)*(1-wz)*Temp[idx+1   ] +    wx *(1-wy)*wz*Temp[idx+1   +nx*ny] +
                            wx *   wy *(1-wz)*Temp[idx+1+nx] +    wx *   wy *wz*Temp[idx+1+nx+nx*ny];
    }
    
	#ifdef _WIN32
    bool useAllCPUs       = mxIsLogicalScalarTrue(mxGetField(prhs[2],0,"useAllCPUs"));
    #pragma omp parallel num_threads(useAllCPUs || omp_get_num_procs() == 1? omp_get_num_procs(): omp_get_num_procs()-1)
    #endif
    {
        long i,ix,iy,iz,ib,n;
        float dTdt,w;
		float b[nx>ny? nx: ny];
		float bz[nx*nz];
		if(heatSinked) {
			b[0] = 1;
			for(ix=0; ix<nx; ix++) bz[ix] = 1;
		}
		
		if(calcDamage) {
			// Initialize omega array
			#ifdef _WIN32
			#pragma omp for schedule(auto)
			#endif
			for(long idx=0; idx<nx*ny*nz; idx++) outputOmega[idx] = Omega[idx];
		}
		
        for(n=0; n<nt; n++) {
            if(ctrlc_caught) break;
			
            // In the following, the Douglas-Gunn Alternating Direction Implicit (DG-ADI) method is used to propagate the heat. This is based on equations (3.23a) and (3.23b) in me690-lctr-nts.pdf.
            
            // Explicit part of substep 1 out of 3
			#ifdef _WIN32
            #pragma omp for schedule(auto)
            #endif
            for(iz=0; iz<nz; iz++) {
                for(iy=0; iy<ny; iy++) {
                    for(ix=0; ix<nx; ix++) {
                        i = ix + iy*nx + iz*nx*ny;
						
                        if(!(heatSinked && (ix == 0 || ix == nx-1 || iy == 0 || iy == ny-1 || iz == 0 || iz == nz-1))) { // If not constant-temperature boundary voxel
							dTdt = lightsOn? dTdt_abs[i]: 0; // This is where the heat is deposited
                            if(ix != 0   ) dTdt += (T1[i-1]     - T1[i])*c[M[i] + M[i-1    ]*nM          ]/2; // Factor ½ is because the first substep splits the x heat flow over an explicit and an implicit part
                            if(ix != nx-1) dTdt += (T1[i+1]     - T1[i])*c[M[i] + M[i+1    ]*nM          ]/2; // Factor ½ is because the first substep splits the x heat flow over an explicit and an implicit part
                            if(iy != 0   ) dTdt += (T1[i-nx]    - T1[i])*c[M[i] + M[i-nx   ]*nM +   nM*nM];
                            if(iy != ny-1) dTdt += (T1[i+nx]    - T1[i])*c[M[i] + M[i+nx   ]*nM +   nM*nM];
                            if(iz != 0   ) dTdt += (T1[i-nx*ny] - T1[i])*c[M[i] + M[i-nx*ny]*nM + 2*nM*nM];
                            if(iz != nz-1) dTdt += (T1[i+nx*ny] - T1[i])*c[M[i] + M[i+nx*ny]*nM + 2*nM*nM];
                        } else dTdt = 0;
                        T2[i] = T1[i] + dt*dTdt;
                    }
                }
            }
			
			// Implicit part of substep 1 out of 3
			#ifdef _WIN32
            #pragma omp for schedule(auto)
            #endif
			for(iz=0; iz<nz; iz++) {
                for(iy=0; iy<ny; iy++) {
					// Thomson algorithm, sweeps up from 0 to nx-1 and then down from nx-1 to 0:
					// Algorithm is taken from https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
					for(ix=0; ix<nx; ix++) {
						if(!(heatSinked && (ix == 0 || ix == nx-1 || iy == 0 || iy == ny-1 || iz == 0 || iz == nz-1))) {
							b[ix] = 1;
							i = ix + iy*nx + iz*nx*ny;
							if(ix < nx-1) b[ix] += dt/2*c[M[i] + nM*M[i+1]];
							if(ix > 0) {
								b[ix] +=  dt/2*c[M[i] + nM*M[i-1]];
								w      = -dt/2*c[M[i] + nM*M[i-1]]/b[ix-1];
								if(ix > 1 || !heatSinked) b[ix] += w*dt/2*c[M[i-1] + nM*M[i]];
								T2[i] -= w*T2[i-1];
							}
						}
					}
					
					for(ix=nx-1; ix>=0; ix--) {
						if(!(heatSinked && (ix == 0 || ix == nx-1 || iy == 0 || iy == ny-1 || iz == 0 || iz == nz-1))) {
							i = ix + iy*nx + iz*nx*ny;
							T2[i] = (T2[i] + (ix == nx-1? 0: dt/2*c[M[i] + nM*M[i+1]]*T2[i+1]))/b[ix];
						}
					}
				}
			}
			
            // Explicit part of substep 2 out of 3
			#ifdef _WIN32
            #pragma omp for schedule(auto)
            #endif
            for(iz=0; iz<nz; iz++) {
                for(iy=0; iy<ny; iy++) {
                    for(ix=0; ix<nx; ix++) {
                        if(!(heatSinked && (ix == 0 || ix == nx-1 || iy == 0 || iy == ny-1 || iz == 0 || iz == nz-1))) { // If not constant-temperature boundary voxel
							i = ix + iy*nx + iz*nx*ny;
							dTdt = 0;
                            if(iy != 0   ) dTdt -= (T1[i-nx]    - T1[i])*c[M[i] + M[i-nx   ]*nM +   nM*nM]/2;
                            if(iy != ny-1) dTdt -= (T1[i+nx]    - T1[i])*c[M[i] + M[i+nx   ]*nM +   nM*nM]/2;
							T2[i] = T2[i] + dt*dTdt;
                        }
                    }
                }
            }
			
			// Implicit part of substep 2 out of 3
			#ifdef _WIN32
            #pragma omp for schedule(auto)
            #endif
			for(iz=0; iz<nz; iz++) {
                for(ix=0; ix<nx; ix++) {
					for(iy=0; iy<ny; iy++) {
						if(!(heatSinked && (ix == 0 || ix == nx-1 || iy == 0 || iy == ny-1 || iz == 0 || iz == nz-1))) {
							b[iy] = 1;
							i = ix + iy*nx + iz*nx*ny;
							if(iy < ny-1) b[iy] += dt/2*c[M[i] + nM*M[i+nx] + nM*nM];
							if(iy > 0) {
								b[iy] +=  dt/2*c[M[i] + nM*M[i-nx] + nM*nM];
								w      = -dt/2*c[M[i] + nM*M[i-nx] + nM*nM]/b[iy-1];
								if(iy > 1 || !heatSinked) b[iy] += w*dt/2*c[M[i-nx] + nM*M[i] + nM*nM];
								T2[i] -= w*T2[i-nx];
							}
						}
					}

					for(iy=ny-1; iy>=0; iy--) {
						if(!(heatSinked && (ix == 0 || ix == nx-1 || iy == 0 || iy == ny-1 || iz == 0 || iz == nz-1))) {
							i = ix + iy*nx + iz*nx*ny;
							T2[i] = (T2[i] + (iy == ny-1? 0: dt/2*c[M[i] + nM*M[i+nx] + nM*nM]*T2[i+nx]))/b[iy];
						}
					}
				}
			}
			
			// Explicit part of substep 3 out of 3
			#ifdef _WIN32
            #pragma omp for schedule(auto)
            #endif
            for(iz=0; iz<nz; iz++) {
                for(iy=0; iy<ny; iy++) {
                    for(ix=0; ix<nx; ix++) {
                        if(!(heatSinked && (ix == 0 || ix == nx-1 || iy == 0 || iy == ny-1 || iz == 0 || iz == nz-1))) { // If not constant-temperature boundary voxel
							i = ix + iy*nx + iz*nx*ny;
							dTdt = 0;
                            if(iz != 0   ) dTdt -= (T1[i-nx*ny] - T1[i])*c[M[i] + M[i-nx*ny]*nM + 2*nM*nM]/2;
                            if(iz != nz-1) dTdt -= (T1[i+nx*ny] - T1[i])*c[M[i] + M[i+nx*ny]*nM + 2*nM*nM]/2;
							T2[i] = T2[i] + dt*dTdt;
                        }
                    }
                }
            }
			
			// Implicit part of substep 3 out of 3
			#ifdef _WIN32
            #pragma omp for schedule(auto)
            #endif
			for(iy=0; iy<ny; iy++) {
				for(iz=0; iz<nz; iz++) {
					for(ix=0; ix<nx; ix++) { // The x sweep is put inside the z sweep on purpose, since it improves performance a lot in this step. The b coefficients used here for the Thomson algorithm are stored in a 2D (x,z) matrix called bz.
						if(!(heatSinked && (ix == 0 || ix == nx-1 || iy == 0 || iy == ny-1 || iz == 0 || iz == nz-1))) {
							ib = ix + iz*nx;
							bz[ib] = 1;
							i = ix + iy*nx + iz*nx*ny;
							if(iz < nz-1) bz[ib] += dt/2*c[M[i] + nM*M[i+nx*ny] + 2*nM*nM];
							if(iz > 0) {
								bz[ib] +=  dt/2*c[M[i] + nM*M[i-nx*ny] + 2*nM*nM];
								w       = -dt/2*c[M[i] + nM*M[i-nx*ny] + 2*nM*nM]/bz[ib-nx];
								if(iz > 1 || !heatSinked) bz[ib] += w*dt/2*c[M[i-nx*ny] + nM*M[i] + 2*nM*nM];
								T2[i] -= w*T2[i-nx*ny];
							}
						}
					}
				}
				for(iz=nz-1; iz>=0; iz--) {
					for(ix=0; ix<nx; ix++) {
						i = ix + iy*nx + iz*nx*ny;
						if(!(heatSinked && (ix == 0 || ix == nx-1 || iy == 0 || iy == ny-1 || iz == 0 || iz == nz-1))) {
							ib = ix + iz*nx;
							T2[i] = (T2[i] + (iz == nz-1? 0: dt/2*c[M[i] + nM*M[i+nx*ny] + 2*nM*nM]*T2[i+nx*ny]))/bz[ib];
						}
						if(calcDamage && E[M[i]]) { // Arrhenius damage integral evaluation
							if(T2[i] == T1[i]) outputOmega[i] += (float)(dt*A[M[i]]*exp(-E[M[i]]/(R*(T2[i] + CELSIUSZERO)))); // For use when the temperature in a voxel (extremely rarely) is constant to float precision
                            else { // Assume a linear temperature slope from the previous to the current time step and use the mathematical integral, see for example https://www.wolframalpha.com/input/?i=int(A*exp(1%2F(a*t%2Bb)),t)
								float a = -R/E[M[i]]*(T2[i] - T1[i])/dt; // "a" parameter in the Wolfram Alpha evaluation above
								float b = -R/E[M[i]]*(T1[i] + CELSIUSZERO); // "b" parameter in the Wolfram Alpha evaluation above
								outputOmega[i] += (float)(A[M[i]]/a*(exp(1/(a*dt+b))*(a*dt+b+expint_E1_scaled(-1/(a*dt+b))) -
																	 exp(1/      b )*(     b+expint_E1_scaled(-1/      b)))); // The expint_E1_scaled(x) function is equal to -Ei(-x)*exp(x)
							}
						}
                        if(maxMediaTemps[M[i]] < T2[i]) maxMediaTemps[M[i]] = T2[i]; // If the temperature in this voxel is a new record high for this medium, write the value into the maxMediaTemps array
					}
				}
			}
			
			#ifdef _WIN32
			#pragma omp master
			#endif
			{
				#ifdef _WIN32
				if(utIsInterruptPending()) {
					ctrlc_caught = true;
					printf("\nCtrl+C detected, stopping.\n");
				}
				#endif

				for(long j=0; j<nSensors; j++) { // Interpolate to get the new temperatures on the temperature sensors
					long idx = tempSensorCornerIdxs[j];
					double wx = tempSensorInterpWeights[j           ];
					double wy = tempSensorInterpWeights[j+  nSensors];
					double wz = tempSensorInterpWeights[j+2*nSensors];
					sensorTemps[j+(n+1)*nSensors] = (1-wx)*(1-wy)*(1-wz)*T2[idx     ] + (1-wx)*(1-wy)*wz*T2[idx     +nx*ny] +
													(1-wx)*   wy *(1-wz)*T2[idx  +nx] + (1-wx)*   wy *wz*T2[idx  +nx+nx*ny] +
													   wx *(1-wy)*(1-wz)*T2[idx+1   ] +    wx *(1-wy)*wz*T2[idx+1   +nx*ny] +
													   wx *   wy *(1-wz)*T2[idx+1+nx] +    wx *   wy *wz*T2[idx+1+nx+nx*ny];
				}

				if(tgtTemp == &outputTemp) {
					tgtTemp = &tempTemp;
					srcTemp = &outputTemp;
				} else {
					tgtTemp = &outputTemp;
					srcTemp = &tempTemp;
				}
			}
			#ifdef _WIN32
			#pragma omp barrier
			#endif
		}
    }
    free(tempTemp);
    return;
}
