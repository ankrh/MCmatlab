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
#ifdef _WIN32 // This is defined on both win32 and win64 systems. We use this preprocessor condition to avoid loading openmp or libut on, e.g., Mac
#include "omp.h"
extern bool utIsInterruptPending(); // Allows catching ctrl+c while executing the mex function
#endif

#define R 8.3144598 // Gas constant [J/mol/K]
#define CELSIUSZERO 273.15
#define T_xyz (firstStep? T_in: T_xyz_notfirststep)
// #define BLOCKSIZE 256
// 
// void transpose(double *A, double *B, long N, long M) {
// 	// Transposes the MxN matrix A and stores result in B
// 	#ifdef _WIN32
// 	#pragma omp for schedule(auto)
// 	#endif
// 	for(long i=0; i<N; i+=BLOCKSIZE) {
// 		long max_i2 = i+BLOCKSIZE < N ? i + BLOCKSIZE : N;
// 		for(long j=0; j<M; j+=BLOCKSIZE) {
// 			long max_j2 = j+BLOCKSIZE < M ? j + BLOCKSIZE : M;
// 			for(long i2=i; i2<max_i2; i2++) {
// 				for(long j2=j; j2<max_j2; j2++) {
// 					B[i2 + j2*N] = A[j2 + i2*M];
// 				}
// 			}
// 		}
// 	}
// }

void mexFunction( int nlhs, mxArray *plhs[],int nrhs, mxArray const *prhs[] ) {
    bool ctrlc_caught = false;           // Has a ctrl+c been passed from MATLAB?
	bool firstStep = true;
    
    mwSize const *dimPtr = mxGetDimensions(prhs[0]);
    long nx = dimPtr[0];
    long ny = dimPtr[1];
    long nz = dimPtr[2];
    long nM = mxGetM(mxGetField(prhs[2],0,"dTdtperdeltaT")); // Number of media in the simulation.
    
    
    double *T_in          = mxGetPr(prhs[0]); // T_in is an nx*ny*nz array of doubles
    plhs[0]               = mxCreateNumericArray(3,dimPtr,mxDOUBLE_CLASS,mxREAL);
    long nt               = *mxGetPr(mxGetField(prhs[2],0,"steps")); // nt is an integer (number of time steps to perform)
    double dt             = *mxGetPr(mxGetField(prhs[2],0,"dt")); // dt is a double (duration of time step)
    unsigned char *M      = mxGetData(mxGetField(prhs[2],0,"M")); // M is a nx*ny*nz array of uint8 (unsigned char) containing values from 0..nM-1
    double *c             = mxGetPr(mxGetField(prhs[2],0,"dTdtperdeltaT")); // c is an nM*nM*3 array of doubles
    double *dTdt_abs      = mxGetPr(mxGetField(prhs[2],0,"dTdt_abs")); // dTdt_abs is an nx*ny*nz array of doubles
    double *A             = mxGetPr(mxGetField(prhs[2],0,"A")); // A is a nM array of doubles
    double *E             = mxGetPr(mxGetField(prhs[2],0,"E")); // E is a nM array of doubles
    bool lightsOn         = mxIsLogicalScalarTrue(mxGetField(prhs[2],0,"lightsOn"));
    long heatSinked       = *mxGetPr(mxGetField(prhs[2],0,"heatBoundaryType")); // 0: Insulating boundaries, 1: Constant-temperature boundaries
    
    double *T_xyz_notfirststep = mxGetPr(plhs[0]); // T_xyz_notfirststep is a nx*ny*nz array of doubles, to be used both in each iteration and for final output. In most of the code it will be referred to as T_xyz, using a #define
	double *T_yzx = malloc(nx*ny*nz*sizeof(double));
	double *T_zxy = malloc(nx*ny*nz*sizeof(double));

	double *Omega       = mxGetPr(prhs[1]); // Omega is an nx*ny*nz array of doubles if we are supposed to calculate damage, a single NaN element otherwise
    bool calcDamage     = !mxIsNaN(Omega[0]); // If the Omega input is just one NaN element then we shouldn't bother with thermal damage calculation
    plhs[1] = calcDamage? mxCreateNumericArray(3,dimPtr,mxDOUBLE_CLASS,mxREAL): mxCreateDoubleMatrix(1,1,mxREAL); // outputOmega is the same dimensions as Omega
    double *outputOmega = mxGetPr(plhs[1]);
    if(!calcDamage) outputOmega[0] = NAN;
    
    double *tempSensorCornerIdxs = mxGetPr(mxGetField(prhs[2],0,"tempSensorCornerIdxs"));
    long nSensors = mxGetM(mxGetField(prhs[2],0,"tempSensorCornerIdxs"));
    double *tempSensorInterpWeights = mxGetPr(mxGetField(prhs[2],0,"tempSensorInterpWeights"));
    
    plhs[2] = nSensors? mxCreateDoubleMatrix(nSensors,nt+1,mxREAL): mxCreateDoubleMatrix(0,0,mxREAL);
    double *sensorTemps = mxGetPr(plhs[2]);
	
    for(long j=0;j<nSensors;j++) { // Interpolate to get the starting temperatures on the temperature sensors
        long idx = tempSensorCornerIdxs[j];
        double wx = tempSensorInterpWeights[j           ];
        double wy = tempSensorInterpWeights[j+  nSensors];
        double wz = tempSensorInterpWeights[j+2*nSensors];
        sensorTemps[j] = (1-wx)*(1-wy)*(1-wz)*T_in[idx     ] + (1-wx)*(1-wy)*wz*T_in[idx     +nx*ny] +
                         (1-wx)*   wy *(1-wz)*T_in[idx  +nx] + (1-wx)*   wy *wz*T_in[idx  +nx+nx*ny] +
                            wx *(1-wy)*(1-wz)*T_in[idx+1   ] +    wx *(1-wy)*wz*T_in[idx+1   +nx*ny] +
                            wx *   wy *(1-wz)*T_in[idx+1+nx] +    wx *   wy *wz*T_in[idx+1+nx+nx*ny];
    }
    
    #ifdef _WIN32
    bool useAllCPUs       = mxIsLogicalScalarTrue(mxGetField(prhs[2],0,"useAllCPUs"));
    #pragma omp parallel num_threads(useAllCPUs || omp_get_num_procs() == 1? omp_get_num_procs(): omp_get_num_procs()-1)
    #endif
    {
		double w,b[nx>ny && nx>nz? nx: (ny>nz? ny: nz)]; // Length of b array is maximum of nx, ny and nz
		double *dT_y_array = malloc(nx*sizeof(double));
		double *T_xline = malloc(nx*sizeof(double));
        
		if(calcDamage) {
			// Initialize omega array
			#ifdef _WIN32
			#pragma omp for schedule(auto)
			#endif
			for(long i_xyz=0;i_xyz<nx*ny*nz;i_xyz++) outputOmega[i_xyz] = Omega[i_xyz];
		}
		
        for(long n=0; n<nt; n++) {
            if(ctrlc_caught) break;
			
            // Substep 1 out of 3, includes explicit x, implicit x and explicit y
			#ifdef _WIN32
            #pragma omp for schedule(auto)
            #endif
            for(long iz=0; iz<nz; iz++) {
                for(long iy=0; iy<ny; iy++) {
					if(!(heatSinked && (iy == 0 || iy == ny-1 || iz == 0 || iz == nz-1))) {
						for(long ix=0; ix<nx; ix++) { // Sweep up in x
							long i_xyz = ix + iy*nx + iz*nx*ny; // xyz index
							long i_zxy = iz + ix*nz + iy*nz*nx; // zxy index

							// Explicit x
							double dTdt_x = 0, dTdt_y = 0, dTdt_z = 0;
							bool heatSinkedVoxel = heatSinked && (ix == 0 || ix == nx-1);
							if(!heatSinkedVoxel) { // If not constant-temperature boundary voxel
								if(ix != 0   ) dTdt_x += (T_xyz[i_xyz-1]     - T_xyz[i_xyz])*c[M[i_xyz] + M[i_xyz-1    ]*nM          ];
								if(ix != nx-1) dTdt_x += (T_xyz[i_xyz+1]     - T_xyz[i_xyz])*c[M[i_xyz] + M[i_xyz+1    ]*nM          ];
								if(iy != 0   ) dTdt_y += (T_xyz[i_xyz-nx]    - T_xyz[i_xyz])*c[M[i_xyz] + M[i_xyz-nx   ]*nM +   nM*nM];
								if(iy != ny-1) dTdt_y += (T_xyz[i_xyz+nx]    - T_xyz[i_xyz])*c[M[i_xyz] + M[i_xyz+nx   ]*nM +   nM*nM];
								if(iz != 0   ) dTdt_z += (T_xyz[i_xyz-nx*ny] - T_xyz[i_xyz])*c[M[i_xyz] + M[i_xyz-nx*ny]*nM + 2*nM*nM];
								if(iz != nz-1) dTdt_z += (T_xyz[i_xyz+nx*ny] - T_xyz[i_xyz])*c[M[i_xyz] + M[i_xyz+nx*ny]*nM + 2*nM*nM];
							}
							dT_y_array[ix] = dt*dTdt_y/2;
							T_zxy[i_zxy] = -dt*dTdt_z/2; // Preparation for the z explicit step

							T_xline[ix] = T_xyz[i_xyz] + dt*((lightsOn && !heatSinkedVoxel? dTdt_abs[i_xyz]: 0) + dTdt_x/2 + dTdt_y + dTdt_z);

							// Implicit x up sweep
							// Thomson algorithm, sweeps up from 0 to nx-1 and then down from nx-1 to 0:
							// Algorithm is taken from https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
							if(ix==0) {
								b[0] = heatSinked? 1: 1 + dt/2*c[M[i_xyz] + nM*M[i_xyz+1]];
							} else if(ix < nx-1) {
								w = -dt/2*c[M[i_xyz] + nM*M[i_xyz-1]]/b[ix-1];
								b[ix] = 1 + dt/2*c[M[i_xyz] + nM*M[i_xyz-1]] + dt/2*c[M[i_xyz] + nM*M[i_xyz+1]];
								if(ix > 1 || !heatSinked) b[ix] += w*dt/2*c[M[i_xyz-1] + nM*M[i_xyz]];
								T_xline[ix] -= w*T_xline[ix-1];
							} else { // ix == nx-1
								w = heatSinked? 0: -dt/2*c[M[i_xyz] + nM*M[i_xyz-1]]/b[ix-1];
								b[ix] = heatSinked? 1: 1 + dt/2*c[M[i_xyz] + nM*M[i_xyz-1]];
								b[ix] += w*dt/2*c[M[i_xyz-1] + nM*M[i_xyz]];
								T_xline[ix] -= w*T_xline[ix-1];
							}
						}
						
						double T_temp = 0;
						for(long ix=nx-1;ix>=0;ix--) { // Sweep down in x, implicit x and explicit y, store result in T_yzx
							long i_xyz = ix + iy*nx + iz*nx*ny; // xyz index
							long i_yzx = iy + iz*ny + ix*ny*nz; // yzx index
							if(ix==nx-1) T_temp = T_xline[ix]/b[ix];
							else if(ix > 0 || !heatSinked) T_temp = (T_xline[ix] + dt/2*c[M[i_xyz] + nM*M[i_xyz+1]]*T_temp)/b[ix]; // The T_temp on the right hand side is the one from the previous iteration
							else T_temp = T_xline[ix];
							T_yzx[i_yzx] = T_temp - dT_y_array[ix];
						}
					} else for(long ix=nx-1;ix>=0;ix--) { // Heatsinked x line
						long i_xyz = ix + iy*nx + iz*nx*ny; // xyz index
						long i_yzx = iy + iz*ny + ix*ny*nz; // yzx index
						long i_zxy = iz + ix*nz + iy*nz*nx; // zxy index
						T_zxy[i_zxy] = 0; // Prepare the zxy matrix with an explicit component (dTdt_z) of 0
						T_yzx[i_yzx] = T_xyz[i_xyz];
					}
                }
            }

            // Substep 2 out of 3, includes implicit y and explicit z
			#ifdef _WIN32
            #pragma omp for schedule(auto)
            #endif
            for(long iz=0; iz<nz; iz++) {
				for(long ix=0; ix<nx; ix++) {
					if(!(heatSinked && (ix == 0 || ix == nx-1 || iz == 0 || iz == nz-1))) {
						for(long iy=0; iy<ny; iy++) { // Sweep up in y
							long i_xyz = ix + iy*nx + iz*nx*ny; // xyz index
							long i_yzx = iy + iz*ny + ix*ny*nz; // yzx index
							if(iy==0) {
								b[0] = heatSinked? 1: 1 + dt/2*c[nM*nM + M[i_xyz] + nM*M[i_xyz+nx]];
							} else if(iy < ny-1) {
								w = -dt/2*c[nM*nM + M[i_xyz] + nM*M[i_xyz-nx]]/b[iy-1];
								b[iy] = 1 + dt/2*c[nM*nM + M[i_xyz] + nM*M[i_xyz-nx]] + dt/2*c[nM*nM + M[i_xyz] + nM*M[i_xyz+nx]];
								if(iy > 1 || !heatSinked) b[iy] += w*dt/2*c[nM*nM + M[i_xyz-nx] + nM*M[i_xyz]];
								T_yzx[i_yzx] -= w*T_yzx[i_yzx-1];
							} else { // iy == ny-1
								w = heatSinked? 0: -dt/2*c[nM*nM + M[i_xyz] + nM*M[i_xyz-nx]]/b[iy-1];
								b[iy] = heatSinked? 1: 1 + dt/2*c[nM*nM + M[i_xyz] + nM*M[i_xyz-nx]];
								b[iy] += w*dt/2*c[nM*nM + M[i_xyz-nx] + nM*M[i_xyz]];
								T_yzx[i_yzx] -= w*T_yzx[i_yzx-1];
							}
						}

						double T_temp = 0;
						for(long iy=ny-1;iy>=0;iy--) { // Sweep down in y, store result in T_zxy, which already contains -dTdt_z/2, the z explicit correction
							long i_xyz = ix + iy*nx + iz*nx*ny; // xyz index
							long i_yzx = iy + iz*ny + ix*ny*nz; // yzx index
							long i_zxy = iz + ix*nz + iy*nz*nx; // zxy index
							
							if(iy==ny-1) T_temp = T_yzx[i_yzx]/b[iy];
							else if(iy > 0 || !heatSinked) T_temp = (T_yzx[i_yzx] + dt/2*c[nM*nM + M[i_xyz] + nM*M[i_xyz+nx]]*T_temp)/b[iy];
							else T_temp = T_yzx[i_yzx];
							T_zxy[i_zxy] += T_temp;
						}
					} else for(long iy=ny-1;iy>=0;iy--) { // Heatsinked y line
						long i_yzx = iy + iz*ny + ix*ny*nz; // yzx index
						long i_zxy = iz + ix*nz + iy*nz*nx; // zxy index
						T_zxy[i_zxy] += T_yzx[i_yzx];
					}
                }
            }
			
            // Substep 3 out of 3, includes implicit z
			#ifdef _WIN32
            #pragma omp for schedule(auto)
            #endif
            for(long iy=0; iy<ny; iy++) {
				for(long ix=0; ix<nx; ix++) {
					if(!(heatSinked && (ix == 0 || ix == nx-1 || iy == 0 || iy == ny-1))) {
						for(long iz=0; iz<nz; iz++) { // Sweep up in z
							long i_xyz = ix + iy*nx + iz*nx*ny; // xyz index
							long i_zxy = iz + ix*nz + iy*nz*nx; // zxy index
							if(iz==0) {
								b[0] = heatSinked? 1: 1 + dt/2*c[2*nM*nM + M[i_xyz] + nM*M[i_xyz+nx*ny]];
							} else if(iz < nz-1) {
								w = -dt/2*c[2*nM*nM + M[i_xyz] + nM*M[i_xyz-nx*ny]]/b[iz-1];
								b[iz] = 1 + dt/2*c[2*nM*nM + M[i_xyz] + nM*M[i_xyz-nx*ny]] + dt/2*c[2*nM*nM + M[i_xyz] + nM*M[i_xyz+nx*ny]];
								if(iz > 1 || !heatSinked) b[iz] += w*dt/2*c[2*nM*nM + M[i_xyz-nx*ny] + nM*M[i_xyz]];
								T_zxy[i_zxy] -= w*T_zxy[i_zxy-1];
							} else { // iz == nz-1
								w = heatSinked? 0: -dt/2*c[2*nM*nM + M[i_xyz] + nM*M[i_xyz-nx*ny]]/b[iz-1];
								b[iz] = heatSinked? 1: 1 + dt/2*c[2*nM*nM + M[i_xyz] + nM*M[i_xyz-nx*ny]];
								b[iz] += w*dt/2*c[2*nM*nM + M[i_xyz-nx*ny] + nM*M[i_xyz]];
								T_zxy[i_zxy] -= w*T_zxy[i_zxy-1];
							}
						}

						for(long iz=nz-1;iz>=0;iz--) { // Sweep down in z, store result in T_xyz
							long i_xyz = ix + iy*nx + iz*nx*ny; // xyz index
							long i_zxy = iz + ix*nz + iy*nz*nx; // zxy index
							
							double T_new;
							if(iz==nz-1) T_new = T_zxy[i_zxy]/b[iz];
							else if(iz > 0 || !heatSinked) T_new = (T_zxy[i_zxy] + dt/2*c[2*nM*nM + M[i_xyz] + nM*M[i_xyz+nx*ny]]*T_xyz_notfirststep[i_xyz+nx*ny])/b[iz];
							else T_new = T_zxy[i_zxy];

							if(calcDamage) outputOmega[i_xyz] += dt*A[M[i_xyz]]*exp(-E[M[i_xyz]]/(R*((T_xyz[i_xyz] + T_new)/2 + CELSIUSZERO))); // Arrhenius damage integral evaluation
							T_xyz_notfirststep[i_xyz] = T_new; // In the first step, T_xyz will be the matlab input, which we are not allowed to write to. Therefore we use T_xyz_notfirststep here. See the #define for T_xyz.
						}
					} else for(long iz=nz-1;iz>=0;iz--) { // Heatsinked z line
						long i_xyz = ix + iy*nx + iz*nx*ny; // xyz index
						long i_zxy = iz + ix*nz + iy*nz*nx; // zxy index
						if(calcDamage) outputOmega[i_xyz] += dt*A[M[i_xyz]]*exp(-E[M[i_xyz]]/(R*(T_xyz[i_xyz] + CELSIUSZERO))); // Arrhenius damage integral evaluation
						T_xyz_notfirststep[i_xyz] = T_zxy[i_zxy]; // In the first step, T_xyz will be the matlab input, which we are not allowed to write to. Therefore we use T_xyz_notfirststep here. See the #define for T_xyz.
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

				firstStep = false;

				for(long j=0;j<nSensors;j++) { // Interpolate to get the new temperatures on the temperature sensors
					long idx = tempSensorCornerIdxs[j];
					double wx = tempSensorInterpWeights[j           ];
					double wy = tempSensorInterpWeights[j+  nSensors];
					double wz = tempSensorInterpWeights[j+2*nSensors];
					sensorTemps[j+(n+1)*nSensors] = (1-wx)*(1-wy)*(1-wz)*T_xyz[idx     ] + (1-wx)*(1-wy)*wz*T_xyz[idx     +nx*ny] +
													(1-wx)*   wy *(1-wz)*T_xyz[idx  +nx] + (1-wx)*   wy *wz*T_xyz[idx  +nx+nx*ny] +
													   wx *(1-wy)*(1-wz)*T_xyz[idx+1   ] +    wx *(1-wy)*wz*T_xyz[idx+1   +nx*ny] +
													   wx *   wy *(1-wz)*T_xyz[idx+1+nx] +    wx *   wy *wz*T_xyz[idx+1+nx+nx*ny];
				}
			}
			#ifdef _WIN32
			#pragma omp barrier
			#endif
		}
		free(dT_y_array);
		free(T_xline);
    }
	free(T_yzx);
	free(T_zxy);
    return;
}
