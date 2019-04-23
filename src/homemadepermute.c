/********************************************
 * Permute test code
 * mex COPTIMFLAGS='$COPTIMFLAGS -msse -Ofast -fopenmp -std=c11 -Wall -pedantic' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall -pedantic' -outdir .\ .\src\homemadepermute.c ".\src\libut.lib"
******/

// Code for fast matrix transpose written by "user2088790", posted on https://stackoverflow.com/questions/16737298/what-is-the-fastest-way-to-transpose-a-matrix-in-c


#include <math.h>
#include "mex.h"
#include <x86intrin.h>
#ifdef _WIN32 // This is defined on both win32 and win64 systems. We use this preprocessor condition to avoid loading openmp or libut on, e.g., Mac
#include "omp.h"
// extern bool utIsInterruptPending(); // Allows catching ctrl+c while executing the mex function
#endif

// inline void transpose4x4_SSE(float *A, float *B, const int M, const int N) {
// 	__m128 row1 = _mm_load_ps(&A[0*M]);
// 	__m128 row2 = _mm_load_ps(&A[1*M]);
// 	__m128 row3 = _mm_load_ps(&A[2*M]);
// 	__m128 row4 = _mm_load_ps(&A[3*M]);
// 	_MM_TRANSPOSE4_PS(row1, row2, row3, row4);
// 	_mm_store_ps(&B[0*N], row1);
// 	_mm_store_ps(&B[1*N], row2);
// 	_mm_store_ps(&B[2*N], row3);
// 	_mm_store_ps(&B[3*N], row4);
// }

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[] ) {
//     bool ctrlc_caught = false;           // Has a ctrl+c been passed from MATLAB?
	
    float *A = mxGetData(prhs[0]); // Input matrix of floats. MATLAB seems to ensure that all matrix data is 32-byte aligned. SSE requires that data is 16-byte aligned (which is of course also satisfied).
    long M = mxGetM(prhs[0]);
    long N = mxGetN(prhs[0]);

    plhs[0] = mxCreateNumericArray(2,(mwSize[]){N, M},mxSINGLE_CLASS,mxREAL); // Like the input data, the allocated memory for output data is (apparently) also 32-byte aligned.
	float *B = mxGetData(plhs[0]);
	
	long block_size = *mxGetPr(prhs[1]);
	bool useSSE = mxIsLogicalScalarTrue(prhs[2]);
	
	#pragma omp parallel for
	for(long i=0; i<N; i+=block_size) {
		long max_i2 = i+block_size < N ? i + block_size : N;
		for(long j=0; j<M; j+=block_size) {
			long max_j2 = j+block_size < M ? j + block_size : M;
			if(useSSE) {
				for(long i2=i; i2<max_i2; i2+=4) {
					for(long j2=j; j2<max_j2; j2+=4) {
	// 					Transpose a 4x4 block using SSE intrinsics
						__m128 row1 = _mm_load_ps(A +  i2   *M + j2);
						__m128 row2 = _mm_load_ps(A + (i2+1)*M + j2);
						__m128 row3 = _mm_load_ps(A + (i2+2)*M + j2);
						__m128 row4 = _mm_load_ps(A + (i2+3)*M + j2);
						_MM_TRANSPOSE4_PS(row1, row2, row3, row4);
						_mm_store_ps(B +  j2   *N + i2, row1);
						_mm_store_ps(B + (j2+1)*N + i2, row2);
						_mm_store_ps(B + (j2+2)*N + i2, row3);
						_mm_store_ps(B + (j2+3)*N + i2, row4);
					}
				}
			} else {
				for(long i2=i; i2<max_i2; i2++) {
					for(long j2=j; j2<max_j2; j2++) {
						B[i2 + j2*N] = A[j2 + i2*M];
					}
				}
			}
		}
	}
	return;
}

