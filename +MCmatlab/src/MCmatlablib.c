
struct debug {
  double             dbls[3];
  unsigned long long ulls[3];
};

struct depositionCriteria { // Struct type for the deposition criteria, that determine whether a photon's weight loss should be deposited in any output arrays
  unsigned long  minS;
  unsigned long  maxS;
  unsigned long  minRefr;
  unsigned long  maxRefr;
  unsigned long  minRefl;
  unsigned long  maxRefl;
  unsigned long  minI;
  unsigned long  maxI;
  unsigned long  minIdx;
  unsigned long  maxIdx;
  bool           onlyCollected;
  bool           evaluateCriteriaAtEndOfLife;
};

struct geometry { // Struct type for the constant geometry definitions, including the wavelength-dependent optical properties and the boundary type
  FLOATORDBL     d[3];
  long           n[3];
  long           farFieldRes;
  int            boundaryType;
  FLOATORDBL     *muav,*musv,*gv,*RIv;
  unsigned char  *CDFidxv;
  FLOATORDBL     *CDFs;
  unsigned char  *M;
  float          *interfaceNormals;
};

struct source { // Struct type for the constant beam definitions
  int            beamType;
  FLOATORDBL     emitterLength;
  FLOATORDBL     *FPIDdist1; // Radial or X
  long           L_FPID1;
  FLOATORDBL     FPIDwidth1;
  FLOATORDBL     *FPIDdist2; // Azimuthal or Y
  long           L_FPID2;
  FLOATORDBL     FPIDwidth2;
  FLOATORDBL     *AIDdist1; // Radial or X
  long           L_AID1;
  FLOATORDBL     AIDwidth1;
  FLOATORDBL     *AIDdist2; // Azimuthal or Y
  long           L_AID2;
  FLOATORDBL     AIDwidth2;
  FLOATORDBL     *S;
  FLOATORDBL     power;
  FLOATORDBL     focus[3];
  FLOATORDBL     u[3];
  FLOATORDBL     v[3];
  FLOATORDBL     w[3];
};

struct lightCollector { // Struct type for the constant light collector definitions. It can be either an objective lens (for 0<f<INFINITY) or a fiber tip or simple aperture (for f=INFINITY)
  FLOATORDBL     r[3]; // Position of center of objective focal plane (not the objective itself) or position of center of the fiber tip
  FLOATORDBL     theta; // Polar angle that the objective or fiber is facing
  FLOATORDBL     phi; // Azimuthal angle of objective or fiber orientation
  FLOATORDBL     f; // Focal length of the objective. If the light collector is a fiber tip, this will be INFINITY.
  FLOATORDBL     diam; // Diameter of the objective aperture or core diameter of the fiber. For an ideal thin lens objective, this is 2*tan(arcsin(lensNA/f)).
  FLOATORDBL     FSorNA; // For an objective lens: Field Size of the imaging system (diameter of area in object plane that gets imaged). For a fiber tip: The fiber's NA.
  long           res[2]; // Resolution of image plane in pixels along the spatial and time axes. For a fiber, spatial resolution is 1.
  FLOATORDBL     tStart; // Start time for the interval used for binned time-resolved detection
  FLOATORDBL     tEnd; // End time for the interval used for binned time-resolved detection
};

struct photon { // Struct type for parameters describing the thread-specific current state of a photon
  FLOATORDBL     i[3],u[3],D[3]; // Fractional position indices i, ray trajectory unit vector u and distances D to next voxel boundary (yz, xz or xy) along current trajectory
  long           j; // Linear index of current voxel (or closest defined voxel if photon outside cuboid)
  FLOATORDBL     mua,mus,g,RI; // Absorption, scattering, anisotropy, refractive index values and a pointer to the phase function CDF array at current photon position
  unsigned char  CDFidx;
  FLOATORDBL     stepLeft,weight,time;
  bool           insideVolume,alive,sameVoxel;
  PRNG_t         PRNGstate; // "State" of the Mersenne Twister pseudo-random number generator
  long           recordSize; // Current size of the list of voxels in which power has been deposited, used only if depositionCriteria.evaluateCriteriaAtEndOfLife is true
  long           recordElems; // Current number of elements used of the record. Starts at 0 every photon launch, used only if depositionCriteria.evaluateCriteriaAtEndOfLife is true
  long           *j_record; // List of the indices of the voxels in which the current photon has deposited power, used only if depositionCriteria.evaluateCriteriaAtEndOfLife is true
  FLOATORDBL     *weight_record; // List of the weights that have been deposited into the voxels, used only if depositionCriteria.evaluateCriteriaAtEndOfLife is true
  unsigned long  scatterings;
  unsigned long  refractions;
  unsigned long  reflections;
  unsigned long  interfaceTransitions;
  char           killed_escaped_collected;
};

struct paths { // Struct type for storing the paths taken by the nExamplePaths first photons simulated by the master thread
  long           pathsElems; // Current number of elements used
  long           nExamplePaths; // Number of photons to store the path of
  long           nExamplePhotonPathsFinished; // Number of example photon paths the master thread has started
  bool           pathStartedThisPhoton;
  long           pathsSize; // Current size of the list
  FLOATORDBL     *data; // Array containing x, y, z and weight data for the photons, in which the paths of different photons are separated by four NaNs
};

struct MATLABoutputs {
  float * NFR;
  float * image;
  float * FF;
  float * NI_xpos;
  float * NI_xneg;
  float * NI_ypos;
  float * NI_yneg;
  float * NI_zpos;
  float * NI_zneg;
};

struct outputs {
  unsigned long long nPhotons;
  unsigned long long nPhotonsCollected;
  FLOATORDBL * NFR;
  FLOATORDBL * image;
  FLOATORDBL * FF;
  FLOATORDBL * NI_xpos;
  FLOATORDBL * NI_xneg;
  FLOATORDBL * NI_ypos;
  FLOATORDBL * NI_yneg;
  FLOATORDBL * NI_zpos;
  FLOATORDBL * NI_zneg;
};

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
bool depositionCriteriaMet(struct photon *P,struct depositionCriteria *DC) {
  return P->scatterings                                                                      >= DC->minS &&
         (!P->alive && P->killed_escaped_collected == 0? ULONG_MAX: P->scatterings)          <= DC->maxS &&
         P->refractions                                                                      >= DC->minRefr &&
         (!P->alive && P->killed_escaped_collected == 0? ULONG_MAX: P->refractions)          <= DC->maxRefr &&
         P->reflections                                                                      >= DC->minRefl &&
         (!P->alive && P->killed_escaped_collected == 0? ULONG_MAX: P->reflections)          <= DC->maxRefl &&
         P->interfaceTransitions                                                             >= DC->minI &&
         (!P->alive && P->killed_escaped_collected == 0? ULONG_MAX: P->interfaceTransitions) <= DC->maxI &&
         (DC->onlyCollected? P->killed_escaped_collected == 2: true);
}

unsigned long infCast(double x) {return mxIsInf(x)? ULONG_MAX: (unsigned long)x;}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
double sqr(double x) {return x*x;}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
float sqrf(float x) {return x*x;}

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600 && !defined(USEFLOATSFORGPU)
__device__ double atomicAdd(double* address, double val) {
  unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;

  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,__double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}
#endif

#ifdef __NVCC__ // If compiling for CUDA
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line) {
   if (code != cudaSuccess) {
      printf("GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      mexEvalString("drawnow;");
      while(true) {;}
   }
}
#endif

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void atomicAddWrapperULL(unsigned long long * ptr, unsigned long long val) {
  #ifdef __NVCC__ // If compiling for CUDA
    atomicAdd(ptr,val);
  #else
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    *ptr += val;
  #endif
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void atomicAddWrapper(FLOATORDBL *ptr, double val) {
  #ifdef __NVCC__ // If compiling for CUDA
    atomicAdd(ptr,val);
  #else
    #ifdef _OPENMP
    #pragma omp atomic
    #endif
    *ptr += val;
  #endif
}

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
void * reallocWrapper(void *ptr_old, unsigned long long oldsize, unsigned long long newsize) {
  #ifdef __NVCC__ // If compiling for CUDA
  void *ptr_new = malloc(newsize); // GPUs don't have a realloc function
  memcpy(ptr_new,ptr_old,oldsize);
  free(ptr_old);
  #else
  void *ptr_new = realloc(ptr_old,newsize);
  if(!ptr_new) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");
  #endif
  return ptr_new;
}

void unitcrossprod(FLOATORDBL *a, FLOATORDBL *b, FLOATORDBL *c) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
  FLOATORDBL norm = SQRT(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
  for(int idx=0;idx<3;idx++) c[idx] /= norm;
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void xyztoXYZ(FLOATORDBL * const r, FLOATORDBL const theta, FLOATORDBL const phi, FLOATORDBL * const r_out) {
  // If input r is given in the simulation (x,y,z) frame, then output r_out is that point's coordinates in
  // the (X,Y,Z) light collector frame, where Z is pointing in the direction denoted by the spherical
  // coordinates theta and phi and X lies in the xy plane.
  r_out[0] =            SIN(phi)*r[0] -            COS(phi)*r[1];
  r_out[1] = COS(theta)*COS(phi)*r[0] + COS(theta)*SIN(phi)*r[1] - SIN(theta)*r[2];
  r_out[2] = SIN(theta)*COS(phi)*r[0] + SIN(theta)*SIN(phi)*r[1] + COS(theta)*r[2];
}

#ifdef __NVCC__ // If compiling for CUDA
__device__ __host__
#endif
void axisrotate(FLOATORDBL const * const r, FLOATORDBL const * const u, FLOATORDBL const theta, FLOATORDBL * const r_out) {
  // Rotate the point r an angle theta around vector u, storing result in r_out
  FLOATORDBL st = SIN(theta), ct = COS(theta);
  
  r_out[0] = (u[0]*u[0]*(1-ct) +      ct)*r[0] + (u[0]*u[1]*(1-ct) - u[2]*st)*r[1] + (u[0]*u[2]*(1-ct) + u[1]*st)*r[2];
  r_out[1] = (u[1]*u[0]*(1-ct) + u[2]*st)*r[0] + (u[1]*u[1]*(1-ct) +      ct)*r[1] + (u[1]*u[2]*(1-ct) - u[0]*st)*r[2];
  r_out[2] = (u[2]*u[0]*(1-ct) - u[1]*st)*r[0] + (u[2]*u[1]*(1-ct) + u[0]*st)*r[1] + (u[2]*u[2]*(1-ct) +      ct)*r[2];
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
long binaryTreeSearch(FLOATORDBL rand, long N, FLOATORDBL *D) {
  long d = N-1; // Number of elements in the current section to be searched, minus one
  long j = d/2; // Index of the middle element of the current section
  while(!(D[j] < rand && D[j+1] >= rand)) { // Binary tree search
    if(D[j] >= rand) {
      j = j + (d-2)/4 - d/2;
      d = d/2 - 1;
    } else {
      d = d - d/2 - 1;
      j = j + d/2 + 1;
    }
  }
  return j;
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void getNewj(struct geometry const *G, struct photon *P) {
  P->j = ((P->i[2] < 0)? 0: ((P->i[2] >= G->n[2])? G->n[2]-1: (long)FLOOR(P->i[2])))*G->n[0]*G->n[1] +
         ((P->i[1] < 0)? 0: ((P->i[1] >= G->n[1])? G->n[1]-1: (long)FLOOR(P->i[1])))*G->n[0]         +
         ((P->i[0] < 0)? 0: ((P->i[0] >= G->n[0])? G->n[0]-1: (long)FLOOR(P->i[0]))); // Index values are restrained to integers in the interval [0,n-1]
}

#ifdef __NVCC__ // If compiling for CUDA
void createDeviceStructs(struct geometry const *G, struct geometry **G_devptr,
                         struct source const *B, struct source **B_devptr,
                         struct lightCollector const *LC, struct lightCollector **LC_devptr,
                         struct paths *Pa, struct paths **Pa_devptr,
                         struct outputs *O, struct outputs **O_devptr,
                         struct depositionCriteria *DC, struct depositionCriteria **DC_devptr,
                         int nM, long L,
                         size_t size_smallArrays,
                         struct debug *D, struct debug **D_devptr) {
  bool matchedInterfaces = true;
  for(long iM=1;iM<nM;iM++) matchedInterfaces &= G->RIv[iM] == G->RIv[iM-1];

  // Allocate and copy geometry struct, including smallArrays
  struct geometry G_tempvar = *G;
  gpuErrchk(cudaMalloc(&G_tempvar.muav,size_smallArrays)); // This is to allocate all of smallArrays
  gpuErrchk(cudaMemcpy( G_tempvar.muav,G->muav,size_smallArrays,cudaMemcpyHostToDevice)); // And for copying all of it to global device memory
  gpuErrchk(cudaMalloc(&G_tempvar.M, L*sizeof(unsigned char)));
  gpuErrchk(cudaMemcpy( G_tempvar.M, G->M, L*sizeof(unsigned char),cudaMemcpyHostToDevice));
  gpuErrchk(cudaMalloc(&G_tempvar.interfaceNormals, (matchedInterfaces? 1:2*L)*sizeof(float)));
  gpuErrchk(cudaMemcpy( G_tempvar.interfaceNormals, G->interfaceNormals, (matchedInterfaces? 1:2*L)*sizeof(float),cudaMemcpyHostToDevice));

  gpuErrchk(cudaMalloc(G_devptr, sizeof(struct geometry)));
  gpuErrchk(cudaMemcpy(*G_devptr,&G_tempvar,sizeof(struct geometry),cudaMemcpyHostToDevice));

  // Allocate and copy beam struct
  struct source B_tempvar = *B;
  if(B->S) {
    gpuErrchk(cudaMalloc(&B_tempvar.S, (L+1)*sizeof(FLOATORDBL)));
    gpuErrchk(cudaMemcpy(B_tempvar.S, B->S, (L+1)*sizeof(FLOATORDBL),cudaMemcpyHostToDevice));
  }
  gpuErrchk(cudaMalloc(B_devptr, sizeof(struct source)));
  gpuErrchk(cudaMemcpy(*B_devptr,&B_tempvar,sizeof(struct source),cudaMemcpyHostToDevice));

  // Allocate and copy deposition criteria struct
  gpuErrchk(cudaMalloc(DC_devptr, sizeof(struct depositionCriteria)));
  gpuErrchk(cudaMemcpy(*DC_devptr,DC,sizeof(struct depositionCriteria),cudaMemcpyHostToDevice));

  // Allocate and copy lightCollector struct
  gpuErrchk(cudaMalloc(LC_devptr, sizeof(struct lightCollector)));
  gpuErrchk(cudaMemcpy(*LC_devptr,LC,sizeof(struct lightCollector),cudaMemcpyHostToDevice));

  // Allocate and copy paths struct
  struct paths Pa_tempvar = *Pa;
  if(Pa->pathsSize) gpuErrchk(cudaMalloc(&Pa_tempvar.data, 4*Pa->pathsSize*sizeof(FLOATORDBL)));
  gpuErrchk(cudaMalloc(Pa_devptr, sizeof(struct paths)));
  gpuErrchk(cudaMemcpy(*Pa_devptr,&Pa_tempvar,sizeof(struct paths),cudaMemcpyHostToDevice));

  // Allocate and copy outputs struct
  struct outputs O_tempvar = *O;
  // Note in the following that memset and cudaMemset is for integers (really, signed chars) and only works here because four signed char zeros after each other have the same bit representation as a double- or single-precision floating point zero
  if(O->NFR) {
    gpuErrchk(cudaMalloc(&O_tempvar.NFR, L*sizeof(double)));
    gpuErrchk(cudaMemset(O_tempvar.NFR,0,L*sizeof(double)));
  }
  if(O->image) {
    gpuErrchk(cudaMalloc(&O_tempvar.image, LC->res[0]*LC->res[0]*LC->res[1]*sizeof(double)));
    gpuErrchk(cudaMemset(O_tempvar.image,0,LC->res[0]*LC->res[0]*LC->res[1]*sizeof(double)));
  }
  if(O->FF) {
    gpuErrchk(cudaMalloc(&O_tempvar.FF, G->farFieldRes*G->farFieldRes*sizeof(double)));
    gpuErrchk(cudaMemset(O_tempvar.FF,0,G->farFieldRes*G->farFieldRes*sizeof(double)));
  }
  if(O->NI_xpos) {
    gpuErrchk(cudaMalloc(&O_tempvar.NI_xpos, G->n[1]*G->n[2]*sizeof(double)));
    gpuErrchk(cudaMemset(O_tempvar.NI_xpos,0,G->n[1]*G->n[2]*sizeof(double)));
  }
  if(O->NI_xneg) {
    gpuErrchk(cudaMalloc(&O_tempvar.NI_xneg, G->n[1]*G->n[2]*sizeof(double)));
    gpuErrchk(cudaMemset(O_tempvar.NI_xneg,0,G->n[1]*G->n[2]*sizeof(double)));
  }
  if(O->NI_ypos) {
    gpuErrchk(cudaMalloc(&O_tempvar.NI_ypos, G->n[0]*G->n[2]*sizeof(double)));
    gpuErrchk(cudaMemset(O_tempvar.NI_ypos,0,G->n[0]*G->n[2]*sizeof(double)));
  }
  if(O->NI_yneg) {
    gpuErrchk(cudaMalloc(&O_tempvar.NI_yneg, G->n[0]*G->n[2]*sizeof(double)));
    gpuErrchk(cudaMemset(O_tempvar.NI_yneg,0,G->n[0]*G->n[2]*sizeof(double)));
  }
  if(O->NI_zpos) {
    gpuErrchk(cudaMalloc(&O_tempvar.NI_zpos, G->n[0]*G->n[1]*sizeof(double)));
    gpuErrchk(cudaMemset(O_tempvar.NI_zpos,0,G->n[0]*G->n[1]*sizeof(double)));
  }
  if(O->NI_zneg) {
    gpuErrchk(cudaMalloc(&O_tempvar.NI_zneg, (G->boundaryType == 2? KILLRANGE*KILLRANGE: 1)*G->n[0]*G->n[1]*sizeof(double)));
    gpuErrchk(cudaMemset(O_tempvar.NI_zneg,0,(G->boundaryType == 2? KILLRANGE*KILLRANGE: 1)*G->n[0]*G->n[1]*sizeof(double)));
  }
  gpuErrchk(cudaMalloc(O_devptr, sizeof(struct outputs)));
  gpuErrchk(cudaMemcpy(*O_devptr,&O_tempvar,sizeof(struct outputs),cudaMemcpyHostToDevice));
  
  // Allocate and copy debug struct
  struct debug D_tempvar = *D;
  gpuErrchk(cudaMalloc(D_devptr, sizeof(struct debug)));
  gpuErrchk(cudaMemcpy(*D_devptr,&D_tempvar,sizeof(struct debug),cudaMemcpyHostToDevice));
}
#endif

#ifdef __NVCC__ // If compiling for CUDA
void retrieveAndFreeDeviceStructs(struct geometry const *G, struct geometry *G_dev,
                                  struct source const *B, struct source *B_dev,
                                  struct lightCollector const *LC, struct lightCollector *LC_dev,
                                  struct paths *Pa, struct paths *Pa_dev,
                                  struct outputs *O, struct outputs *O_dev,
                                  struct depositionCriteria *DC_dev,
                                  long L,
                                  struct debug *D, struct debug *D_dev) {
  struct geometry G_temp; gpuErrchk(cudaMemcpy(&G_temp, G_dev, sizeof(struct geometry),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaFree(G_temp.muav)); // This frees all of smallArrays in the global memory on the device
  gpuErrchk(cudaFree(G_temp.M));
  gpuErrchk(cudaFree(G_temp.interfaceNormals));
  gpuErrchk(cudaFree(G_dev));

  struct source B_temp; gpuErrchk(cudaMemcpy(&B_temp, B_dev, sizeof(struct source),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaFree(B_temp.S));
  gpuErrchk(cudaFree(B_dev));
  
  gpuErrchk(cudaFree(DC_dev));

  gpuErrchk(cudaFree(LC_dev));
  
  struct paths Pa_temp = *Pa;
  gpuErrchk(cudaMemcpy(Pa, Pa_dev, sizeof(struct paths),cudaMemcpyDeviceToHost));
  if(Pa->data) {
    gpuErrchk(cudaMemcpy(Pa_temp.data, Pa->data, 4*Pa->pathsElems*sizeof(FLOATORDBL),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(Pa->data));
    Pa->data = Pa_temp.data;
  }
  gpuErrchk(cudaFree(Pa_dev));

  struct outputs O_temp; gpuErrchk(cudaMemcpy(&O_temp, O_dev, sizeof(struct outputs),cudaMemcpyDeviceToHost));
  O->nPhotons = O_temp.nPhotons;
  O->nPhotonsCollected = O_temp.nPhotonsCollected;
  if(O->NFR) {
    gpuErrchk(cudaMemcpy(O->NFR, O_temp.NFR, L*sizeof(FLOATORDBL),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.NFR));
  }
  if(O->image) {
    gpuErrchk(cudaMemcpy(O->image, O_temp.image, LC->res[0]*LC->res[0]*LC->res[1]*sizeof(FLOATORDBL),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.image));
  }
  if(O->FF) {
    gpuErrchk(cudaMemcpy(O->FF, O_temp.FF, G->farFieldRes*G->farFieldRes*sizeof(FLOATORDBL),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.FF));
  }
  if(O->NI_xpos) {
    gpuErrchk(cudaMemcpy(O->NI_xpos, O_temp.NI_xpos, G->n[1]*G->n[2]*sizeof(FLOATORDBL),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.NI_xpos));
  }
  if(O->NI_xneg) {
    gpuErrchk(cudaMemcpy(O->NI_xneg, O_temp.NI_xneg, G->n[1]*G->n[2]*sizeof(FLOATORDBL),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.NI_xneg));
  }
  if(O->NI_ypos) {
    gpuErrchk(cudaMemcpy(O->NI_ypos, O_temp.NI_ypos, G->n[0]*G->n[2]*sizeof(FLOATORDBL),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.NI_ypos));
  }
  if(O->NI_yneg) {
    gpuErrchk(cudaMemcpy(O->NI_yneg, O_temp.NI_yneg, G->n[0]*G->n[2]*sizeof(FLOATORDBL),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.NI_yneg));
  }
  if(O->NI_zpos) {
    gpuErrchk(cudaMemcpy(O->NI_zpos, O_temp.NI_zpos, G->n[0]*G->n[1]*sizeof(FLOATORDBL),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.NI_zpos));
  }
  if(O->NI_zneg) {
    gpuErrchk(cudaMemcpy(O->NI_zneg, O_temp.NI_zneg, (G->boundaryType == 2? KILLRANGE*KILLRANGE: 1)*G->n[0]*G->n[1]*sizeof(FLOATORDBL),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.NI_zneg));
  }
  gpuErrchk(cudaFree(O_dev));
  
  gpuErrchk(cudaMemcpy(D, D_dev, sizeof(struct debug),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaFree(D_dev));
}
#endif

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void deletePhotonPath(struct paths * const Pa) {
  for(; Pa->pathsElems > 1; Pa->pathsElems--) {
    if(ISNAN(Pa->data[4*(Pa->pathsElems - 1)]) && ISNAN(Pa->data[4*(Pa->pathsElems - 2)])) { // 4 entries per column. Two previous columns NaN means the start of the path.
      Pa->pathsElems--;
      Pa->pathsElems--;
      break;
    }
  }
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void addToPhotonPath(struct photon * const P, struct paths * const Pa, struct geometry const * const G, char NaNcolumns) {
  if(Pa->pathsElems == Pa->pathsSize - NaNcolumns) {
    #ifdef __NVCC__ // If compiling for CUDA
    Pa->pathsElems = -1; // Buffer not large enough, has to be resized by host
    #else
    Pa->pathsSize *= 2; // double the record's size
    Pa->data = (FLOATORDBL *)reallocWrapper(Pa->data,2*Pa->pathsSize*sizeof(FLOATORDBL),4*Pa->pathsSize*sizeof(FLOATORDBL));
    #endif
  }
  if(Pa->pathsElems != -1) {
    for(int i = 0; i < NaNcolumns; i++) {
      Pa->data[4*Pa->pathsElems  ] = NAN;
      Pa->data[4*Pa->pathsElems+1] = NAN;
      Pa->data[4*Pa->pathsElems+2] = NAN;
      Pa->data[4*Pa->pathsElems+3] = NAN;
      Pa->pathsElems++;
    }
    Pa->data[4*Pa->pathsElems  ] = (P->i[0] - G->n[0]/2.0f)*G->d[0]; // Store x
    Pa->data[4*Pa->pathsElems+1] = (P->i[1] - G->n[1]/2.0f)*G->d[1]; // Store y
    Pa->data[4*Pa->pathsElems+2] = (P->i[2]               )*G->d[2]; // Store z
    Pa->data[4*Pa->pathsElems+3] = P->weight;                        // Store photon weight
    Pa->pathsElems++;
  }
  Pa->pathStartedThisPhoton = true;
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void updatePaths(struct photon * const P, struct paths * const Pa, struct geometry const * const G, struct depositionCriteria * DC, bool photonTeleported) {
  if(Pa->nExamplePhotonPathsFinished < Pa->nExamplePaths) {
    if(DC->evaluateCriteriaAtEndOfLife || depositionCriteriaMet(P,DC)) {
      char NaNcolumns = Pa->pathStartedThisPhoton? (photonTeleported? 1: 0): 2;
      addToPhotonPath(P,Pa,G,NaNcolumns);
    }
  }
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void launchPhoton(struct photon * const P, struct source const * const B, struct geometry const * const G, struct paths * const Pa, struct depositionCriteria *DC, bool * abortingPtr, struct debug * D) {
  FLOATORDBL X,Y,r,phi,tanphiX,tanphiY,costheta,sintheta;
  long   j,idx;
  FLOATORDBL target[3]={0},w0[3];
  
  P->sameVoxel = false;
  P->weight = 1;
  P->recordElems = 0;
  P->killed_escaped_collected = 0; // Default state to killed
  long launchAttempts = 0;
  do{
    if(B->S) { // If a 3D source distribution was defined
      // ... then search the cumulative distribution function via binary tree method to find the voxel to start the photon in
      j = binaryTreeSearch(RandomNum,G->n[0]*G->n[1]*G->n[2],B->S);
      P->i[0] = j%G->n[0]         + 1 - RandomNum;
      P->i[1] = j/G->n[0]%G->n[1] + 1 - RandomNum;
      P->i[2] = j/G->n[0]/G->n[1] + 1 - RandomNum;
      costheta = 1 - 2*RandomNum;
      sintheta = SQRT(1 - costheta*costheta);
      phi = 2*PI*RandomNum;
      P->u[0] = sintheta*COS(phi);
      P->u[1] = sintheta*SIN(phi);
      P->u[2] = costheta;
      getNewj(G,P);
      P->RI = G->RIv[G->M[P->j]]; // Set the refractive index
      P->time = 0;
    } else switch (B->beamType) {
      case 0: // pencil beam
        P->i[0] = (B->focus[0] - B->focus[2]*B->u[0]/B->u[2])/G->d[0] + G->n[0]/2.0f;
        P->i[1] = (B->focus[1] - B->focus[2]*B->u[1]/B->u[2])/G->d[1] + G->n[1]/2.0f;
        P->i[2] = 0;
        for(idx=0;idx<3;idx++) P->u[idx] = B->u[idx];
        getNewj(G,P);
        P->RI = G->RIv[G->M[P->j]]; // Set the refractive index
        P->time = -P->RI/C*SQRT(SQR((P->i[0] - G->n[0]/2.0f)*G->d[0] - B->focus[0]) +
                                SQR((P->i[1] - G->n[1]/2.0f)*G->d[1] - B->focus[1]) +
                                SQR((P->i[2]               )*G->d[2] - B->focus[2])); // Starting time is set so that the wave crosses the focal plane at time = 0
        break;
      case 1: // isotropically emitting point source
        r = RandomNum;
        P->i[0] = (B->focus[0] + B->u[0]*(r-0.5)*B->emitterLength)/G->d[0] + G->n[0]/2.0f;
        P->i[1] = (B->focus[1] + B->u[1]*(r-0.5)*B->emitterLength)/G->d[1] + G->n[1]/2.0f;
        P->i[2] = (B->focus[2] + B->u[2]*(r-0.5)*B->emitterLength)/G->d[2];
        costheta = 1 - 2*RandomNum;
        sintheta = SQRT(1 - costheta*costheta);
        phi = 2*PI*RandomNum;
        P->u[0] = sintheta*COS(phi);
        P->u[1] = sintheta*SIN(phi);
        P->u[2] = costheta;
        getNewj(G,P);
        P->RI = G->RIv[G->M[P->j]]; // Set the refractive index
        P->time = 0;
        break;
      case 2: // infinite plane wave
        P->i[0] = ((G->boundaryType==1 || G->boundaryType==3)? 1: KILLRANGE)*G->n[0]*(RandomNum-0.5f) + G->n[0]/2.0f; // Generates a random ix coordinate within the cuboid
        P->i[1] = ((G->boundaryType==1 || G->boundaryType==3)? 1: KILLRANGE)*G->n[1]*(RandomNum-0.5f) + G->n[1]/2.0f; // Generates a random iy coordinate within the cuboid
        P->i[2] = 0;
        for(idx=0;idx<3;idx++) P->u[idx] = B->u[idx];
        getNewj(G,P);
        P->RI = G->RIv[G->M[P->j]]; // Set the refractive index
        P->time = P->RI/C*((P->i[0] - G->n[0]/2.0f)*G->d[0]*B->u[0] +
                           (P->i[1] - G->n[1]/2.0f)*G->d[1]*B->u[1] +
                           (P->i[2]               )*G->d[2]*B->u[2]); // Starting time is set so that the wave crosses (x=0,y=0,z=0) at time = 0
        break;
      case 3: // Laguerre-Gaussian LG01 beam
        phi     = RandomNum*2*PI;
        axisrotate(B->v,B->u,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
        r    = B->FPIDwidth1*SQRT(((FLOATORDBL)gsl_sf_lambert_Wm1(-RandomNum*EXP(-1.0f))+1)/(-2))/1.50087f; // for target calculation
        for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
        phi     = B->AIDwidth1*SQRT(((FLOATORDBL)gsl_sf_lambert_Wm1(-RandomNum*EXP(-1.0f))+1)/(-2))/1.50087f; // for trajectory calculation. The sqrt is valid within paraxial approximation.
        axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
        P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0f; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
        P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0f;
        P->i[2] = 0;
        getNewj(G,P);
        P->RI = G->RIv[G->M[P->j]]; // Set the refractive index
        P->time = -P->RI/C*SQRT(SQR((P->i[0] - G->n[0]/2.0f)*G->d[0] - target[0]) +
                                SQR((P->i[1] - G->n[1]/2.0f)*G->d[1] - target[1]) +
                                SQR((P->i[2]               )*G->d[2] - target[2])); // Starting time is set so that the wave crosses the focal plane at time = 0
        break;
      case 4: // Radial
        // Near Field
        axisrotate(B->v,B->u,RandomNum*2*PI,w0); // w0 unit vector now points in the direction from focus center point to ray target point
        if(*B->FPIDdist1 == -1) { // Top-hat radial distribution
          r = B->FPIDwidth1*SQRT(RandomNum); // for target calculation
        } else if(*B->FPIDdist1 == -2) { // Gaussian radial distribution
          r = B->FPIDwidth1*SQRT(-0.5f*LOG(RandomNum)); // for target calculation
        } else { // Custom distribution
          r = B->FPIDwidth1*(binaryTreeSearch(RandomNum,B->L_FPID1-1,B->FPIDdist1)+RandomNum)/(B->L_FPID1-1);
        }
        for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
        
        // Far Field
        axisrotate(B->v,B->u,RandomNum*2*PI,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.

        if(*B->AIDdist1 == -1) { // Top-hat radial distribution
          phi = ATAN(TAN(B->AIDwidth1)*SQRT(RandomNum)); // for trajectory calculation. The sqrt is valid within paraxial approximation.
        } else if(*B->AIDdist1 == -2) { // Gaussian radial distribution
          phi = ATAN(TAN(B->AIDwidth1)*SQRT(-0.5f*LOG(RandomNum))); // for trajectory calculation. The sqrt is valid within paraxial approximation.
        } else if(*B->AIDdist1 == -3) { // Lambertian
          phi = ASIN(SQRT(RandomNum));
        } else { // Custom distribution
          phi = ATAN(TAN(B->AIDwidth1)*(binaryTreeSearch(RandomNum,B->L_AID1-1,B->AIDdist1)+RandomNum)/(B->L_AID1-1));
        }
        axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
        
        P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0f; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
        P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0f;
        P->i[2] = 0;
        getNewj(G,P);
        P->RI = G->RIv[G->M[P->j]]; // Set the refractive index
        P->time = -P->RI/C*SQRT(SQR((P->i[0] - G->n[0]/2.0f)*G->d[0] - target[0]) +
                                SQR((P->i[1] - G->n[1]/2.0f)*G->d[1] - target[1]) +
                                SQR((P->i[2]               )*G->d[2] - target[2])); // Starting time is set so that the wave crosses the focal plane at time = 0
        break;
      case 5: // X/Y
        // Near Field
        if(*B->FPIDdist1 == -1) { // Top-hat X distribution
          X = B->FPIDwidth1*(RandomNum*2-1); // for target calculation
        } else if(*B->FPIDdist1 == -2) { // Gaussian X distribution
          X = B->FPIDwidth1*SQRT(-0.5f*LOG(RandomNum))*COS(2*PI*RandomNum); // Box-Muller transform, for target calculation
        } else { // Custom X distribution
          X = B->FPIDwidth1*((binaryTreeSearch(RandomNum,B->L_FPID1-1,B->FPIDdist1)+RandomNum)/(B->L_FPID1-1)*2-1);
        }
        if(*B->FPIDdist2 == -1) { // Top-hat Y distribution
          Y = B->FPIDwidth2*(RandomNum*2-1); // for target calculation
        } else if(*B->FPIDdist2 == -2) { // Gaussian Y distribution
          Y = B->FPIDwidth2*SQRT(-0.5f*LOG(RandomNum))*COS(2*PI*RandomNum); // Box-Muller transform, for target calculation
        } else { // Custom distribution
          Y = B->FPIDwidth2*((binaryTreeSearch(RandomNum,B->L_FPID2-1,B->FPIDdist2)+RandomNum)/(B->L_FPID2-1)*2-1);
        }
        for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + X*B->v[idx] + Y*B->w[idx];

        // Far Field
        if(*B->AIDdist1 == -3) { // Lambertian
          axisrotate(B->v,B->u,RandomNum*2*PI,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction
          axisrotate(B->u,w0,ASIN(SQRT(RandomNum)),P->u); // ray propagation direction is found by rotating beam center axis around w0
        } else {
          if(*B->AIDdist1 == -1) { // Top-hat phiX distribution
            tanphiX = TAN(B->AIDwidth1)*(RandomNum*2-1); // for trajectory calculation
          } else if(*B->AIDdist1 == -2) { // Gaussian phiX distribution
            tanphiX = TAN(B->AIDwidth1)*SQRT(-0.5f*LOG(RandomNum))*COS(2*PI*RandomNum); // Box-Muller transform, for trajectory calculation
          } else { // Custom phi_X distribution
            tanphiX = TAN(B->AIDwidth1)*((binaryTreeSearch(RandomNum,B->L_AID1,B->AIDdist1)+RandomNum)/B->L_AID1*2-1);
          }
          if(*B->AIDdist2 == -1) { // Top-hat phiY distribution
            tanphiY = TAN(B->AIDwidth2)*(RandomNum*2-1); // for trajectory calculation
          } else if(*B->AIDdist2 == -2) { // Gaussian phiY distribution
            tanphiY = TAN(B->AIDwidth2)*SQRT(-0.5f*LOG(RandomNum))*COS(2*PI*RandomNum); // Box-Muller transform, for trajectory calculation
          } else { // Custom distribution
            tanphiY = TAN(B->AIDwidth2)*((binaryTreeSearch(RandomNum,B->L_AID2,B->AIDdist2)+RandomNum)/B->L_AID2*2-1);
          }
          axisrotate(B->v,B->u,ATAN2(tanphiX,tanphiY),w0); // w0 is now orthogonal to both beam propagation axis and ray propagation axis
          axisrotate(B->u,w0,ATAN(SQRT(tanphiX*tanphiX + tanphiY*tanphiY)),P->u); // ray propagation direction is found by rotating beam center axis around w0
        }
        P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0f; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
        P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0f;
        P->i[2] = 0;
        getNewj(G,P);
        P->RI = G->RIv[G->M[P->j]]; // Set the refractive index
        P->time = -P->RI/C*SQRT(SQR((P->i[0] - G->n[0]/2.0f)*G->d[0] - target[0]) +
                                SQR((P->i[1] - G->n[1]/2.0f)*G->d[1] - target[1]) +
                                SQR((P->i[2]               )*G->d[2] - target[2])); // Starting time is set so that the wave crosses the focal plane at time = 0
        break;
    }

    switch (G->boundaryType) {
      case 0:
        P->alive = (FABS(P->i[0]/G->n[0] - 1.0f/2) <  KILLRANGE/2.0f &&
                    FABS(P->i[1]/G->n[1] - 1.0f/2) <  KILLRANGE/2.0f &&
                    FABS(P->i[2]/G->n[2] - 1.0f/2) <  KILLRANGE/2.0f);
        break;
      case 1:
        P->alive = P->i[0] < G->n[0] && P->i[0] >= 0 &&
                   P->i[1] < G->n[1] && P->i[1] >= 0 &&
                   P->i[2] < G->n[2] && P->i[2] >= 0;
        break;
      case 2:
        P->alive = (FABS(P->i[0]/G->n[0] - 1.0f/2) <  KILLRANGE/2.0f &&
                    FABS(P->i[1]/G->n[1] - 1.0f/2) <  KILLRANGE/2.0f &&
                         P->i[2]/G->n[2] - 1.0f/2  <  KILLRANGE/2.0f &&
                         P->i[2]                   >= 0);
        break;
      case 3:
        P->alive = P->i[0] < G->n[0] && P->i[0] >= 0 &&
                   P->i[1] < G->n[1] && P->i[1] >= 0 &&
                   P->i[2] < G->n[2] && P->i[2] >= 0;
        break;
    }
  } while(!P->alive && ++launchAttempts < 1000000); // If photon happened to be initialized outside the volume in which it is allowed to travel, we try again unless it's happened a million times in a row.

  if(launchAttempts >= 1000000) *abortingPtr = true;
  
  // Calculate distances to next voxel boundary planes
  for(idx=0;idx<3;idx++) P->D[idx] = P->u[idx]? (FLOOR(P->i[idx]) + (P->u[idx]>0) - P->i[idx])*G->d[idx]/P->u[idx] : INFINITY;

  P->insideVolume = P->i[0] < G->n[0] && P->i[0] >= 0 &&
                    P->i[1] < G->n[1] && P->i[1] >= 0 &&
                    P->i[2] < G->n[2] && P->i[2] >= 0;

  P->stepLeft  = -LOG(RandomNum);
  
  P->scatterings = P->refractions = P->reflections = P->interfaceTransitions = 0;

  #ifdef __NVCC__ // If compiling for CUDA
  if(!threadIdx.x && !blockIdx.x)
  #elif defined(_OPENMP)
  #pragma omp master
  #endif
  {
    Pa->pathStartedThisPhoton = false;
    updatePaths(P,Pa,G,DC,false);
  }
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void formImage(struct photon * const P, struct geometry const * const G, struct lightCollector const * const LC, struct depositionCriteria *DC, struct outputs *O, unsigned long long * nPhotonsCollectedPtr) {
  FLOATORDBL U[3];
  xyztoXYZ(P->u,LC->theta,LC->phi,U); // U is now the photon trajectory in basis of detection frame (X,Y,Z)
  
  if(U[2] < 0) { // If the Z component of U is negative then the photon is moving towards the light collector plane
    FLOATORDBL resc[3] = {(P->i[0] - G->n[0]/2.0f)*G->d[0] - LC->r[0],
                          (P->i[1] - G->n[1]/2.0f)*G->d[1] - LC->r[1],
                          (P->i[2]               )*G->d[2] - LC->r[2]}; // Photon position relative to the light collector focal plane center when it escapes the cuboid, in the (x,y,z) basis
    
    FLOATORDBL Resc[3];
    xyztoXYZ(resc,LC->theta,LC->phi,Resc); // Resc is now the photon position in the light collector frame (X,Y,Z) when it escapes the cuboid
    
    if(Resc[2] > 0) { // If the Z coordinate in the light collector frame is positive (negative means the photon is already on the wrong side)
      FLOATORDBL RLCP[2]; // XY coordinates of the point where the photon crosses the light collector plane
      RLCP[0] = Resc[0] - Resc[2]*U[0]/U[2];
      RLCP[1] = Resc[1] - Resc[2]*U[1]/U[2];
      
      FLOATORDBL distLCP = SQRT(RLCP[0]*RLCP[0] + RLCP[1]*RLCP[1]); // Distance between light collector center and the point where the photon crosses the light collector plane
      
      if(distLCP < LC->diam/2) { // If the distance is less than the radius of the light collector
        if(ISFINITE(LC->f)) { // If the light collector is an objective lens
          FLOATORDBL RImP[2]; // Back-propagated position that the photon would have had in the object plane if propagating freely. This corresponds to where the photon will end up in the image plane for magnification 1x.
          RImP[0] = RLCP[0] + LC->f*U[0]/U[2];
          RImP[1] = RLCP[1] + LC->f*U[1]/U[2];
          FLOATORDBL distImP = SQRT(RImP[0]*RImP[0] + RImP[1]*RImP[1]);
          if(distImP < LC->FSorNA/2) { // If the photon is coming from the area within the Field Size
            long Xindex = (long)(LC->res[0]*(RImP[0]/LC->FSorNA + 1.0f/2));
            long Yindex = (long)(LC->res[0]*(RImP[1]/LC->FSorNA + 1.0f/2));
            long timeindex = LC->res[1] > 1? min(LC->res[1]-1,max(0L,(long)(1+(LC->res[1]-2)*(P->time - (Resc[2] - LC->f)/U[2]*P->RI/C - LC->tStart)/(LC->tEnd - LC->tStart)))): 0; // If we are not measuring time-resolved, LC->res[1] == 1
            P->killed_escaped_collected = 2; // Collected
            if(depositionCriteriaMet(P,DC)) {
              atomicAddWrapper(&O->image[Xindex               +
                                         Yindex   *LC->res[0] +
                                         timeindex*LC->res[0]*LC->res[0]],P->weight);
              atomicAddWrapperULL(nPhotonsCollectedPtr,1);
            }
          }
        } else { // If the light collector is a fiber tip
          FLOATORDBL thetaLCFF = ATAN(-SQRT(U[0]*U[0] + U[1]*U[1])/U[2]); // Light collector far field polar angle
          if(thetaLCFF < ASIN(min(1.0f,LC->FSorNA))) { // If the photon has an angle within the fiber's NA acceptance
            long timeindex = LC->res[1] > 1? min(LC->res[1]-1,max(0L,(long)(1+(LC->res[1]-2)*(P->time - Resc[2]/U[2]*P->RI/C - LC->tStart)/(LC->tEnd - LC->tStart)))): 0; // If we are not measuring time-resolved, LC->res[1] == 1
            P->killed_escaped_collected = 2; // Collected
            if(depositionCriteriaMet(P,DC)) {
              atomicAddWrapper(&O->image[timeindex],P->weight);
              atomicAddWrapperULL(nPhotonsCollectedPtr,1);
            }
          }
        }
      }
    }
  }
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void formFarField(struct photon const * const P, struct geometry const *G, struct outputs const *O) {
  FLOATORDBL theta = (1-FLOATORDBLEPS)*ACOS(P->u[2]); // The (1-EPS) factor is to ensure that photons exiting with theta = PI will be stored correctly
  FLOATORDBL phi_shifted = (1-FLOATORDBLEPS)*(PI + ATAN2(P->u[1],P->u[0])); // Here it's to handle the case of phi = +PI
  atomicAddWrapper(&O->FF[(long)FLOOR(theta/PI*G->farFieldRes) + G->farFieldRes*(long)FLOOR(phi_shifted/(2*PI)*G->farFieldRes)],P->weight);
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void formEdgeFluxes(struct photon const * const P, struct geometry const * const G, struct outputs const *O) {
  if(G->boundaryType == 1) {
    if(P->i[2] < 0)             atomicAddWrapper(&O->NI_zneg[(long)P->i[0] + G->n[0]*(long)P->i[1]],P->weight);
    else if(P->i[2] >= G->n[2]) atomicAddWrapper(&O->NI_zpos[(long)P->i[0] + G->n[0]*(long)P->i[1]],P->weight);
    else if(P->i[1] < 0)        atomicAddWrapper(&O->NI_yneg[(long)P->i[0] + G->n[0]*(long)P->i[2]],P->weight);
    else if(P->i[1] >= G->n[1]) atomicAddWrapper(&O->NI_ypos[(long)P->i[0] + G->n[0]*(long)P->i[2]],P->weight);
    else if(P->i[0] < 0)        atomicAddWrapper(&O->NI_xneg[(long)P->i[1] + G->n[1]*(long)P->i[2]],P->weight);
    else if(P->i[0] >= G->n[0]) atomicAddWrapper(&O->NI_xpos[(long)P->i[1] + G->n[1]*(long)P->i[2]],P->weight);
  } else if(G->boundaryType == 2) {
    if(P->i[2] < 0)             atomicAddWrapper(&O->NI_zneg[(long)(P->i[0] + G->n[0]*(KILLRANGE-1)/2.0f) + (KILLRANGE*G->n[0])*((long)(P->i[1] + G->n[1]*(KILLRANGE-1)/2.0f))],P->weight);
  } else { // boundaryType == 3
    if(P->i[2] < 0)             atomicAddWrapper(&O->NI_zneg[(long)P->i[0] + G->n[0]*(long)P->i[1]],P->weight);
    else if(P->i[2] >= G->n[2]) atomicAddWrapper(&O->NI_zpos[(long)P->i[0] + G->n[0]*(long)P->i[1]],P->weight);
  }
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void checkEscape(struct photon * const P, struct paths *Pa, struct geometry const * const G, struct lightCollector const * const LC,
        struct outputs *O, struct depositionCriteria *DC, unsigned long long * nPhotonsCollectedPtr) {
  bool escaped = false;
  switch (G->boundaryType) {
    case 0:
      P->alive = (FABS(P->i[0]/G->n[0] - 1.0f/2) <  KILLRANGE/2.0f &&
                  FABS(P->i[1]/G->n[1] - 1.0f/2) <  KILLRANGE/2.0f &&
                  FABS(P->i[2]/G->n[2] - 1.0f/2) <  KILLRANGE/2.0f);
      break;
    case 1:
      P->alive = P->i[0] < G->n[0] && P->i[0] >= 0 &&
                 P->i[1] < G->n[1] && P->i[1] >= 0 &&
                 P->i[2] < G->n[2] && P->i[2] >= 0;
      escaped = !P->alive && P->RI == 1;
      break;
    case 2:
      P->alive = (FABS(P->i[0]/G->n[0] - 1.0f/2) <  KILLRANGE/2.0f &&
                  FABS(P->i[1]/G->n[1] - 1.0f/2) <  KILLRANGE/2.0f &&
                       P->i[2]/G->n[2] - 1.0f/2  <  KILLRANGE/2.0f &&
                       P->i[2]                   >= 0);
      escaped = (P->i[2] < 0) && P->RI == 1;
      break;
    case 3:; // A semicolon is necessary here because C is a wonderful language :-)
      bool tooOld = P->time > 1000/C*SQRT(SQR(G->d[0]*G->n[0]) + SQR(G->d[1]*G->n[1]));  // Photon is "too old" if it has traveled more than 1000 times the diagonal xy box distance (may be in an infinite or near-infinite loop of cycling)
      P->alive = P->i[2] < G->n[2] && P->i[2] >= 0 && !tooOld;
      escaped = !P->alive && P->RI == 1 && !tooOld;
      bool photonTeleported = false;
      if(P->i[0] >= G->n[0]) {
        P->i[0] = 0; // Wrap x
        photonTeleported = true;
      }
      if(P->i[0] < 0) {
        P->i[0] = G->n[0]*(1-FLOATORDBLEPS); // Wrap x
        photonTeleported = true;
      }
      if(P->i[1] >= G->n[1]) {
        P->i[1] = 0; // Wrap y
        photonTeleported = true;
      }
      if(P->i[1] < 0) {
        P->i[1] = G->n[1]*(1-FLOATORDBLEPS); // Wrap y
        photonTeleported = true;
      }

      #ifdef __NVCC__ // If compiling for CUDA
      if(!threadIdx.x && !blockIdx.x)
      #elif defined(_OPENMP)
      #pragma omp master
      #endif
      {
        if(photonTeleported && (DC->evaluateCriteriaAtEndOfLife || depositionCriteriaMet(P,DC))) {
          updatePaths(P,Pa,G,DC,photonTeleported);
        }
      }
      break;
  }

  P->insideVolume = P->i[0] < G->n[0] && P->i[0] >= 0 &&
                    P->i[1] < G->n[1] && P->i[1] >= 0 &&
                    P->i[2] < G->n[2] && P->i[2] >= 0;
  
  if(escaped) {
    P->killed_escaped_collected = 1; // Escaped, may be overwritten by collected in formImage
    // We have to check formImage first because that's where we find out if the photon is collected
    if(O->image) formImage(P,G,LC,DC,O,nPhotonsCollectedPtr); // If image is not NULL then that's because useLightCollector was set to true (non-zero)
    if(O->FF && depositionCriteriaMet(P,DC)) formFarField(P,G,O);
  }
  if(!P->alive && G->boundaryType && depositionCriteriaMet(P,DC)) formEdgeFluxes(P,G,O);
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void getNewVoxelProperties(struct photon * const P, struct geometry const * const G, struct debug *D) {
  /* Get optical properties of current voxel.
   * If photon is outside cuboid, properties are those of
   * the closest defined voxel. */
  getNewj(G,P);
  P->mua    = G->muav   [G->M[P->j]];
  P->mus    = G->musv   [G->M[P->j]];
  P->g      = G->gv     [G->M[P->j]];
  P->CDFidx = G->CDFidxv[G->M[P->j]];
  // Refractive index is not retrieved here, but only when we refract in through a surface
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void getNormal(struct geometry const * const G, FLOATORDBL *nxPtr, FLOATORDBL *nyPtr, FLOATORDBL *nzPtr, long j) {
  float theta = G->interfaceNormals[2*j    ];
  float phi   = G->interfaceNormals[2*j + 1];
  *nxPtr = sin(theta)*cos(phi); // Normal vector x composant
  *nyPtr = sin(theta)*sin(phi); // Normal vector y composant
  *nzPtr = cos(theta); // Normal vector z composant
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void getInterpolatedNormal(struct geometry const * const G, struct photon *P, FLOATORDBL *nxPtr, FLOATORDBL *nyPtr, FLOATORDBL *nzPtr, long j) {
  long ix_coerced = (long)((P->i[0] < 0)? 0: ((P->i[0] >= G->n[0])? G->n[0]-1: P->i[0]));
  long iy_coerced = (long)((P->i[1] < 0)? 0: ((P->i[1] >= G->n[1])? G->n[1]-1: P->i[1]));
  long iz_coerced = (long)((P->i[2] < 0)? 0: ((P->i[2] >= G->n[2])? G->n[2]-1: P->i[2]));

  // Determine which set of 8 voxels to interpolate between. Let the corner with lowest indices have indices ix,iy,iz and the corner with highest indices have indices ix+1,iy+1,iz+1.
  long ix = (long)FLOOR(ix_coerced - 0.5);
  long iy = (long)FLOOR(iy_coerced - 0.5);
  long iz = (long)FLOOR(iz_coerced - 0.5);
  
  // Calculate weights based on where in the 8-voxel space the photon is
  FLOATORDBL wx = (FLOATORDBL)(ix_coerced - 0.5 - ix);
  FLOATORDBL wy = (FLOATORDBL)(iy_coerced - 0.5 - iy);
  FLOATORDBL wz = (FLOATORDBL)(iz_coerced - 0.5 - iz);
  
  // Add up contributions from those of the 8 voxels that are of the same refractive index
  *nxPtr = *nyPtr = *nzPtr = 0;
  FLOATORDBL nx,ny,nz;
  if(ix >= 0          && iy >= 0          && iz >= 0         ) {
    long jcorner = ix      +  iy     *G->n[0] +  iz     *G->n[0]*G->n[1];
    if(G->RIv[G->M[jcorner]] == G->RIv[G->M[j]]) {
      getNormal(G,&nx,&ny,&nz,jcorner);
      *nxPtr += (1-wx)*(1-wy)*(1-wz)*nx;
      *nyPtr += (1-wx)*(1-wy)*(1-wz)*ny;
      *nzPtr += (1-wx)*(1-wy)*(1-wz)*nz;
    }
  }
  if(ix >= 0          && iy >= 0          && iz < G->n[2] - 1) {
    long jcorner = ix      +  iy     *G->n[0] + (iz + 1)*G->n[0]*G->n[1];
    if(G->RIv[G->M[jcorner]] == G->RIv[G->M[j]]) {
      getNormal(G,&nx,&ny,&nz,jcorner);
      *nxPtr += (1-wx)*(1-wy)*wz*nx;
      *nyPtr += (1-wx)*(1-wy)*wz*ny;
      *nzPtr += (1-wx)*(1-wy)*wz*nz;
    }
  }
  if(ix >= 0          && iy < G->n[1] - 1 && iz >= 0         ) {
    long jcorner = ix      + (iy + 1)*G->n[0] +  iz     *G->n[0]*G->n[1];
    if(G->RIv[G->M[jcorner]] == G->RIv[G->M[j]]) {
      getNormal(G,&nx,&ny,&nz,jcorner);
      *nxPtr += (1-wx)*wy*(1-wz)*nx;
      *nyPtr += (1-wx)*wy*(1-wz)*ny;
      *nzPtr += (1-wx)*wy*(1-wz)*nz;
    }
  }
  if(ix >= 0          && iy < G->n[1] - 1 && iz < G->n[2] - 1) {
    long jcorner = ix      + (iy + 1)*G->n[0] + (iz + 1)*G->n[0]*G->n[1];
    if(G->RIv[G->M[jcorner]] == G->RIv[G->M[j]]) {
      getNormal(G,&nx,&ny,&nz,jcorner);
      *nxPtr += (1-wx)*wy*wz*nx;
      *nyPtr += (1-wx)*wy*wz*ny;
      *nzPtr += (1-wx)*wy*wz*nz;
    }
  }
  if(ix < G->n[0] - 1 && iy >= 0          && iz >= 0         ) {
    long jcorner = (ix + 1) +  iy     *G->n[0] +  iz     *G->n[0]*G->n[1];
    if(G->RIv[G->M[jcorner]] == G->RIv[G->M[j]]) {
      getNormal(G,&nx,&ny,&nz,jcorner);
      *nxPtr += wx*(1-wy)*(1-wz)*nx;
      *nyPtr += wx*(1-wy)*(1-wz)*ny;
      *nzPtr += wx*(1-wy)*(1-wz)*nz;
    }
  }
  if(ix < G->n[0] - 1 && iy >= 0          && iz < G->n[2] - 1) {
    long jcorner = (ix + 1) +  iy     *G->n[0] + (iz + 1)*G->n[0]*G->n[1];
    if(G->RIv[G->M[jcorner]] == G->RIv[G->M[j]]) {
      getNormal(G,&nx,&ny,&nz,jcorner);
      *nxPtr += wx*(1-wy)*wz*nx;
      *nyPtr += wx*(1-wy)*wz*ny;
      *nzPtr += wx*(1-wy)*wz*nz;
    }
  }
  if(ix < G->n[0] - 1 && iy < G->n[1] - 1 && iz >= 0         ) {
    long jcorner = (ix + 1) + (iy + 1)*G->n[0] +  iz     *G->n[0]*G->n[1];
    if(G->RIv[G->M[jcorner]] == G->RIv[G->M[j]]) {
      getNormal(G,&nx,&ny,&nz,jcorner);
      *nxPtr += wx*wy*(1-wz)*nx;
      *nyPtr += wx*wy*(1-wz)*ny;
      *nzPtr += wx*wy*(1-wz)*nz;
    }
  }
  if(ix < G->n[0] - 1 && iy < G->n[1] - 1 && iz < G->n[2] - 1) {
    long jcorner = (ix + 1) + (iy + 1)*G->n[0] + (iz + 1)*G->n[0]*G->n[1];
    if(G->RIv[G->M[jcorner]] == G->RIv[G->M[j]]) {
      getNormal(G,&nx,&ny,&nz,jcorner);
      *nxPtr += wx*wy*wz*nx;
      *nyPtr += wx*wy*wz*ny;
      *nzPtr += wx*wy*wz*nz;
    }
  }
  FLOATORDBL norm = SQRT(SQR(*nxPtr) + SQR(*nyPtr) + SQR(*nzPtr));
  *nxPtr /= norm;
  *nyPtr /= norm;
  *nzPtr /= norm;
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void propagatePhoton(struct photon * const P, struct geometry const * const G, struct outputs const *O, struct depositionCriteria *DC, struct paths * const Pa, struct debug *D) {
  long idx;

  P->sameVoxel = true;
  
  FLOATORDBL s = min(P->stepLeft/P->mus,min(P->D[0],min(P->D[1],P->D[2])));

  P->stepLeft  = s==P->stepLeft/P->mus? 0: P->stepLeft - s*P->mus; // zero case is to avoid rounding errors
  P->time     += s*P->RI/C;
  
  for(idx=0;idx<3;idx++) { // Propagate photon
    long i_old = (long)FLOOR(P->i[idx]);
    if(s == P->D[idx]) { // If we're supposed to go to the voxel boundary along this dimension
      P->i[idx] = (P->u[idx] > 0)? i_old + 1: i_old - FLOATORDBLEPS*(labs(i_old)+1);
      P->D[idx] = G->d[idx]/FABS(P->u[idx]); // Reset voxel boundary distance
      P->sameVoxel = false;
    } else { // We're supposed to remain in the same voxel along this dimension
      P->i[idx] += s*P->u[idx]/G->d[idx]; // First take the expected step (including various rounding errors)
      if(FLOOR(P->i[idx]) != i_old) P->i[idx] = (P->u[idx] > 0)? i_old + 1 - FLOATORDBLEPS*(labs(i_old)+1): i_old ; // If photon due to rounding errors actually crossed the border, set it to be barely in the original voxel
      P->D[idx] -= s;
    }
  }

  if(!P->sameVoxel) {
    long j_new = ((P->i[2] < 0)? 0: ((P->i[2] >= G->n[2])? G->n[2]-1: (long)FLOOR(P->i[2])))*G->n[0]*G->n[1] +
                 ((P->i[1] < 0)? 0: ((P->i[1] >= G->n[1])? G->n[1]-1: (long)FLOOR(P->i[1])))*G->n[0]         +
                 ((P->i[0] < 0)? 0: ((P->i[0] >= G->n[0])? G->n[0]-1: (long)FLOOR(P->i[0]))); // Index values are restrained to integers in the interval [0,n-1]
    // https://physics.stackexchange.com/questions/435512/snells-law-in-vector-form  or  http://www.starkeffects.com/snells-law-vector.shtml#:~:text=Snell's%20Law%20in%20Vector%20Form&text=Since%20the%20incident%20ray%2C%20the,yourself%20to%20fit%20the%20equation.
    FLOATORDBL mu = P->RI/G->RIv[G->M[j_new]]; // RI ratio
    if(mu != 1) { // If there's a refractive index change
      bool photonReflected = false;
      FLOATORDBL nx,ny,nz;
      getInterpolatedNormal(G,P,&nx,&ny,&nz,j_new);
      FLOATORDBL cos_in = nx*P->u[0] + ny*P->u[1] + nz*P->u[2]; // dot product of n and u
      if(cos_in > 0) { // If cos_in is negative, we're dealing with a case in which the photon is headed in a direction that, according to the surface normal, is actually away from the medium. Then we choose not to do any reflection or refraction.
        FLOATORDBL cos_out_sqr = 1 - SQR(mu)*(1 - SQR(cos_in));
        if(cos_out_sqr > 0) { // If we don't experience total internal reflection
          FLOATORDBL cos_out = SQRT(cos_out_sqr);
          FLOATORDBL R = SQR((mu*cos_in  - cos_out)/(mu*cos_in  + cos_out))/2 +
                         SQR((mu*cos_out - cos_in )/(mu*cos_out + cos_in ))/2; // R is the reflectivity assuming equal probability of p or s polarization (unpolarized light at all times)
          photonReflected = RandomNum <= R;
        } else photonReflected = true;
        if(photonReflected) { // u_refl = u - 2*n*(u dot n) = u - 2*n*cos(theta_in)
          P->u[0] -= 2*nx*cos_in;
          P->u[1] -= 2*ny*cos_in;
          P->u[2] -= 2*nz*cos_in;
          P->D[0] = P->u[0]? (FLOOR(P->i[0]) + (P->u[0]>0) - P->i[0])*G->d[0]/P->u[0]: INFINITY; // Recalculate voxel boundary distance for x
          P->D[1] = P->u[1]? (FLOOR(P->i[1]) + (P->u[1]>0) - P->i[1])*G->d[1]/P->u[1]: INFINITY; // Recalculate voxel boundary distance for y
          P->D[2] = P->u[2]? (FLOOR(P->i[2]) + (P->u[2]>0) - P->i[2])*G->d[2]/P->u[2]: INFINITY; // Recalculate voxel boundary distance for z
          // We deliberately do not get the refractive index of the new voxel here, since a reflection means the photon is effectively still in the same medium, despite perhaps temporarily traveling in the new medium's voxel
          if(G->M[j_new] >= DC->minIdx && G->M[j_new] <= DC->maxIdx)
            P->reflections++;
        } else { // u_refr = sqrt(1 - mu^2*(1 - (u dot n)^2))*n + mu*(u - (u dot n)*n) = sqrt(1 - mu^2*(1 - cos(theta_in)^2))*n + mu*(u - cos(theta_in)*n) = (sqrt(1 - mu^2*(1 - cos(theta_in)^2)) - mu*cos(theta_in))*n + mu*u
          FLOATORDBL ncoeff = SQRT(cos_out_sqr) - mu*cos_in;
          P->u[0] = ncoeff*nx + mu*P->u[0];
          P->u[1] = ncoeff*ny + mu*P->u[1];
          P->u[2] = ncoeff*nz + mu*P->u[2];
          P->D[0] = P->u[0]? (FLOOR(P->i[0]) + (P->u[0]>0) - P->i[0])*G->d[0]/P->u[0]: INFINITY; // Recalculate voxel boundary distance for x
          P->D[1] = P->u[1]? (FLOOR(P->i[1]) + (P->u[1]>0) - P->i[1])*G->d[1]/P->u[1]: INFINITY; // Recalculate voxel boundary distance for y
          P->D[2] = P->u[2]? (FLOOR(P->i[2]) + (P->u[2]>0) - P->i[2])*G->d[2]/P->u[2]: INFINITY; // Recalculate voxel boundary distance for z
          P->RI = G->RIv[G->M[j_new]]; // Since we have refracted into the new medium, we retrieve the new refractive index
          if(G->M[j_new] >= DC->minIdx && G->M[j_new] <= DC->maxIdx) {
            P->refractions++;
            P->interfaceTransitions++;
          }
        }
      }
    } else if(G->M[j_new] != G->M[P->j] && G->RIv[G->M[P->j]] == P->RI) { // No refraction or reflection, but we are in a new medium. The G->RIv[G->M[P->j]] == P->RI check is to ensure that we are not just coming out from a reflection event, because in that case the RI of the old voxel will not correspond to the P->RI value.
      if(G->M[j_new] >= DC->minIdx && G->M[j_new] <= DC->maxIdx)
        P->interfaceTransitions++;
    }
  }
  
  FLOATORDBL absorb = -P->weight*EXPM1(-P->mua*s);   // photon weight absorbed at this step. expm1(x) = exp(x) - 1, accurate even for very small x 
  P->weight -= absorb;             // decrement WEIGHT by amount absorbed

  if(P->insideVolume) {  // only save data if the photon is inside simulation cuboid
    if(!DC->evaluateCriteriaAtEndOfLife) {
      if(O->NFR && depositionCriteriaMet(P,DC)) {
        atomicAddWrapper(&O->NFR[P->j],absorb);
      }
    } else { // store indices and weights in pseudosparse array, to later add to NFR if photon ends up on the light collector
      if(P->recordElems == P->recordSize) {
        P->recordSize *= 2; // double the record's size
        P->j_record = (long *)reallocWrapper(P->j_record,P->recordSize/2*sizeof(long),P->recordSize*sizeof(long));
        P->weight_record = (FLOATORDBL *)reallocWrapper(P->weight_record,P->recordSize/2*sizeof(FLOATORDBL),P->recordSize*sizeof(FLOATORDBL));
      }
      P->j_record[P->recordElems] = P->j;
      P->weight_record[P->recordElems] = absorb;
      P->recordElems++;
    }
  }
  #ifdef __NVCC__ // If compiling for CUDA
  if(!threadIdx.x && !blockIdx.x)
  #elif defined(_OPENMP)
  #pragma omp master
  #endif
  {
    updatePaths(P,Pa,G,DC,false);
  }
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void checkRoulette(struct photon * const P) {
  /**** CHECK ROULETTE
   * If photon weight below THRESHOLD, then terminate photon using Roulette technique.
   * Photon has CHANCE probability of having its weight increased by factor of 1/CHANCE,
   * and 1-CHANCE probability of terminating. */
  if(P->weight < THRESHOLD) {
    if(RandomNum <= CHANCE) P->weight /= CHANCE;
    else P->alive = false;
  }
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void scatterPhoton(struct photon * const P, struct geometry const * const G, struct paths *Pa, struct depositionCriteria *DC, struct debug *D) {
  FLOATORDBL costheta;
  if(ISNAN(P->g)) {
    // Sample for theta using the cumulative distribution function (CDF)
    long jTheta = binaryTreeSearch(RandomNum,CDFSIZE,G->CDFs + P->CDFidx*CDFSIZE);
    costheta = COS((jTheta + RandomNum)*PI/(CDFSIZE-1));
  } else {
    // Sample for costheta using Henyey-Greenstein scattering
    costheta = FABS(P->g) == 1.0? P->g:
                          FABS(P->g) <= SQRT(FLOATORDBLEPS)? 2*RandomNum - 1:
                          (1 + P->g*P->g - SQR((1 - P->g*P->g)/(1 - P->g + 2*P->g*RandomNum)))/(2*P->g);
  }

  FLOATORDBL sintheta = SQRT(1 - costheta*costheta);
  FLOATORDBL phi = 2*PI*RandomNum;
  FLOATORDBL cosphi = COS(phi);
  FLOATORDBL sinphi = SIN(phi);
  
  if(FABS(P->u[2]) < 1) {
    FLOATORDBL ux_temp =  sintheta*(P->u[0]*P->u[2]*cosphi - P->u[1]*sinphi)/SQRT(P->u[0]*P->u[0] + P->u[1]*P->u[1]) + P->u[0]*costheta;
    FLOATORDBL uy_temp =  sintheta*(P->u[1]*P->u[2]*cosphi + P->u[0]*sinphi)/SQRT(P->u[0]*P->u[0] + P->u[1]*P->u[1]) + P->u[1]*costheta;
    P->u[2]            = -sintheta*(                cosphi                 )*SQRT(P->u[0]*P->u[0] + P->u[1]*P->u[1]) + P->u[2]*costheta;
    P->u[1]            = uy_temp;
    P->u[0]            = ux_temp;
  } else {
    P->u[0] = sintheta*cosphi;
    P->u[1] = sintheta*sinphi;
    P->u[2] = costheta*SIGN(P->u[2]);
  }

  // Calculate distances to next voxel boundary planes
  for(long idx=0;idx<3;idx++) P->D[idx] = P->u[idx]? (FLOOR(P->i[idx]) + (P->u[idx]>0) - P->i[idx])*G->d[idx]/P->u[idx] : INFINITY;
  
  P->stepLeft  = -LOG(RandomNum);

  if(G->M[P->j] >= DC->minIdx && G->M[P->j] <= DC->maxIdx)
    P->scatterings++;

  #ifdef __NVCC__ // If compiling for CUDA
  if(!threadIdx.x && !blockIdx.x)
  #elif defined(_OPENMP)
  #pragma omp master
  #endif
  {
    updatePaths(P,Pa,G,DC,false);
  }
}

void normalizeDepositionAndResetO(struct source const * const B, struct geometry const * const G, struct lightCollector const * const LC,
        struct outputs *O, struct MATLABoutputs *O_MATLAB, long iWavelength, double Pfraction) {
  long j;
  double V = G->d[0]*G->d[1]*G->d[2]; // Voxel volume
  long L = G->n[0]*G->n[1]*G->n[2]; // Total number of voxels in cuboid
  long L_LC = LC->res[0]*LC->res[0]; // Total number of spatial pixels in light collector planes
  long L_FF = G->farFieldRes*G->farFieldRes; // Total number of pixels in the far field array
  // Normalize deposition to yield normalized fluence rate (NFR). For fluorescence, the result is relative to
  // the incident excitation power (not emitted fluorescence power).
  double normfactor = (double)O->nPhotons/Pfraction;
  if(B->S) { // For a 3D source distribution (e.g., fluorescence)
    normfactor /= B->power;
  } else if(B->beamType == 2 && (G->boundaryType == 0 || G->boundaryType == 2)) { // For infinite plane wave launched into volume without absorbing walls
    normfactor /= KILLRANGE*KILLRANGE;
  }
  
  O->nPhotons = 0;
  O->nPhotonsCollected = 0;
  if(O->NFR) for(j=0;j<L   ;j++) {
    O_MATLAB->NFR[j + (unsigned long long)iWavelength*G->n[0]*G->n[1]*G->n[2]] = (float)(O->NFR[j]/(V*normfactor*G->muav[G->M[j]]));
    O->NFR[j] = 0;
  }
  if(O->FF) for(j=0;j<L_FF;j++) {
    O_MATLAB->FF[j + iWavelength*G->farFieldRes*G->farFieldRes] = (float)(O->FF[j]/normfactor);
    O->FF[j] = 0;
  }
  if(O->image) {
    if(L_LC > 1) for(j=0;j<L_LC*LC->res[1];j++) {
      O_MATLAB->image[j + iWavelength*LC->res[0]*LC->res[0]*LC->res[1]] = (float)(O->image[j]/(LC->FSorNA*LC->FSorNA/L_LC*normfactor));
      O->image[j] = 0;
    }
    else         for(j=0;j<     LC->res[1];j++) {
      O_MATLAB->image[j + iWavelength*LC->res[0]*LC->res[0]*LC->res[1]] = (float)(O->image[j]/normfactor);
      O->image[j] = 0;
    }
  }
  if(G->boundaryType == 1) {
    for(j=0;j<G->n[1]*G->n[2];j++) {
      O_MATLAB->NI_xpos[j + iWavelength*G->n[1]*G->n[2]] = (float)(O->NI_xpos[j]/(G->d[1]*G->d[2]*normfactor));
      O->NI_xpos[j] = 0;
      O_MATLAB->NI_xneg[j + iWavelength*G->n[1]*G->n[2]] = (float)(O->NI_xneg[j]/(G->d[1]*G->d[2]*normfactor));
      O->NI_xneg[j] = 0;
    }
    for(j=0;j<G->n[0]*G->n[2];j++) {
      O_MATLAB->NI_ypos[j + iWavelength*G->n[0]*G->n[2]] = (float)(O->NI_ypos[j]/(G->d[0]*G->d[2]*normfactor));
      O->NI_ypos[j] = 0;
      O_MATLAB->NI_yneg[j + iWavelength*G->n[0]*G->n[2]] = (float)(O->NI_yneg[j]/(G->d[0]*G->d[2]*normfactor));
      O->NI_yneg[j] = 0;
    }
    for(j=0;j<G->n[0]*G->n[1];j++) {
      O_MATLAB->NI_zpos[j + iWavelength*G->n[0]*G->n[1]] = (float)(O->NI_zpos[j]/(G->d[0]*G->d[1]*normfactor));
      O->NI_zpos[j] = 0;
      O_MATLAB->NI_zneg[j + iWavelength*G->n[0]*G->n[1]] = (float)(O->NI_zneg[j]/(G->d[0]*G->d[1]*normfactor));
      O->NI_zneg[j] = 0;
    }
  } else if(G->boundaryType == 2) for(j=0;j<KILLRANGE*G->n[0]*KILLRANGE*G->n[1];j++) {
    O_MATLAB->NI_zneg[j + iWavelength*G->n[0]*G->n[1]*(G->boundaryType == 2? KILLRANGE*KILLRANGE: 1)] = (float)(O->NI_zneg[j]/(G->d[0]*G->d[1]*normfactor));
    O->NI_zneg[j] = 0;
  } else if(G->boundaryType == 3) for(j=0;j<G->n[0]*G->n[1];j++) {
    O_MATLAB->NI_zpos[j + iWavelength*G->n[0]*G->n[1]] = (float)(O->NI_zpos[j]/(G->d[0]*G->d[1]*normfactor));
    O->NI_zpos[j] = 0;
    O_MATLAB->NI_zneg[j + iWavelength*G->n[0]*G->n[1]] = (float)(O->NI_zneg[j]/(G->d[0]*G->d[1]*normfactor));
    O->NI_zneg[j] = 0;
  }
}
