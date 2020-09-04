
struct debug {
  double             dbls[3];
  unsigned long long ulls[3];
};

struct geometry { // Struct type for the constant geometry definitions, including the wavelength-dependent optical properties and the boundary type
  FLOATORDBL     d[3];
  long           n[3];
  long           farFieldRes;
  int            boundaryType;
  FLOATORDBL     *muav,*musv,*gv;
  unsigned char  *M;
  FLOATORDBL     *RIv;    // Refractive index of each z slice
};

struct beam { // Struct type for the constant beam definitions
  int            beamType;
  FLOATORDBL     *NFdist1; // Radial or X
  long           L_NF1;
  FLOATORDBL     NFwidth1;
  FLOATORDBL     *NFdist2; // Azimuthal or Y
  long           L_NF2;
  FLOATORDBL     NFwidth2;
  FLOATORDBL     *FFdist1; // Radial or X
  long           L_FF1;
  FLOATORDBL     FFwidth1;
  FLOATORDBL     *FFdist2; // Azimuthal or Y
  long           L_FF2;
  FLOATORDBL     FFwidth2;
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
  FLOATORDBL     i[3],u[3],D[3]; // Fractional position indices i, ray trajectory unit vector u and distances D to next voxel boundary (yz, xz or yz) along current trajectory
  long           j; // Linear index of current voxel (or closest defined voxel if photon outside cuboid)
  FLOATORDBL     mua,mus,g,RI; // Absorption, scattering, anisotropy and refractive index values at current photon position
  FLOATORDBL     stepLeft,weight,time;
  bool           insideVolume,alive,sameVoxel;
  PRNG_t         PRNGstate; // "State" of the Mersenne Twister pseudo-random number generator
  long           recordSize; // Current size of the list of voxels in which power has been deposited, used only if calcNFRdet is true
  long           recordElems; // Current number of elements used of the record. Starts at 0 every photon launch, used only if calcNFRdet is true
  long           *j_record; // List of the indices of the voxels in which the current photon has deposited power, used only if calcNFRdet is true
  FLOATORDBL     *weight_record; // List of the weights that have been deposited into the voxels, used only if calcNFRdet is true
};

struct paths { // Struct type for storing the paths taken by the nExamplePaths first photons simulated by the master thread
  long           pathsElems; // Current number of elements used
  long           threadNum;
  long           nExamplePaths; // Number of photons to store the path of
  long           nMasterPhotonsLaunched; // Number of photons the master thread has launched
  long           pathsSize; // Current size of the list
  FLOATORDBL     *data; // Array containing x, y, z and weight data for the photons, in which the paths of different photons are separated by four NaNs
};

struct outputs {
  unsigned long long nPhotons;
  double * NFR;
  double * NFRdet;
  double * image;
  double * FF;
  double * NI_xpos;
  double * NI_xneg;
  double * NI_ypos;
  double * NI_yneg;
  double * NI_zpos;
  double * NI_zneg;
};

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
double sqr(double x) {return x*x;}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
float sqrf(float x) {return x*x;}

#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ < 600
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
void atomicAddWrapperULL(unsigned long long *ptr, unsigned long long val) {
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
void atomicAddWrapperDBL(double *ptr, double val) {
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
void createDeviceStructs(struct geometry const *G, struct geometry **G_devptr,
                         struct beam const *B, struct beam **B_devptr,
                         struct lightCollector const *LC, struct lightCollector **LC_devptr,
                         struct paths *Pa, struct paths **Pa_devptr,
                         struct outputs *O, struct outputs **O_devptr,
                         int nM, long L,
                         struct debug *D, struct debug **D_devptr) {
  // Allocate and copy geometry struct, including smallArrays
  struct geometry G_tempvar = *G;
  long size_smallArrays = (nM*3 + G->n[2] + B->L_NF1 + B->L_FF1 + B->L_NF2 + B->L_FF2)*sizeof(FLOATORDBL);
  gpuErrchk(cudaMalloc(&G_tempvar.muav,size_smallArrays)); // This is to allocate all of smallArrays
  gpuErrchk(cudaMemcpy(G_tempvar.muav,G->muav,size_smallArrays,cudaMemcpyHostToDevice)); // And for copying all of it to global device memory
  gpuErrchk(cudaMalloc(&G_tempvar.M   ,       L*sizeof(unsigned char)));
  gpuErrchk(cudaMemcpy(G_tempvar.M   , G->M   ,       L*sizeof(unsigned char),cudaMemcpyHostToDevice));

  gpuErrchk(cudaMalloc(G_devptr, sizeof(struct geometry)));
  gpuErrchk(cudaMemcpy(*G_devptr,&G_tempvar,sizeof(struct geometry),cudaMemcpyHostToDevice));

  // Allocate and copy beam struct
  struct beam B_tempvar = *B;
  if(B->S) {
    gpuErrchk(cudaMalloc(&B_tempvar.S, (L+1)*sizeof(FLOATORDBL)));
    gpuErrchk(cudaMemcpy(B_tempvar.S, B->S, (L+1)*sizeof(FLOATORDBL),cudaMemcpyHostToDevice));
  }
  gpuErrchk(cudaMalloc(B_devptr, sizeof(struct beam)));
  gpuErrchk(cudaMemcpy(*B_devptr,&B_tempvar,sizeof(struct beam),cudaMemcpyHostToDevice));

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
  if(O->NFRdet) {
    gpuErrchk(cudaMalloc(&O_tempvar.NFRdet, L*sizeof(double)));
    gpuErrchk(cudaMemset(O_tempvar.NFRdet,0,L*sizeof(double)));
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
    gpuErrchk(cudaMalloc(&O_tempvar.NI_xneg, G->n[1]*G->n[2]*sizeof(double)));
    gpuErrchk(cudaMemset(O_tempvar.NI_xneg,0,G->n[1]*G->n[2]*sizeof(double)));
    gpuErrchk(cudaMalloc(&O_tempvar.NI_ypos, G->n[0]*G->n[2]*sizeof(double)));
    gpuErrchk(cudaMemset(O_tempvar.NI_ypos,0,G->n[0]*G->n[2]*sizeof(double)));
    gpuErrchk(cudaMalloc(&O_tempvar.NI_yneg, G->n[0]*G->n[2]*sizeof(double)));
    gpuErrchk(cudaMemset(O_tempvar.NI_yneg,0,G->n[0]*G->n[2]*sizeof(double)));
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
                                  struct beam const *B, struct beam *B_dev,
                                  struct lightCollector const *LC, struct lightCollector *LC_dev,
                                  struct paths *Pa, struct paths *Pa_dev,
                                  struct outputs *O, struct outputs *O_dev,
                                  long L,
                                  struct debug *D, struct debug *D_dev) {
  struct geometry G_temp; gpuErrchk(cudaMemcpy(&G_temp, G_dev, sizeof(struct geometry),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaFree(G_temp.muav)); // This frees all of smallArrays in the global memory on the device
  gpuErrchk(cudaFree(G_temp.M));
  gpuErrchk(cudaFree(G_dev));

  struct beam B_temp; gpuErrchk(cudaMemcpy(&B_temp, B_dev, sizeof(struct beam),cudaMemcpyDeviceToHost));
  gpuErrchk(cudaFree(B_temp.S));
  gpuErrchk(cudaFree(B_dev));
  
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
  if(O->NFR) {
    gpuErrchk(cudaMemcpy(O->NFR, O_temp.NFR, L*sizeof(double),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.NFR));
  }
  if(O->NFRdet) {
    gpuErrchk(cudaMemcpy(O->NFRdet, O_temp.NFRdet, L*sizeof(double),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.NFRdet));
  }
  if(O->image) {
    gpuErrchk(cudaMemcpy(O->image, O_temp.image, LC->res[0]*LC->res[0]*LC->res[1]*sizeof(double),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.image));
  }
  if(O->FF) {
    gpuErrchk(cudaMemcpy(O->FF, O_temp.FF, G->farFieldRes*G->farFieldRes*sizeof(double),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.FF));
  }
  if(O->NI_xpos) {
    gpuErrchk(cudaMemcpy(O->NI_xpos, O_temp.NI_xpos, G->n[1]*G->n[2]*sizeof(double),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.NI_xpos));
    gpuErrchk(cudaMemcpy(O->NI_xneg, O_temp.NI_xneg, G->n[1]*G->n[2]*sizeof(double),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.NI_xneg));
    gpuErrchk(cudaMemcpy(O->NI_ypos, O_temp.NI_ypos, G->n[0]*G->n[2]*sizeof(double),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.NI_ypos));
    gpuErrchk(cudaMemcpy(O->NI_yneg, O_temp.NI_yneg, G->n[0]*G->n[2]*sizeof(double),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.NI_yneg));
    gpuErrchk(cudaMemcpy(O->NI_zpos, O_temp.NI_zpos, G->n[0]*G->n[1]*sizeof(double),cudaMemcpyDeviceToHost));
    gpuErrchk(cudaFree(O_temp.NI_zpos));
  }
  if(O->NI_zneg) {
    gpuErrchk(cudaMemcpy(O->NI_zneg, O_temp.NI_zneg, (G->boundaryType == 2? KILLRANGE*KILLRANGE: 1)*G->n[0]*G->n[1]*sizeof(double),cudaMemcpyDeviceToHost));
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
void launchPhoton(struct photon * const P, struct beam const * const B, struct geometry const * const G, struct paths * const Pa) {
  FLOATORDBL X,Y,r,phi,tanphiX,tanphiY,costheta,sintheta;
  long   j,idx;
  FLOATORDBL target[3]={0},w0[3];
  
  P->sameVoxel = false;
  P->weight = 1;
  P->recordElems = 0;
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
      P->time = 0;
    } else switch (B->beamType) {
      case 0: // pencil beam
        P->i[0] = (B->focus[0] - B->focus[2]*B->u[0]/B->u[2])/G->d[0] + G->n[0]/2.0f;
        P->i[1] = (B->focus[1] - B->focus[2]*B->u[1]/B->u[2])/G->d[1] + G->n[1]/2.0f;
        P->i[2] = 0;
        for(idx=0;idx<3;idx++) P->u[idx] = B->u[idx];
        P->time = -G->RIv[0]/C*SQRT(SQR((P->i[0] - G->n[0]/2.0f)*G->d[0] - B->focus[0]) +
                                    SQR((P->i[1] - G->n[1]/2.0f)*G->d[1] - B->focus[1]) +
                                    SQR((P->i[2]               )*G->d[2] - B->focus[2])); // Starting time is set so that the wave crosses the focal plane at time = 0
        break;
      case 1: // isotropically emitting point source
        P->i[0] = B->focus[0]/G->d[0] + G->n[0]/2.0f;
        P->i[1] = B->focus[1]/G->d[1] + G->n[1]/2.0f;
        P->i[2] = B->focus[2]/G->d[2];
        costheta = 1 - 2*RandomNum;
        sintheta = SQRT(1 - costheta*costheta);
        phi = 2*PI*RandomNum;
        P->u[0] = sintheta*COS(phi);
        P->u[1] = sintheta*SIN(phi);
        P->u[2] = costheta;
        P->time = 0;
        break;
      case 2: // infinite plane wave
        P->i[0] = ((G->boundaryType==1)? 1: KILLRANGE)*G->n[0]*(RandomNum-0.5f) + G->n[0]/2.0f; // Generates a random ix coordinate within the cuboid
        P->i[1] = ((G->boundaryType==1)? 1: KILLRANGE)*G->n[1]*(RandomNum-0.5f) + G->n[1]/2.0f; // Generates a random iy coordinate within the cuboid
        P->i[2] = 0;
        for(idx=0;idx<3;idx++) P->u[idx] = B->u[idx];
        P->time = G->RIv[0]/C*((P->i[0] - G->n[0]/2.0f)*G->d[0]*B->u[0] +
                               (P->i[1] - G->n[1]/2.0f)*G->d[1]*B->u[1] +
                               (P->i[2]               )*G->d[2]*B->u[2]); // Starting time is set so that the wave crosses (x=0,y=0,z=0) at time = 0
        break;
      case 3: // Laguerre-Gaussian LG01 beam
        phi     = RandomNum*2*PI;
        axisrotate(B->v,B->u,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
        r    = B->NFwidth1*SQRT(((FLOATORDBL)gsl_sf_lambert_Wm1(-RandomNum*EXP(-1.0f))+1)/(-2))/1.50087f; // for target calculation
        for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
        phi     = B->FFwidth1*SQRT(((FLOATORDBL)gsl_sf_lambert_Wm1(-RandomNum*EXP(-1.0f))+1)/(-2))/1.50087f; // for trajectory calculation. The sqrt is valid within paraxial approximation.
        axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
        P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0f; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
        P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0f;
        P->i[2] = 0;
        P->time = -G->RIv[0]/C*SQRT(SQR((P->i[0] - G->n[0]/2.0f)*G->d[0] - target[0]) +
                                    SQR((P->i[1] - G->n[1]/2.0f)*G->d[1] - target[1]) +
                                    SQR((P->i[2]               )*G->d[2] - target[2])); // Starting time is set so that the wave crosses the focal plane at time = 0
        break;
      case 4: // Radial
        // Near Field
        axisrotate(B->v,B->u,RandomNum*2*PI,w0); // w0 unit vector now points in the direction from focus center point to ray target point
        
        if(*B->NFdist1 == -1) { // Top-hat radial distribution
          r = B->NFwidth1*SQRT(RandomNum); // for target calculation
        } else if(*B->NFdist1 == -2) { // Gaussian radial distribution
          r = B->NFwidth1*SQRT(-0.5f*LOG(RandomNum)); // for target calculation
        } else { // Custom distribution
          r = B->NFwidth1*(binaryTreeSearch(RandomNum,B->L_NF1-1,B->NFdist1)+RandomNum)/(B->L_NF1-1);
        }
        for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
        
        // Far Field
        axisrotate(B->v,B->u,RandomNum*2*PI,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.

        if(*B->FFdist1 == -1) { // Top-hat radial distribution
          phi = ATAN(TAN(B->FFwidth1)*SQRT(RandomNum)); // for trajectory calculation. The sqrt is valid within paraxial approximation.
        } else if(*B->FFdist1 == -2) { // Gaussian radial distribution
          phi = ATAN(TAN(B->FFwidth1)*SQRT(-0.5f*LOG(RandomNum))); // for trajectory calculation. The sqrt is valid within paraxial approximation.
        } else if(*B->FFdist1 == -3) { // Lambertian
          phi = ASIN(SQRT(RandomNum));
        } else { // Custom distribution
          phi = ATAN(TAN(B->FFwidth1)*(binaryTreeSearch(RandomNum,B->L_FF1-1,B->FFdist1)+RandomNum)/(B->L_FF1-1));
        }
        axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
        
        P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0f; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
        P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0f;
        P->i[2] = 0;
        P->time = -G->RIv[0]/C*SQRT(SQR((P->i[0] - G->n[0]/2.0f)*G->d[0] - target[0]) +
                                    SQR((P->i[1] - G->n[1]/2.0f)*G->d[1] - target[1]) +
                                    SQR((P->i[2]               )*G->d[2] - target[2])); // Starting time is set so that the wave crosses the focal plane at time = 0
        break;
      case 5: // X/Y
        // Near Field
        if(*B->NFdist1 == -1) { // Top-hat X distribution
          X = B->NFwidth1*(RandomNum*2-1); // for target calculation
        } else if(*B->NFdist1 == -2) { // Gaussian X distribution
          X = B->NFwidth1*SQRT(-0.5f*LOG(RandomNum))*COS(2*PI*RandomNum); // Box-Muller transform, for target calculation
        } else { // Custom X distribution
          X = B->NFwidth1*((binaryTreeSearch(RandomNum,B->L_NF1-1,B->NFdist1)+RandomNum)/(B->L_NF1-1)*2-1);
        }
        if(*B->NFdist2 == -1) { // Top-hat Y distribution
          Y = B->NFwidth2*(RandomNum*2-1); // for target calculation
        } else if(*B->NFdist1 == -2) { // Gaussian Y distribution
          Y = B->NFwidth2*SQRT(-0.5f*LOG(RandomNum))*COS(2*PI*RandomNum); // Box-Muller transform, for target calculation
        } else { // Custom distribution
          Y = B->NFwidth2*((binaryTreeSearch(RandomNum,B->L_NF2-1,B->NFdist2)+RandomNum)/(B->L_NF2-1)*2-1);
        }
        for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + X*B->v[idx] + Y*B->w[idx];

        // Far Field
        if(*B->FFdist1 == -3) { // Lambertian
          axisrotate(B->v,B->u,RandomNum*2*PI,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction
          axisrotate(B->u,w0,ASIN(SQRT(RandomNum)),P->u); // ray propagation direction is found by rotating beam center axis around w0
        } else {
          if(*B->FFdist1 == -1) { // Top-hat phiX distribution
            tanphiX = TAN(B->FFwidth1)*(RandomNum*2-1); // for trajectory calculation
          } else if(*B->FFdist1 == -2) { // Gaussian phiX distribution
            tanphiX = TAN(B->FFwidth1)*SQRT(-0.5f*LOG(RandomNum))*COS(2*PI*RandomNum); // Box-Muller transform, for trajectory calculation
          } else { // Custom phi_X distribution
            tanphiX = TAN(B->FFwidth1)*((binaryTreeSearch(RandomNum,B->L_FF1,B->FFdist1)+RandomNum)/B->L_FF1*2-1);
          }
          if(*B->FFdist2 == -1) { // Top-hat phiY distribution
            tanphiY = TAN(B->FFwidth2)*(RandomNum*2-1); // for trajectory calculation
          } else if(*B->FFdist2 == -2) { // Gaussian phiY distribution
            tanphiY = TAN(B->FFwidth2)*SQRT(-0.5f*LOG(RandomNum))*COS(2*PI*RandomNum); // Box-Muller transform, for trajectory calculation
          } else { // Custom distribution
            tanphiY = TAN(B->FFwidth2)*((binaryTreeSearch(RandomNum,B->L_FF2,B->FFdist2)+RandomNum)/B->L_FF2*2-1);
          }
          axisrotate(B->v,B->u,ATAN2(tanphiX,tanphiY),w0); // w0 is now orthogonal to both beam propagation axis and ray propagation axis
          axisrotate(B->u,w0,ATAN(SQRT(tanphiX*tanphiX + tanphiY*tanphiY)),P->u); // ray propagation direction is found by rotating beam center axis around w0
        }
        P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0f; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
        P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0f;
        P->i[2] = 0;
        P->time = -G->RIv[0]/C*SQRT(SQR((P->i[0] - G->n[0]/2.0f)*G->d[0] - target[0]) +
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
    }
  } while(!P->alive); // If photon happened to be initialized outside the volume in which it is allowed to travel, we try again

  // Calculate distances to next voxel boundary planes
  for(idx=0;idx<3;idx++) P->D[idx] = P->u[idx]? (FLOOR(P->i[idx]) + (P->u[idx]>0) - P->i[idx])*G->d[idx]/P->u[idx] : INFINITY;
  
  P->stepLeft  = -LOG(RandomNum);
  
  #ifdef __NVCC__ // If compiling for CUDA
  if(!threadIdx.x && !blockIdx.x)
  #elif defined(_OPENMP)
  #pragma omp master
  #endif
  {
    Pa->nMasterPhotonsLaunched++;
    if(Pa->nMasterPhotonsLaunched <= Pa->nExamplePaths) {
      if(Pa->pathsElems + 2 > Pa->pathsSize) {
        #ifdef __NVCC__ // If compiling for CUDA
        Pa->pathsElems = -1; // Buffer not large enough, has to be resized by host
        #else
        Pa->pathsSize *= 2; // double the record's size
        Pa->data = (FLOATORDBL *)reallocWrapper(Pa->data,2*Pa->pathsSize*sizeof(FLOATORDBL),4*Pa->pathsSize*sizeof(FLOATORDBL));
        #endif
      }
      if(Pa->pathsElems != -1) {
        for(long i=0;i<4;i++) Pa->data[4*Pa->pathsElems+i] = NAN;
        Pa->pathsElems++;
        Pa->data[4*Pa->pathsElems  ] = (P->i[0] - G->n[0]/2.0f)*G->d[0]; // Store x
        Pa->data[4*Pa->pathsElems+1] = (P->i[1] - G->n[1]/2.0f)*G->d[1]; // Store y
        Pa->data[4*Pa->pathsElems+2] = (P->i[2]               )*G->d[2]; // Store z
        Pa->data[4*Pa->pathsElems+3] = P->weight;                       // Store photon weight
        Pa->pathsElems++;
      }
    }
  }
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void formImage(struct photon * const P, struct geometry const * const G, struct lightCollector const * const LC, struct outputs const *O) {
  FLOATORDBL U[3];
  xyztoXYZ(P->u,LC->theta,LC->phi,U); // U is now the photon trajectory in basis of detection frame (X,Y,Z)
  
  if(U[2] < 0) { // If the Z component of U is negative then the photon is moving towards the light collector plane
    FLOATORDBL resc[3] = {(P->i[0] - G->n[0]/2.0f)*G->d[0] - LC->r[0],
                          (P->i[1] - G->n[1]/2.0f)*G->d[1] - LC->r[1],
                          (P->i[2]               )*G->d[2] - LC->r[2]}; // Photon position relative to the light collector focal plane center when it escapes the cuboid, in the (x,y,z) basis
    
    FLOATORDBL Resc[3];
    xyztoXYZ(resc,LC->theta,LC->phi,Resc); // Resc is now the photon position in the light collector frame (X,Y,Z) when it escapes the cuboid
    
    FLOATORDBL RLCP[2]; // XY coordinates of the point where the photon crosses the light collector plane
    RLCP[0] = Resc[0] - Resc[2]*U[0]/U[2];
    RLCP[1] = Resc[1] - Resc[2]*U[1]/U[2];
    
    FLOATORDBL distLCP = SQRT(RLCP[0]*RLCP[0] + RLCP[1]*RLCP[1]); // Distance between light collector center and the point where the photon crosses the light collector plane
    
    if(distLCP < LC->diam/2) { // If the distance is less than the radius of the light collector
      if(isfinite(LC->f)) { // If the light collector is an objective lens
        FLOATORDBL RImP[2]; // Back-propagated position that the photon would have had in the object plane if propagating freely. This corresponds to where the photon will end up in the image plane for magnification 1x.
        RImP[0] = RLCP[0] + LC->f*U[0]/U[2];
        RImP[1] = RLCP[1] + LC->f*U[1]/U[2];
        FLOATORDBL distImP = SQRT(RImP[0]*RImP[0] + RImP[1]*RImP[1]);
        if(distImP < LC->FSorNA/2) { // If the photon is coming from the area within the Field Size
          long Xindex = (long)(LC->res[0]*(RImP[0]/LC->FSorNA + 1.0f/2));
          long Yindex = (long)(LC->res[0]*(RImP[1]/LC->FSorNA + 1.0f/2));
          long timeindex = LC->res[1] > 1? min(LC->res[1]-1,max(0L,(long)(1+(LC->res[1]-2)*(P->time - (Resc[2] - LC->f)/U[2]*P->RI/C - LC->tStart)/(LC->tEnd - LC->tStart)))): 0; // If we are not measuring time-resolved, LC->res[1] == 1
          atomicAddWrapperDBL(&O->image[Xindex               +
                                        Yindex   *LC->res[0] +
                                        timeindex*LC->res[0]*LC->res[0]],P->weight);
          if(O->NFRdet) for(long i=0;i<P->recordElems;i++) {
            atomicAddWrapperDBL(&O->NFRdet[P->j_record[i]],P->weight*P->weight_record[i]);
          }
        }
      } else { // If the light collector is a fiber tip
        FLOATORDBL thetaLCFF = ATAN(-SQRT(U[0]*U[0] + U[1]*U[1])/U[2]); // Light collector far field polar angle
        if(thetaLCFF < ASIN(min(1.0f,LC->FSorNA))) { // If the photon has an angle within the fiber's NA acceptance
          long timeindex = LC->res[1] > 1? min(LC->res[1]-1,max(0L,(long)(1+(LC->res[1]-2)*(P->time - Resc[2]/U[2]*P->RI/C - LC->tStart)/(LC->tEnd - LC->tStart)))): 0; // If we are not measuring time-resolved, LC->res[1] == 1
          atomicAddWrapperDBL(&O->image[timeindex],P->weight);
          if(O->NFRdet) for(long i=0;i<P->recordElems;i++) {
            atomicAddWrapperDBL(&O->NFRdet[P->j_record[i]],P->weight*P->weight_record[i]);
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
  atomicAddWrapperDBL(&O->FF[(long)FLOOR(theta/PI*G->farFieldRes) + G->farFieldRes*(long)FLOOR(phi_shifted/(2*PI)*G->farFieldRes)],P->weight);
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void formEdgeFluxes(struct photon const * const P, struct geometry const * const G, struct outputs const *O) {
  if(G->boundaryType == 1) {
    if(P->i[2] < 0)             atomicAddWrapperDBL(&O->NI_zneg[(long)P->i[0] + G->n[0]*(long)P->i[1]],P->weight);
    else if(P->i[2] >= G->n[2]) atomicAddWrapperDBL(&O->NI_zpos[(long)P->i[0] + G->n[0]*(long)P->i[1]],P->weight);
    else if(P->i[1] < 0)        atomicAddWrapperDBL(&O->NI_yneg[(long)P->i[0] + G->n[0]*(long)P->i[2]],P->weight);
    else if(P->i[1] >= G->n[1]) atomicAddWrapperDBL(&O->NI_ypos[(long)P->i[0] + G->n[0]*(long)P->i[2]],P->weight);
    else if(P->i[0] < 0)        atomicAddWrapperDBL(&O->NI_xneg[(long)P->i[1] + G->n[1]*(long)P->i[2]],P->weight);
    else if(P->i[0] >= G->n[0]) atomicAddWrapperDBL(&O->NI_xpos[(long)P->i[1] + G->n[1]*(long)P->i[2]],P->weight);
  } else { // boundaryType == 2
    if(P->i[2] < 0)             atomicAddWrapperDBL(&O->NI_zneg[(long)(P->i[0] + G->n[0]*(KILLRANGE-1)/2.0f) + (KILLRANGE*G->n[0])*((long)(P->i[1] + G->n[1]*(KILLRANGE-1)/2.0f))],P->weight);
  }
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void checkEscape(struct photon * const P, struct geometry const * const G, struct lightCollector const * const LC,
        struct outputs const *O) {
  bool escaped = false;
  P->insideVolume = P->i[0] < G->n[0] && P->i[0] >= 0 &&
                    P->i[1] < G->n[1] && P->i[1] >= 0 &&
                    P->i[2] < G->n[2] && P->i[2] >= 0;
  
  switch (G->boundaryType) {
    case 0:
      P->alive = (FABS(P->i[0]/G->n[0] - 1.0f/2) <  KILLRANGE/2.0f &&
                  FABS(P->i[1]/G->n[1] - 1.0f/2) <  KILLRANGE/2.0f &&
                  FABS(P->i[2]/G->n[2] - 1.0f/2) <  KILLRANGE/2.0f);
      break;
    case 1:
      P->alive = P->insideVolume;
      escaped = !P->insideVolume && P->RI == 1;
      break;
    case 2:
      P->alive = (FABS(P->i[0]/G->n[0] - 1.0f/2) <  KILLRANGE/2.0f &&
                  FABS(P->i[1]/G->n[1] - 1.0f/2) <  KILLRANGE/2.0f &&
                       P->i[2]/G->n[2] - 1.0f/2  <  KILLRANGE/2.0f &&
                       P->i[2]                   >= 0);
      escaped = (P->i[2] < 0) && P->RI == 1;
      break;
  }
  
  if(!P->alive && G->boundaryType) formEdgeFluxes(P,G,O);
  if(escaped && O->FF)             formFarField(P,G,O);
  if(escaped && O->image)          formImage(P,G,LC,O); // If image is not NULL then that's because useLightCollector was set to true (non-zero)
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void getNewVoxelProperties(struct photon * const P, struct geometry const * const G) {
  /* Get optical properties of current voxel.
   * If photon is outside cuboid, properties are those of
   * the closest defined voxel. */
  P->j = ((P->i[2] < 0)? 0: ((P->i[2] >= G->n[2])? G->n[2]-1: (long)FLOOR(P->i[2])))*G->n[0]*G->n[1] +
         ((P->i[1] < 0)? 0: ((P->i[1] >= G->n[1])? G->n[1]-1: (long)FLOOR(P->i[1])))*G->n[0]         +
         ((P->i[0] < 0)? 0: ((P->i[0] >= G->n[0])? G->n[0]-1: (long)FLOOR(P->i[0]))); // Index values are restrained to integers in the interval [0,n-1]
  P->mua = G->muav[G->M[P->j]];
  P->mus = G->musv[G->M[P->j]];
  P->g   = G->gv  [G->M[P->j]];
  P->RI  = G->RIv[(P->i[2] < 0)? 0: ((P->i[2] >= G->n[2])? G->n[2]-1: (long)FLOOR(P->i[2]))];
}

#ifdef __NVCC__ // If compiling for CUDA
__device__
#endif
void propagatePhoton(struct photon * const P, struct geometry const * const G, struct outputs const *O, struct paths * const Pa, struct debug *D) {
  long idx;
  
  P->sameVoxel = true;
  
  FLOATORDBL s = min(P->stepLeft/P->mus,min(P->D[0],min(P->D[1],P->D[2])));
  P->stepLeft  = s==P->stepLeft/P->mus? 0: P->stepLeft - s*P->mus; // zero case is to avoid rounding errors
  P->time     += s*P->RI/C;
  
  for(idx=0;idx<3;idx++) {
    long i_old = (long)FLOOR(P->i[idx]);
    long i_new = i_old + SIGN(P->u[idx]);
    if(s == P->D[idx]) { // If we're supposed to go to the voxel boundary along this dimension
      int photondeflection = 0; // Switch that can be 0, 1 or 2. 0 means photon travels straight, 1 means photon is refracted at the voxel boundary, 2 means reflected
      FLOATORDBL cos_new = 0;
      if(idx == 2) { // If we're entering a new z slice
        FLOATORDBL RI_ratio = P->RI/G->RIv[i_new>=0? (i_new<G->n[2]? i_new:G->n[2]-1):0];
        if(RI_ratio != 1) { // If there's a refractive index change
          FLOATORDBL sin_newsq = (P->u[0]*P->u[0] + P->u[1]*P->u[1])*RI_ratio*RI_ratio;
          if(sin_newsq < 1) { // If we don't experience total internal reflection
            cos_new = SIGN(P->u[2])*SQRT(1 - sin_newsq);
            FLOATORDBL R = SQR((RI_ratio*P->u[2] - cos_new)/(RI_ratio*P->u[2] + cos_new))/2 +
                           SQR((RI_ratio*cos_new - P->u[2])/(RI_ratio*cos_new + P->u[2]))/2; // R is the reflectivity assuming equal probability of p or s polarization (unpolarized light at all times)
            photondeflection = RandomNum > R? (FABS(P->u[2]) == 1? 0: 1):2;
          } else photondeflection = 2;
        }
      }
      
      switch(photondeflection) {
        case 0: // Travel straight
          P->i[idx] = (P->u[idx] > 0)? i_old + 1: i_old - FLOATORDBLEPS*(labs(i_old)+1);
          P->sameVoxel = false;
          P->D[idx] = G->d[idx]/FABS(P->u[idx]); // Reset voxel boundary distance
          break;
        case 1: // Refract, can only happen for idx = 2
          P->i[2] = (P->u[2] > 0)? i_old + 1: i_old - FLOATORDBLEPS*(labs(i_old)+1);
          P->sameVoxel = false;
          P->u[0] *= SQRT((1 - cos_new*cos_new)/(1 - P->u[2]*P->u[2]));
          P->D[0] = P->u[0]? (FLOOR(P->i[0]) + (P->u[0]>0) - P->i[0])*G->d[0]/P->u[0]: INFINITY; // Recalculate voxel boundary distance for x
          P->u[1] *= SQRT((1 - cos_new*cos_new)/(1 - P->u[2]*P->u[2]));
          P->D[1] = P->u[1]? (FLOOR(P->i[1]) + (P->u[1]>0) - P->i[1])*G->d[1]/P->u[1]: INFINITY; // Recalculate voxel boundary distance for y
          P->u[2] = cos_new;
          P->D[2] = G->d[2]/FABS(P->u[2]); // Reset voxel boundary distance for z. Note that here it is important that x and y propagations have already been applied before we here modify their D during the z propagation.
          break;
        case 2: // Reflect, can only happen for idx = 2
          P->i[2] = (P->u[2] > 0)? i_old + 1 - FLOATORDBLEPS*(labs(i_old)+1): i_old;
          P->u[2] *= -1;
          P->D[2] = G->d[2]/FABS(P->u[2]); // Reset voxel boundary distance
          break;
      }
    } else { // We're supposed to remain in the same voxel along this dimension
      P->i[idx] += s*P->u[idx]/G->d[idx]; // First take the expected step (including various rounding errors)
      if(FLOOR(P->i[idx]) != i_old) P->i[idx] = (P->u[idx] > 0)? i_old + 1 - FLOATORDBLEPS*(labs(i_old)+1): i_old ; // If photon due to rounding errors actually crossed the border, set it to be barely in the original voxel
      P->D[idx] -= s;
    }
  }
  
  FLOATORDBL absorb = -P->weight*EXPM1(-P->mua*s);   // photon weight absorbed at this step. expm1(x) = exp(x) - 1, accurate even for very small x 
  P->weight -= absorb;             // decrement WEIGHT by amount absorbed

  if(P->insideVolume) {  // only save data if the photon is inside simulation cuboid
    if(O->NFR) {
      atomicAddWrapperDBL(&O->NFR[P->j],absorb);
    }
    if(P->recordSize) { // store indices and weights in pseudosparse array, to later add to NFRdet if photon ends up on the light collector
      if(P->recordElems == P->recordSize) {
        P->recordSize *= 2; // double the record's size
        P->j_record = (long *)reallocWrapper(P->j_record,P->recordSize/2*sizeof(long),P->recordSize*sizeof(long));
        P->weight_record = (FLOATORDBL *)reallocWrapper(P->weight_record,P->recordSize/2*sizeof(FLOATORDBL),P->recordSize*sizeof(FLOATORDBL));
      }
      P->j_record[P->recordElems] = P->j;
      P->weight_record[P->recordElems] = s;
      P->recordElems++;
    }
  }
  #ifdef __NVCC__ // If compiling for CUDA
  if(!threadIdx.x && !blockIdx.x)
  #elif defined(_OPENMP)
  #pragma omp master
  #endif
  {
    if(Pa->nMasterPhotonsLaunched <= Pa->nExamplePaths) {
      if(Pa->pathsElems == Pa->pathsSize) {
        #ifdef __NVCC__ // If compiling for CUDA
        Pa->pathsElems = -1; // Buffer not large enough, has to be resized by host
        #else
        Pa->pathsSize *= 2; // double the record's size
        Pa->data = (FLOATORDBL *)reallocWrapper(Pa->data,2*Pa->pathsSize*sizeof(FLOATORDBL),4*Pa->pathsSize*sizeof(FLOATORDBL));
        #endif
      }
      if(Pa->pathsElems != -1) {
        Pa->data[4*Pa->pathsElems  ] = (P->i[0] - G->n[0]/2.0f)*G->d[0]; // Store x
        Pa->data[4*Pa->pathsElems+1] = (P->i[1] - G->n[1]/2.0f)*G->d[1]; // Store y
        Pa->data[4*Pa->pathsElems+2] = (P->i[2]               )*G->d[2]; // Store z
        Pa->data[4*Pa->pathsElems+3] = P->weight;                       // Store photon weight
        Pa->pathsElems++;
      }
    }
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
void scatterPhoton(struct photon * const P, struct geometry const * const G) {
  // Sample for costheta using Henyey-Greenstein scattering
  FLOATORDBL costheta = FABS(P->g) > SQRT(FLOATORDBLEPS)? (1 + P->g*P->g - SQR((1 - P->g*P->g)/(1 - P->g + 2*P->g*RandomNum)))/(2*P->g) : 2*RandomNum - 1;
  FLOATORDBL sintheta = SQRT(1 - costheta*costheta);
  FLOATORDBL phi = 2*PI*RandomNum;
  FLOATORDBL cosphi = COS(phi);
  FLOATORDBL sinphi = SIN(phi);
  long idx;
  
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
  for(idx=0;idx<3;idx++) P->D[idx] = P->u[idx]? (FLOOR(P->i[idx]) + (P->u[idx]>0) - P->i[idx])*G->d[idx]/P->u[idx] : INFINITY;
  
  P->stepLeft  = -LOG(RandomNum);
}

void normalizeDeposition(struct beam const * const B, struct geometry const * const G, struct lightCollector const * const LC,
        struct outputs const *O) {
  long j;
  double V = G->d[0]*G->d[1]*G->d[2]; // Voxel volume
  long L = G->n[0]*G->n[1]*G->n[2]; // Total number of voxels in cuboid
  long L_LC = LC->res[0]*LC->res[0]; // Total number of spatial pixels in light collector planes
  long L_FF = G->farFieldRes*G->farFieldRes; // Total number of pixels in the far field array
  // Normalize deposition to yield normalizes fluence rate (NFR). For fluorescence, the result is relative to
  // the incident excitation power (not emitted fluorescence power).
  double normfactor = (double)O->nPhotons;
  if(B->S) { // For a 3D source distribution (e.g., fluorescence)
    normfactor /= B->power;
  } else if(B->beamType == 2 && G->boundaryType != 1) { // For infinite plane wave launched into volume without absorbing walls
    normfactor /= KILLRANGE*KILLRANGE;
  }
  
  if(O->NFR)        for(j=0;j<L   ;j++) O->NFR[j]    /= V*normfactor*G->muav[G->M[j]];
  if(O->NFRdet)     for(j=0;j<L   ;j++) O->NFRdet[j] /= V*normfactor;
  if(O->FF)         for(j=0;j<L_FF;j++) O->FF[j]     /=   normfactor;
  if(O->image) {
    if(L_LC > 1) for(j=0;j<L_LC*LC->res[1];j++) O->image[j] /= LC->FSorNA*LC->FSorNA/L_LC*normfactor;
    else         for(j=0;j<     LC->res[1];j++) O->image[j] /= normfactor;
  }
  if(G->boundaryType == 1) {
    for(j=0;j<G->n[1]*G->n[2];j++) O->NI_xpos[j] /= G->d[1]*G->d[2]*normfactor;
    for(j=0;j<G->n[1]*G->n[2];j++) O->NI_xneg[j] /= G->d[1]*G->d[2]*normfactor;
    for(j=0;j<G->n[0]*G->n[2];j++) O->NI_ypos[j] /= G->d[0]*G->d[2]*normfactor;
    for(j=0;j<G->n[0]*G->n[2];j++) O->NI_yneg[j] /= G->d[0]*G->d[2]*normfactor;
    for(j=0;j<G->n[0]*G->n[1];j++) O->NI_zpos[j] /= G->d[0]*G->d[1]*normfactor;
    for(j=0;j<G->n[0]*G->n[1];j++) O->NI_zneg[j] /= G->d[0]*G->d[1]*normfactor;
  } else if(G->boundaryType == 2) {
    for(j=0;j<KILLRANGE*G->n[0]*KILLRANGE*G->n[1];j++) O->NI_zneg[j] /= G->d[0]*G->d[1]*normfactor;
  }
}
