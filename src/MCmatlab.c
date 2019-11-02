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
 * Can be compiled in MATLAB with "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall -pedantic' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -fopenmp -std=c11 -Wall -pedantic' -outdir helperfuncs\private .\src\MCmatlab.c ".\src\libut.lib""
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Go to MATLAB's addon manager and tell it to install the "Support for MinGW-w64 compiler"
 * 2. Type "mex -setup" in the MATLAB command window and ensure that MATLAB has set the C compiler to MinGW64
 * 3. mex should now be able to compile the code using the above command
 *
 ** COMPILING ON MAC
 * As of June 2017, the macOS compiler doesn't support libut (for ctrl+c 
 * breaking) or openmp (for multithreading).
 * Compile in MATLAB with "mex COPTIMFLAGS='$COPTIMFLAGS -Ofast -std=c11 -Wall -pedantic' LDOPTIMFLAGS='$LDOPTIMFLAGS -Ofast -std=c11 -Wall -pedantic' -outdir helperfuncs/private ./src/MCmatlab.c"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Install XCode from the App Store
 * 2. Type "mex -setup" in the MATLAB command window
 ********************************************/

#include "mex.h"
#include <math.h>
#include <time.h>
#include "lambert.c" // For calculating the Lambert W function, originally part of the GNU Scientific Library, created by K. Briggs, G. Jungman and B. Gough and slightly modified by A. Hansen for easier MCmatlab integration
#define DSFMT_MEXP 19937 // Mersenne exponent for dSFMT
#include "dSFMT-src-2.2.3/dSFMT.c" // Double precision SIMD oriented Fast Mersenne Twister(dSFMT)
#ifdef _WIN32 // This is defined on both win32 and win64 systems. We use this preprocessor condition to avoid loading openmp or libut on, e.g., Mac
#include <omp.h>
extern bool utIsInterruptPending(); // Allows catching ctrl+c while executing the mex function
#endif

#define PI          acos(-1)
#define C           29979245800 // speed of light in vacuum in cm/s
#define THRESHOLD   0.01    // used in roulette
#define CHANCE      0.1      // used in roulette
#define SIGN(x)     ((x)>=0 ? 1:-1)
#define RandomNum   dsfmt_genrand_open_close(&P->dsfmt) // Calls for a random number in (0,1]
#define KILLRANGE   5 // Must be odd integer
    /* KILLRANGE determines the region that photons are allowed to stay 
     * alive in in multiples of the cuboid size (if outside, the probability
     * of returning to the region of interest is judged as too low). When
     * launching an infinite plane wave without boundaries, photons will be
     * launched in this whole extended region. */

struct geometry { // Struct type for the constant geometry definitions, including the wavelength-dependent optical properties and the boundary type
  double         d[3];
  long           n[3];
  int            boundaryType;
  double         *muav,*musv,*gv;
  unsigned char  *M;
  double         *RIv;    // Refractive index of each z slice
};

struct beam { // Struct type for the constant beam definitions
  int            beamType;
  int            nearFieldType;
  int            farFieldType;
  double         *S;
  double         power;
  double         waist,divergence;
  double         focus[3];
  double         u[3];
  double         v[3];
};

struct lightCollector { // Struct type for the constant light collector definitions. It can be either an objective lens (for 0<f<INFINITY) or a fiber tip or simple aperture (for f=INFINITY)
  double         r[3]; // Position of center of objective focal plane (not the objective itself) or position of center of the fiber tip
  double         theta; // Polar angle that the objective or fiber is facing
  double         phi; // Azimuthal angle of objective or fiber orientation
  double         f; // Focal length of the objective. If the light collector is a fiber tip, this will be INFINITY.
  double         diam; // Diameter of the objective aperture or core diameter of the fiber. For an ideal thin lens objective, this is 2*tan(arcsin(lensNA/f)).
  double         FSorNA; // For an objective lens: Field Size of the imaging system (diameter of area in object plane that gets imaged). For a fiber tip: The fiber's NA.
  mwSize         res[2]; // Resolution of image plane in pixels along the spatial and time axes. For a fiber, spatial resolution is 1.
  double         tStart; // Start time for the interval used for binned time-resolved detection
  double         tEnd; // End time for the interval used for binned time-resolved detection
};

struct photon { // Struct type for parameters describing the thread-specific current state of a photon
  double         i[3],u[3],D[3]; // Fractional position indices i, ray trajectory unit vector u and distances D to next voxel boundary (yz, xz or yz) along current trajectory
  long           j; // Linear index of current voxel (or closest defined voxel if photon outside cuboid)
  double         mua,mus,g,RI; // Absorption, scattering, anisotropy and refractive index values at current photon position
  double         stepLeft,weight,time;
  bool           insideVolume,alive,sameVoxel;
  dsfmt_t        dsfmt; // "State" of the Mersenne Twister pseudo-random number generator
  long           recordSize; // Current size of the list of voxels in which power has been deposited, used only if calcNFRdet is true
  long           recordElems; // Current number of elements used of the record. Starts at 0 every photon launch, used only if calcNFRdet is true
  long           *j_record; // List of the indices of the voxels in which the current photon has deposited power, used only if calcNFRdet is true
  double         *weight_record; // List of the weights that have been deposited into the voxels, used only if calcNFRdet is true
};

struct paths { // Struct type for storing the paths taken by the nExamplePaths first photons simulated by the master thread
  long           nExamplePaths; // Number of photons to store the path of
  long           nMasterPhotonsLaunched; // Number of photons the master thread has launched
  long           pathsSize; // Current size of the list
  long           pathsElems; // Current number of elements used
  double         *data; // Array containing x, y, z and weight data for the photons, in which the paths of different photons are separated by four NaNs
};

void unitcrossprod(double *a, double *b, double *c) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
  double norm = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);
  for(int idx=0;idx<3;idx++) c[idx] /= norm;
}

void xyztoXYZ(double const * const r, double const theta, double const phi, double * const r_out) {
  // If input r is given in the simulation (x,y,z) frame, then output r_out is that point's coordinates in
  // the (X,Y,Z) light collector frame, where Z is pointing in the direction denoted by the spherical
  // coordinates theta and phi and X lies in the xy plane.
  r_out[0] =            sin(phi)*r[0] -            cos(phi)*r[1];
  r_out[1] = cos(theta)*cos(phi)*r[0] + cos(theta)*sin(phi)*r[1] - sin(theta)*r[2];
  r_out[2] = sin(theta)*cos(phi)*r[0] + sin(theta)*sin(phi)*r[1] + cos(theta)*r[2];
}

void axisrotate(double const * const r, double const * const u, double const theta, double * const r_out) {
  // Rotate the point r an angle theta around vector u, storing result in r_out
  double st = sin(theta), ct = cos(theta);
  
  r_out[0] = (u[0]*u[0]*(1-ct) +      ct)*r[0] + (u[0]*u[1]*(1-ct) - u[2]*st)*r[1] + (u[0]*u[2]*(1-ct) + u[1]*st)*r[2];
  r_out[1] = (u[1]*u[0]*(1-ct) + u[2]*st)*r[0] + (u[1]*u[1]*(1-ct) +      ct)*r[1] + (u[1]*u[2]*(1-ct) - u[0]*st)*r[2];
  r_out[2] = (u[2]*u[0]*(1-ct) - u[1]*st)*r[0] + (u[2]*u[1]*(1-ct) + u[0]*st)*r[1] + (u[2]*u[2]*(1-ct) +      ct)*r[2];
}

void launchPhoton(struct photon * const P, struct beam const * const B, struct geometry const * const G, struct paths * const Pa) {
  double r,s,phi,rand,costheta,sintheta;
  long   d,j,idx;
  double target[3]={0},w0[3];
  
  P->sameVoxel = false;
  P->weight = 1;
  P->recordElems = 0;
  
  do{
    if(B->S) { // If a 3D source distribution was defined
      // ... then search the cumulative distribution function via binary tree method to find the voxel to start the photon in
      rand = RandomNum;
      d = G->n[0]*G->n[1]*G->n[2]-1; // Number of elements in the current section to be searched, minus one
      j = d/2; // Index of the middle element of the current section
      while(!(B->S[j] < rand && B->S[j+1] >= rand)) { // Binary tree search
        if(B->S[j] >= rand) {
          j = j + (d-2)/4 - d/2;
          d = d/2 - 1;
        } else {
          d = d - d/2 - 1;
          j = j + d/2 + 1;
        }
      }
      P->i[0] = j%G->n[0]         + 1 - RandomNum;
      P->i[1] = j/G->n[0]%G->n[1] + 1 - RandomNum;
      P->i[2] = j/G->n[0]/G->n[1] + 1 - RandomNum;
      costheta = 1 - 2*RandomNum;
      sintheta = sqrt(1 - costheta*costheta);
      phi = 2*PI*RandomNum;
      P->u[0] = sintheta*cos(phi);
      P->u[1] = sintheta*sin(phi);
      P->u[2] = costheta;
      P->time = 0;
    } else switch (B->beamType) {
      case -1: // Custom near/far field beam
        switch (B->nearFieldType) {
          case 0: // Gaussian
            phi     = RandomNum*2*PI;
            axisrotate(B->v,B->u,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
            r    = B->waist*sqrt(-0.5*log(RandomNum)); // for target calculation
            for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
            break;
          case 1: // Circular top-hat
            phi     = RandomNum*2*PI;
            axisrotate(B->v,B->u,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
            r    = B->waist*sqrt(RandomNum); // for target calculation
            for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
            break;
          case 2: // Square top-hat
            axisrotate(B->v,B->u,PI/2,w0); // w0 unit vector is now orthogonal to both u and v
            r    = B->waist*(2*RandomNum-1); // for target calculation, displacement along v
            s    = B->waist*(2*RandomNum-1); // for target calculation, displacement along w0
            for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*B->v[idx] + s*w0[idx];
            break;
        }
        switch (B->farFieldType) {
          case 0: // Gaussian
            phi     = RandomNum*2*PI;
            axisrotate(B->v,B->u,phi,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.
            phi     = B->divergence*sqrt(-0.5*log(RandomNum)); // for trajectory calculation. The sqrt is valid within paraxial approximation.
            axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
            break;
          case 1: // Circular top-hat
            phi     = RandomNum*2*PI;
            axisrotate(B->v,B->u,phi,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.
            phi     = B->divergence*sqrt(RandomNum); // for trajectory calculation. The sqrt is valid within paraxial approximation.
            axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
            break;
          case 2: // Cosine distribution (Lambertian)
            phi     = RandomNum*2*PI;
            axisrotate(B->v,B->u,phi,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.
            phi     = asin(sqrt(RandomNum));
            axisrotate(B->u,w0,phi,P->u);
            break;
        }
        P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
        P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
        P->i[2] = 0;
        P->time = -G->RIv[0]/C*sqrt(pow((P->i[0] - G->n[0]/2.0)*G->d[0] - target[0],2) +
                                    pow((P->i[1] - G->n[1]/2.0)*G->d[1] - target[1],2) +
                                    pow((P->i[2]              )*G->d[2] - target[2],2)); // Starting time is set so that the wave crosses the focal plane at time = 0
        break;
      case 0: // pencil beam
        P->i[0] = (B->focus[0] - B->focus[2]*B->u[0]/B->u[2])/G->d[0] + G->n[0]/2.0;
        P->i[1] = (B->focus[1] - B->focus[2]*B->u[1]/B->u[2])/G->d[1] + G->n[1]/2.0;
        P->i[2] = 0;
        for(idx=0;idx<3;idx++) P->u[idx] = B->u[idx];
        P->time = -G->RIv[0]/C*sqrt(pow((P->i[0] - G->n[0]/2.0)*G->d[0] - B->focus[0],2) +
                                    pow((P->i[1] - G->n[1]/2.0)*G->d[1] - B->focus[1],2) +
                                    pow((P->i[2]              )*G->d[2] - B->focus[2],2)); // Starting time is set so that the wave crosses the focal plane at time = 0
        break;
      case 1: // isotropically emitting point source
        P->i[0] = B->focus[0]/G->d[0] + G->n[0]/2.0;
        P->i[1] = B->focus[1]/G->d[1] + G->n[1]/2.0;
        P->i[2] = B->focus[2]/G->d[2];
        costheta = 1 - 2*RandomNum;
        sintheta = sqrt(1 - costheta*costheta);
        phi = 2*PI*RandomNum;
        P->u[0] = sintheta*cos(phi);
        P->u[1] = sintheta*sin(phi);
        P->u[2] = costheta;
        P->time = 0;
        break;
      case 2: // infinite plane wave
        P->i[0] = ((G->boundaryType==1)? 1: KILLRANGE)*G->n[0]*(RandomNum-0.5) + G->n[0]/2.0; // Generates a random ix coordinate within the cuboid
        P->i[1] = ((G->boundaryType==1)? 1: KILLRANGE)*G->n[1]*(RandomNum-0.5) + G->n[1]/2.0; // Generates a random iy coordinate within the cuboid
        P->i[2] = 0;
        for(idx=0;idx<3;idx++) P->u[idx] = B->u[idx];
        P->time = G->RIv[0]/C*((P->i[0] - G->n[0]/2.0)*G->d[0]*B->u[0] +
                               (P->i[1] - G->n[1]/2.0)*G->d[1]*B->u[1] +
                               (P->i[2]              )*G->d[2]*B->u[2]); // Starting time is set so that the wave crosses (x=0,y=0,z=0) at time = 0
        break;
      case 3: // Gaussian focus, Gaussian far field beam
        phi     = RandomNum*2*PI;
        axisrotate(B->v,B->u,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
        r    = B->waist*sqrt(-0.5*log(RandomNum)); // for target calculation
        for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
        phi     = RandomNum*2*PI;
        axisrotate(B->v,B->u,phi,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.
        phi     = B->divergence*sqrt(-0.5*log(RandomNum)); // for trajectory calculation. The sqrt is valid within paraxial approximation.
        axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
        P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
        P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
        P->i[2] = 0;
        P->time = -G->RIv[0]/C*sqrt(pow((P->i[0] - G->n[0]/2.0)*G->d[0] - target[0],2) +
                                    pow((P->i[1] - G->n[1]/2.0)*G->d[1] - target[1],2) +
                                    pow((P->i[2]              )*G->d[2] - target[2],2)); // Starting time is set so that the wave crosses the focal plane at time = 0
        break;
      case 4: // Gaussian focus, top-hat far field beam
        phi     = RandomNum*2*PI;
        axisrotate(B->v,B->u,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
        r    = B->waist*sqrt(-0.5*log(RandomNum)); // for target calculation
        for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
        phi     = RandomNum*2*PI;
        axisrotate(B->v,B->u,phi,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.
        phi     = B->divergence*sqrt(RandomNum); // for trajectory calculation. The sqrt is valid within paraxial approximation.
        axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
        P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
        P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
        P->i[2] = 0;
        P->time = -G->RIv[0]/C*sqrt(pow((P->i[0] - G->n[0]/2.0)*G->d[0] - target[0],2) +
                                    pow((P->i[1] - G->n[1]/2.0)*G->d[1] - target[1],2) +
                                    pow((P->i[2]              )*G->d[2] - target[2],2)); // Starting time is set so that the wave crosses the focal plane at time = 0
        break;
      case 5: // top-hat focus, Gaussian far field beam
        phi     = RandomNum*2*PI;
        axisrotate(B->v,B->u,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
        r    = B->waist*sqrt(RandomNum); // for target calculation
        for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
        phi     = RandomNum*2*PI;
        axisrotate(B->v,B->u,phi,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.
        phi     = B->divergence*sqrt(-0.5*log(RandomNum)); // for trajectory calculation. The sqrt is valid within paraxial approximation.
        axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
        P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
        P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
        P->i[2] = 0;
        P->time = -G->RIv[0]/C*sqrt(pow((P->i[0] - G->n[0]/2.0)*G->d[0] - target[0],2) +
                                    pow((P->i[1] - G->n[1]/2.0)*G->d[1] - target[1],2) +
                                    pow((P->i[2]              )*G->d[2] - target[2],2)); // Starting time is set so that the wave crosses the focal plane at time = 0
        break;
      case 6: // top-hat focus, top-hat far field beam
        phi    = RandomNum*2*PI;
        axisrotate(B->v,B->u,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
        r    = B->waist*sqrt(RandomNum); // for target calculation
        for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
        phi     = RandomNum*2*PI;
        axisrotate(B->v,B->u,phi,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.
        phi     = B->divergence*sqrt(RandomNum); // for trajectory calculation. The sqrt is valid within paraxial approximation.
        axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
        P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
        P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
        P->i[2] = 0;
        P->time = -G->RIv[0]/C*sqrt(pow((P->i[0] - G->n[0]/2.0)*G->d[0] - target[0],2) +
                                    pow((P->i[1] - G->n[1]/2.0)*G->d[1] - target[1],2) +
                                    pow((P->i[2]              )*G->d[2] - target[2],2)); // Starting time is set so that the wave crosses the focal plane at time = 0
        break;
      case 7: // Laguerre-Gaussian LG01 beam
        phi     = RandomNum*2*PI;
        axisrotate(B->v,B->u,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
        r    = B->waist*sqrt((gsl_sf_lambert_Wm1(-RandomNum*exp(-1))+1)/(-2))/1.50087; // for target calculation
        for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
        phi     = B->divergence*sqrt((gsl_sf_lambert_Wm1(-RandomNum*exp(-1))+1)/(-2))/1.50087; // for trajectory calculation. The sqrt is valid within paraxial approximation.
        axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
        P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
        P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
        P->i[2] = 0;
        P->time = -G->RIv[0]/C*sqrt(pow((P->i[0] - G->n[0]/2.0)*G->d[0] - target[0],2) +
                                    pow((P->i[1] - G->n[1]/2.0)*G->d[1] - target[1],2) +
                                    pow((P->i[2]              )*G->d[2] - target[2],2)); // Starting time is set so that the wave crosses the focal plane at time = 0
        break;
    }
    
    switch (G->boundaryType) {
      case 0:
        P->alive = (fabs(P->i[0]/G->n[0] - 1.0/2) <  KILLRANGE/2.0 &&
                    fabs(P->i[1]/G->n[1] - 1.0/2) <  KILLRANGE/2.0 &&
                    fabs(P->i[2]/G->n[2] - 1.0/2) <  KILLRANGE/2.0);
        break;
      case 1:
        P->alive = P->i[0] < G->n[0] && P->i[0] >= 0 &&
                   P->i[1] < G->n[1] && P->i[1] >= 0 &&
                   P->i[2] < G->n[2] && P->i[2] >= 0;
        break;
      case 2:
        P->alive = (fabs(P->i[0]/G->n[0] - 1.0/2) <  KILLRANGE/2.0 &&
                    fabs(P->i[1]/G->n[1] - 1.0/2) <  KILLRANGE/2.0 &&
                         P->i[2]/G->n[2] - 1.0/2  <  KILLRANGE/2.0 &&
                         P->i[2]                  >= 0);
        break;
    }
  } while(!P->alive); // If photon happened to be initialized outside the volume in which it is allowed to travel, we try again
  
  // Calculate distances to next voxel boundary planes
  for(idx=0;idx<3;idx++) P->D[idx] = P->u[idx]? (floor(P->i[idx]) + (P->u[idx]>0) - P->i[idx])*G->d[idx]/P->u[idx] : INFINITY;
  
  P->stepLeft  = -log(RandomNum);
  
  #ifdef _WIN32
  #pragma omp master
  #endif
  {
    Pa->nMasterPhotonsLaunched++;
    if(Pa->nMasterPhotonsLaunched <= Pa->nExamplePaths) {
      if(Pa->pathsElems + 1 >= Pa->pathsSize) {
        Pa->pathsSize *= 2; // double the record's size
        Pa->data = realloc(Pa->data,4*Pa->pathsSize*sizeof(double));
        if(!Pa->data) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");
      }
      for(long i=0;i<4;i++) Pa->data[4*Pa->pathsElems+i] = NAN;
      Pa->pathsElems++;
      Pa->data[4*Pa->pathsElems  ] = (P->i[0] - G->n[0]/2.0)*G->d[0]; // Store x
      Pa->data[4*Pa->pathsElems+1] = (P->i[1] - G->n[1]/2.0)*G->d[1]; // Store y
      Pa->data[4*Pa->pathsElems+2] = (P->i[2]              )*G->d[2]; // Store z
      Pa->data[4*Pa->pathsElems+3] = P->weight;                       // Store photon weight
      Pa->pathsElems++;
    }
  }
}

void formImage(struct photon * const P, struct geometry const * const G, struct lightCollector const * const LC, double * const NFRdet, double * const image) {
  double U[3];
  xyztoXYZ(P->u,LC->theta,LC->phi,U); // U is now the photon trajectory in basis of detection frame (X,Y,Z)
  
  if(U[2] < 0) { // If the Z component of U is negative then the photon is moving towards the light collector plane
    double resc[3] = {(P->i[0] - G->n[0]/2.0)*G->d[0] - LC->r[0],
                      (P->i[1] - G->n[1]/2.0)*G->d[1] - LC->r[1],
                      (P->i[2]              )*G->d[2] - LC->r[2]}; // Photon position relative to the light collector focal plane center when it escapes the cuboid, in the (x,y,z) basis
    
    double Resc[3];
    xyztoXYZ(resc,LC->theta,LC->phi,Resc); // Resc is now the photon position in the light collector frame (X,Y,Z) when it escapes the cuboid
    
    double RLCP[2]; // XY coordinates of the point where the photon crosses the light collector plane
    RLCP[0] = Resc[0] - Resc[2]*U[0]/U[2];
    RLCP[1] = Resc[1] - Resc[2]*U[1]/U[2];
    
    double distLCP = sqrt(RLCP[0]*RLCP[0] + RLCP[1]*RLCP[1]); // Distance between light collector center and the point where the photon crosses the light collector plane
    
    if(distLCP < LC->diam/2) { // If the distance is less than the radius of the light collector
      if(isfinite(LC->f)) { // If the light collector is an objective lens
        double RImP[2]; // Back-propagated position that the photon would have had in the object plane if propagating freely. This corresponds to where the photon will end up in the image plane for magnification 1x.
        RImP[0] = RLCP[0] + LC->f*U[0]/U[2];
        RImP[1] = RLCP[1] + LC->f*U[1]/U[2];
        double distImP = sqrt(RImP[0]*RImP[0] + RImP[1]*RImP[1]);
        if(distImP < LC->FSorNA/2) { // If the photon is coming from the area within the Field Size
          long Xindex = LC->res[0]*(RImP[0]/LC->FSorNA + 1.0/2);
          long Yindex = LC->res[0]*(RImP[1]/LC->FSorNA + 1.0/2);
          long timeindex = LC->res[1] > 1? fmin(LC->res[1]-1,fmax(0,1+(LC->res[1]-2)*(P->time - (Resc[2] - LC->f)/U[2]*P->RI/C - LC->tStart)/(LC->tEnd - LC->tStart))): 0; // If we are not measuring time-resolved, LC->res[1] == 1
          #ifdef _WIN32
          #pragma omp atomic
          #endif
          image[Xindex               +
                Yindex   *LC->res[0] +
                timeindex*LC->res[0]*LC->res[0]] += P->weight;
          if(NFRdet) for(long i=0;i<P->recordElems;i++) {
            #ifdef _WIN32
            #pragma omp atomic
            #endif
            NFRdet[P->j_record[i]] += P->weight_record[i];
          }
        }
      } else { // If the light collector is a fiber tip
        double thetaLCFF = atan(-sqrt(U[0]*U[0] + U[1]*U[1])/U[2]); // Light collector far field polar angle
        if(thetaLCFF < asin(fmin(1,LC->FSorNA))) { // If the photon has an angle within the fiber's NA acceptance
          long timeindex = LC->res[1] > 1? fmin(LC->res[1]-1,fmax(0,1+(LC->res[1]-2)*(P->time - Resc[2]/U[2]*P->RI/C - LC->tStart)/(LC->tEnd - LC->tStart))): 0; // If we are not measuring time-resolved, LC->res[1] == 1
          #ifdef _WIN32
          #pragma omp atomic
          #endif
          image[timeindex] += P->weight;
          if(NFRdet) for(long i=0;i<P->recordElems;i++) {
            #ifdef _WIN32
            #pragma omp atomic
            #endif
            NFRdet[P->j_record[i]] += P->weight_record[i];
          }
        }
      }
    }
  }
}

void formFarField(struct photon const * const P, long const farFieldRes, double * const FF) {
  double theta = (1-DBL_EPSILON)*acos(P->u[2]); // The (1-DBL_EPS) factor is to ensure that photons exiting with theta = PI will be stored correctly
  double phi_shifted = (1-DBL_EPSILON)*(PI + atan2(P->u[1],P->u[0])); // Here it's to handle the case of phi = +PI
  #ifdef _WIN32
  #pragma omp atomic
  #endif
  FF[(long)floor(theta/PI*farFieldRes) + farFieldRes*(long)floor(phi_shifted/(2*PI)*farFieldRes)] += P->weight;
}

void formEdgeFluxes(struct photon const * const P, struct geometry const * const G,
        double * const NI_xpos, double * const NI_xneg, double * const NI_ypos,
        double * const NI_yneg, double * const NI_zpos, double * const NI_zneg) {
  if(G->boundaryType == 1) {
    if(P->i[2] < 0) {
      #ifdef _WIN32
      #pragma omp atomic
      #endif
      NI_zneg[(long)P->i[0] + G->n[0]*(long)P->i[1]] += P->weight;
    } else if(P->i[2] >= G->n[2]) {
      #ifdef _WIN32
      #pragma omp atomic
      #endif
      NI_zpos[(long)P->i[0] + G->n[0]*(long)P->i[1]] += P->weight;
    } else if(P->i[1] < 0) {
      #ifdef _WIN32
      #pragma omp atomic
      #endif
      NI_yneg[(long)P->i[0] + G->n[0]*(long)P->i[2]] += P->weight;
    } else if(P->i[1] >= G->n[1]) {
      #ifdef _WIN32
      #pragma omp atomic
      #endif
      NI_ypos[(long)P->i[0] + G->n[0]*(long)P->i[2]] += P->weight;
    } else if(P->i[0] < 0) {
      #ifdef _WIN32
      #pragma omp atomic
      #endif
      NI_xneg[(long)P->i[1] + G->n[1]*(long)P->i[2]] += P->weight;
    } else if(P->i[0] >= G->n[0]) {
      #ifdef _WIN32
      #pragma omp atomic
      #endif
      NI_xpos[(long)P->i[1] + G->n[1]*(long)P->i[2]] += P->weight;
    }
  } else { // boundaryType == 2
    if(P->i[2] < 0) {
      #ifdef _WIN32
      #pragma omp atomic
      #endif
      NI_zneg[(long)(P->i[0] + G->n[0]*(KILLRANGE-1)/2.0) + (KILLRANGE*G->n[0])*((long)(P->i[1] + G->n[1]*(KILLRANGE-1)/2.0))] += P->weight;
    }
  }
}

void checkEscape(struct photon * const P, struct geometry const * const G, struct lightCollector const * const LC,
        double * const NFRdet, double * const image, long const farFieldRes, double * const FF,
        double * const NI_xpos, double * const NI_xneg, double * const NI_ypos,
        double * const NI_yneg, double * const NI_zpos, double * const NI_zneg) {
  bool escaped = false;
  P->insideVolume = P->i[0] < G->n[0] && P->i[0] >= 0 &&
                    P->i[1] < G->n[1] && P->i[1] >= 0 &&
                    P->i[2] < G->n[2] && P->i[2] >= 0;
  
  switch (G->boundaryType) {
    case 0:
      P->alive = (fabs(P->i[0]/G->n[0] - 1.0/2) <  KILLRANGE/2.0 &&
                  fabs(P->i[1]/G->n[1] - 1.0/2) <  KILLRANGE/2.0 &&
                  fabs(P->i[2]/G->n[2] - 1.0/2) <  KILLRANGE/2.0);
      break;
    case 1:
      P->alive = P->insideVolume;
      escaped = !P->insideVolume && P->RI == 1;
      break;
    case 2:
      P->alive = (fabs(P->i[0]/G->n[0] - 1.0/2) <  KILLRANGE/2.0 &&
                  fabs(P->i[1]/G->n[1] - 1.0/2) <  KILLRANGE/2.0 &&
                       P->i[2]/G->n[2] - 1.0/2  <  KILLRANGE/2.0 &&
                       P->i[2]                  >= 0);
      escaped = (P->i[2] < 0) && P->RI == 1;
      break;
  }
  
  if(!P->alive && G->boundaryType) formEdgeFluxes(P,G,NI_xpos,NI_xneg,NI_ypos,NI_yneg,NI_zpos,NI_zneg);
  if(escaped && FF)                formFarField(P,farFieldRes,FF);
  if(escaped && image)             formImage(P,G,LC,NFRdet,image); // If image is not NULL then that's because useLightCollector was set to true (non-zero)
}

void getNewVoxelProperties(struct photon * const P, struct geometry const * const G) {
  /* Get optical properties of current voxel.
   * If photon is outside cuboid, properties are those of
   * the closest defined voxel. */
  P->j = ((P->i[2] < 0)? 0: ((P->i[2] >= G->n[2])? G->n[2]-1: floor(P->i[2])))*G->n[0]*G->n[1] +
         ((P->i[1] < 0)? 0: ((P->i[1] >= G->n[1])? G->n[1]-1: floor(P->i[1])))*G->n[0]         +
         ((P->i[0] < 0)? 0: ((P->i[0] >= G->n[0])? G->n[0]-1: floor(P->i[0]))); // Index values are restrained to integers in the interval [0,n-1]
  P->mua = G->muav[G->M[P->j]];
  P->mus = G->musv[G->M[P->j]];
  P->g   = G->gv  [G->M[P->j]];
  P->RI  = G->RIv[(P->i[2] < 0)? 0: ((P->i[2] >= G->n[2])? G->n[2]-1: (long)floor(P->i[2]))];
}

void propagatePhoton(struct photon * const P, struct geometry const * const G, double * const NFR, struct paths * const Pa) {
  long idx;
  
  P->sameVoxel = true;
  
  double s = fmin(P->stepLeft/P->mus,fmin(P->D[0],fmin(P->D[1],P->D[2])));
  P->stepLeft  = s==P->stepLeft/P->mus? 0: P->stepLeft - s*P->mus; // zero case is to avoid rounding errors
  P->time     += s*P->RI/C;
  
  for(idx=0;idx<3;idx++) {
    long i_old = floor(P->i[idx]);
    long i_new = i_old + SIGN(P->u[idx]);
    if(s == P->D[idx]) { // If we're supposed to go to the voxel boundary along this dimension
      int photondeflection = 0; // Switch that can be 0, 1 or 2. 0 means photon travels straight, 1 means photon is refracted at the voxel boundary, 2 means reflected
      double cos_new = 0;
      if(idx == 2) { // If we're entering a new z slice
        double RI_ratio = P->RI/G->RIv[i_new>=0? (i_new<G->n[2]? i_new:G->n[2]-1):0];
        if(RI_ratio != 1) { // If there's a refractive index change
          double sin_newsq = (P->u[0]*P->u[0] + P->u[1]*P->u[1])*RI_ratio*RI_ratio;
          if(sin_newsq < 1) { // If we don't experience total internal reflection
            cos_new = SIGN(P->u[2])*sqrt(1 - sin_newsq);
            double R = pow((RI_ratio*P->u[2] - cos_new)/(RI_ratio*P->u[2] + cos_new),2)/2 +
                       pow((RI_ratio*cos_new - P->u[2])/(RI_ratio*cos_new + P->u[2]),2)/2; // R is the reflectivity assuming equal probability of p or s polarization (unpolarized light at all times)
            photondeflection = RandomNum > R? (fabs(P->u[2]) == 1? 0: 1):2;
          } else photondeflection = 2;
        }
      }
      
      switch(photondeflection) {
        case 0: // Travel straight
          P->i[idx] = (P->u[idx] > 0)? i_old + 1: i_old - DBL_EPSILON*(labs(i_old)+1);
          P->sameVoxel = false;
          P->D[idx] = G->d[idx]/fabs(P->u[idx]); // Reset voxel boundary distance
          break;
        case 1: // Refract, can only happen for idx = 2
          P->i[2] = (P->u[2] > 0)? i_old + 1: i_old - DBL_EPSILON*(labs(i_old)+1);
          P->sameVoxel = false;
          double scalefactor = sqrt((1 - cos_new*cos_new)/(1 - P->u[2]*P->u[2])); // Here we already know that P->u[2] is not +-1
          if(P->u[0]) {
            P->u[0] *= scalefactor;
            P->D[0] = (floor(P->i[0]) + (P->u[0]>0) - P->i[0])*G->d[0]/P->u[0]; // Recalculate voxel boundary distance for x
          }
          if(P->u[1]) {
            P->u[1] *= scalefactor;
            P->D[1] = (floor(P->i[1]) + (P->u[1]>0) - P->i[1])*G->d[1]/P->u[1]; // Recalculate voxel boundary distance for y
          }
          P->u[2] = cos_new;
          P->D[2] = G->d[2]/fabs(P->u[2]); // Reset voxel boundary distance for z. Note that here it is important that x and y propagations have already been applied before we here modify their D during the z propagation.
          break;
        case 2: // Reflect, can only happen for idx = 2
          P->i[2] = (P->u[2] > 0)? i_old + 1 - DBL_EPSILON*(labs(i_old)+1): i_old;
          P->u[2] *= -1;
          P->D[2] = G->d[2]/fabs(P->u[2]); // Reset voxel boundary distance
          break;
      }
    } else { // We're supposed to remain in the same voxel along this dimension
      P->i[idx] += s*P->u[idx]/G->d[idx]; // First take the expected step (including various rounding errors)
      if(floor(P->i[idx]) != i_old) P->i[idx] = (P->u[idx] > 0)? i_old + 1 - DBL_EPSILON*(labs(i_old)+1): i_old ; // If photon due to rounding errors actually crossed the border, set it to be barely in the original voxel
      P->D[idx] -= s;
    }
  }
  
  double absorb = P->weight*(1 - exp(-P->mua*s));   // photon weight absorbed at this step
  P->weight -= absorb;             // decrement WEIGHT by amount absorbed
  if(P->insideVolume) {  // only save data if the photon is inside simulation cuboid
    if(NFR) {
      #ifdef _WIN32
      #pragma omp atomic
      #endif
      NFR[P->j] += absorb;
    }
    if(P->recordSize) { // store indices and weights in pseudosparse array, to later add to NFRdet if photon ends up on the light collector
      if(P->recordElems == P->recordSize) {
        P->recordSize *= 2; // double the record's size
        P->j_record = realloc(P->j_record,P->recordSize*sizeof(long));
        P->weight_record = realloc(P->weight_record,P->recordSize*sizeof(double));
        if(!P->j_record || !P->weight_record) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");
      }
      P->j_record[P->recordElems] = P->j;
      P->weight_record[P->recordElems] = absorb;
      P->recordElems++;
    }
  }
  #ifdef _WIN32
  #pragma omp master
  #endif
  {
    if(Pa->nMasterPhotonsLaunched <= Pa->nExamplePaths) {
      if(Pa->pathsElems == Pa->pathsSize) {
        Pa->pathsSize *= 2; // double the record's size
        Pa->data = realloc(Pa->data,4*Pa->pathsSize*sizeof(double));
        if(!Pa->data) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");
      }
      Pa->data[4*Pa->pathsElems  ] = (P->i[0] - G->n[0]/2.0)*G->d[0]; // Store x
      Pa->data[4*Pa->pathsElems+1] = (P->i[1] - G->n[1]/2.0)*G->d[1]; // Store y
      Pa->data[4*Pa->pathsElems+2] = (P->i[2]              )*G->d[2]; // Store z
      Pa->data[4*Pa->pathsElems+3] = P->weight;                       // Store photon weight
      Pa->pathsElems++;
    }
  }
}

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

void scatterPhoton(struct photon * const P, struct geometry const * const G) {
  // Sample for costheta using Henyey-Greenstein scattering
  double costheta = fabs(P->g) > sqrt(DBL_EPSILON)? (1 + P->g*P->g - pow((1 - P->g*P->g)/(1 - P->g + 2*P->g*RandomNum),2))/(2*P->g) : 2*RandomNum - 1;
  double sintheta = sqrt(1 - costheta*costheta);
  double phi = 2*PI*RandomNum;
  double cosphi = cos(phi);
  double sinphi = sin(phi);
  long idx;
  
  if(fabs(P->u[2]) < 1) {
    double ux_temp =  sintheta*(P->u[0]*P->u[2]*cosphi - P->u[1]*sinphi)/sqrt(P->u[0]*P->u[0] + P->u[1]*P->u[1]) + P->u[0]*costheta;
    double uy_temp =  sintheta*(P->u[1]*P->u[2]*cosphi + P->u[0]*sinphi)/sqrt(P->u[0]*P->u[0] + P->u[1]*P->u[1]) + P->u[1]*costheta;
    P->u[2]        = -sintheta*(                cosphi                 )*sqrt(P->u[0]*P->u[0] + P->u[1]*P->u[1]) + P->u[2]*costheta;
    P->u[1]        = uy_temp;
    P->u[0]        = ux_temp;
  } else {
    P->u[0] = sintheta*cosphi;
    P->u[1] = sintheta*sinphi;
    P->u[2] = costheta*SIGN(P->u[2]);
  }
  
  // Calculate distances to next voxel boundary planes
  for(idx=0;idx<3;idx++) P->D[idx] = P->u[idx]? (floor(P->i[idx]) + (P->u[idx]>0) - P->i[idx])*G->d[idx]/P->u[idx] : INFINITY;
  
  P->stepLeft  = -log(RandomNum);
}

void normalizeDeposition(struct beam const * const B, struct geometry const * const G, struct lightCollector const * const LC,
        double * const NFR, double * const NFRdet, double * const image, double * const FF, long const farFieldRes,
        double * const NI_xpos, double * const NI_xneg, double * const NI_ypos,
        double * const NI_yneg, double * const NI_zpos, double * const NI_zneg, double nPhotons) {
  long j;
  double V = G->d[0]*G->d[1]*G->d[2]; // Voxel volume
  long L = G->n[0]*G->n[1]*G->n[2]; // Total number of voxels in cuboid
  long L_LC = LC->res[0]*LC->res[0]; // Total number of spatial pixels in light collector planes
  long L_FF = farFieldRes*farFieldRes; // Total number of pixels in the far field array
  // Normalize deposition to yield normalizes fluence rate (NFR). For fluorescence, the result is relative to
  // the incident excitation power (not emitted fluorescence power).
  double normfactor = nPhotons;
  if(B->S) { // For a 3D source distribution (e.g., fluorescence)
    normfactor /= B->power;
    free(B->S);
  } else if(B->beamType == 2 && G->boundaryType != 1) { // For infinite plane wave launched into volume without absorbing walls
    normfactor /= KILLRANGE*KILLRANGE;
  }
  
  if(NFR)        for(j=0;j<L   ;j++) NFR[j]    /= V*normfactor*G->muav[G->M[j]];
  if(NFRdet)     for(j=0;j<L   ;j++) NFRdet[j] /= V*normfactor*G->muav[G->M[j]];
  if(FF)         for(j=0;j<L_FF;j++) FF[j]   /=   normfactor;
  if(image) {
    if(L_LC > 1) for(j=0;j<L_LC*LC->res[1];j++) image[j] /= LC->FSorNA*LC->FSorNA/L_LC*normfactor;
    else         for(j=0;j<     LC->res[1];j++) image[j] /= normfactor;
  }
  if(G->boundaryType == 1) {
    for(j=0;j<G->n[1]*G->n[2];j++) NI_xpos[j] /= G->d[1]*G->d[2]*normfactor;
    for(j=0;j<G->n[1]*G->n[2];j++) NI_xneg[j] /= G->d[1]*G->d[2]*normfactor;
    for(j=0;j<G->n[0]*G->n[2];j++) NI_ypos[j] /= G->d[0]*G->d[2]*normfactor;
    for(j=0;j<G->n[0]*G->n[2];j++) NI_yneg[j] /= G->d[0]*G->d[2]*normfactor;
    for(j=0;j<G->n[0]*G->n[1];j++) NI_zpos[j] /= G->d[0]*G->d[1]*normfactor;
    for(j=0;j<G->n[0]*G->n[1];j++) NI_zneg[j] /= G->d[0]*G->d[1]*normfactor;
  } else if(G->boundaryType == 2) {
    for(j=0;j<KILLRANGE*G->n[0]*KILLRANGE*G->n[1];j++) NI_zneg[j] /= G->d[0]*G->d[1]*normfactor;
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
  long idx;                  // General-purpose non-thread-specific index variable
  int  pctProgress = 0;      // Simulation progress in percent
  int  newPctProgress;
  bool ctrlc_caught = false; // Has a ctrl+c been passed from MATLAB?
  
  bool simFluorescence = *mxGetPr(prhs[1]) == 2;
  mxArray *MatlabMC = mxGetField(prhs[0],0,simFluorescence? "FMC": "MC");
  
  bool silentMode = mxIsLogicalScalarTrue(mxGetField(MatlabMC,0,"silentMode"));
  bool calcNFR    = mxIsLogicalScalarTrue(mxGetField(MatlabMC,0,"calcNFR")); // Are we supposed to calculate the NFR matrix?
  bool calcNFRdet = mxIsLogicalScalarTrue(mxGetField(MatlabMC,0,"calcNFRdet")); // Are we supposed to calculate the NFRdet matrix?
  
  // To know if we are simulating fluorescence, we check if a "sourceDistribution" field exists. If so, we will use it later in the beam definition.
  mxArray *MatlabBeam = mxGetField(MatlabMC,0,"beam");
  double *S_PDF       = mxGetData(mxGetField(MatlabBeam,0,"sourceDistribution")); // Power emitted by the individual voxels per unit volume. Can be percieved as an unnormalized probability density function of the 3D source distribution

  // Variables for timekeeping and number of photons
  bool            simulationTimed = mxIsNaN(*mxGetPr(mxGetField(MatlabMC,0,"nPhotonsRequested")));
  double          simulationTimeRequested = simulationTimed? *mxGetPr(mxGetField(MatlabMC,0,"simulationTime")): INFINITY;
  double          simulationTimeSpent;
  struct timespec simulationTimeStart, simulationTimeCurrent;
  long long       nPhotonsRequested = simulationTimed? -1: *mxGetPr(mxGetField(MatlabMC,0,"nPhotonsRequested"));

  // Geometry struct definition
  mxArray *mediaProperties = mxGetField(MatlabMC,0,"mediaProperties");
  int     nM = mxGetN(mediaProperties);
  double   muav[nM];            // muav[0:nM-1], absorption coefficient of ith medium
  double   musv[nM];            // scattering coeff.
  double   gv[nM];              // anisotropy of scattering
  for(idx=0;idx<nM;idx++) {
    muav[idx] = *mxGetPr(mxGetField(mediaProperties,idx,"mua"));
    musv[idx] = *mxGetPr(mxGetField(mediaProperties,idx,"mus"));
    gv[idx]   = *mxGetPr(mxGetField(mediaProperties,idx,"g"));
  }
  
  mwSize const *dimPtr = mxGetDimensions(mxGetField(MatlabMC,0,"M"));
  long L = dimPtr[0]*dimPtr[1]*dimPtr[2]; // Total number of voxels in cuboid
  mxArray *MatlabG         = mxGetField(prhs[0],0,"G");
  struct geometry const G_var = (struct geometry) {
    .d = {*mxGetPr(mxGetField(MatlabG,0,"dx")),
          *mxGetPr(mxGetField(MatlabG,0,"dy")),
          *mxGetPr(mxGetField(MatlabG,0,"dz"))},
    .n = {dimPtr[0], dimPtr[1], dimPtr[2]},
    .boundaryType = *mxGetPr(mxGetField(MatlabMC,0,"boundaryType")),
    .muav = muav,
    .musv = musv,
    .gv = gv,
    .M = malloc(L*sizeof(unsigned char)),
    .RIv = mxGetData(mxGetField(MatlabMC,0,"RI"))
  };
  struct geometry const *G = &G_var;
  unsigned char *M_matlab = mxGetData(mxGetField(MatlabMC,0,"M"));
  for(idx=0;idx<L;idx++) G->M[idx] = M_matlab[idx] - 1; // Convert from MATLAB 1-based indexing to C 0-based indexing
  
  // Beam struct definition
  double power        = 0;
  double *S           = NULL; // Cumulative distribution function
  if(S_PDF) {
    S = malloc((L+1)*sizeof(double));
    if(!S) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");
    S[0] = 0;
    for(idx=1;idx<(L+1);idx++) S[idx] = S[idx-1] + S_PDF[idx-1];
    power = S[L]*G->d[0]*G->d[1]*G->d[2];
    for(idx=1;idx<(L+1);idx++) S[idx] /= S[L];
  }
  
  double tb = S? 0: *mxGetPr(mxGetField(MatlabBeam,0,"theta"));
  double pb = S? 0: *mxGetPr(mxGetField(MatlabBeam,0,"phi"));

  double u[3] = {sin(tb)*cos(pb),
  sin(tb)*sin(pb),
  cos(tb)}; // Temporary array
  
  double v[3] = {1,0,0};
  if(u[2]!=1) unitcrossprod(u,(double[]){0,0,1},v);
  
  bool customNearFarField = S? 0: mxIsNaN(*mxGetPr(mxGetField(MatlabBeam,0,"beamType")));
  struct beam const B_var = (struct beam) {
    .beamType      =  S? 0: ( customNearFarField? -1: *mxGetPr(mxGetField(MatlabBeam,0,"beamType"))),
    .nearFieldType =  S? 0: (!customNearFarField? -1: *mxGetPr(mxGetField(MatlabBeam,0,"nearFieldType"))),
    .farFieldType  =  S? 0: (!customNearFarField? -1: *mxGetPr(mxGetField(MatlabBeam,0,"farFieldType"))),
    .S             =  S,
    .power         =  power,
    .waist         =  S? 0: *mxGetPr(mxGetField(MatlabBeam,0,"waist")),
    .divergence    =  S? 0: *mxGetPr(mxGetField(MatlabBeam,0,"divergence")),
    .focus         = {S? 0: *mxGetPr(mxGetField(MatlabBeam,0,"xFocus")),
                      S? 0: *mxGetPr(mxGetField(MatlabBeam,0,"yFocus")),
                      S? 0: *mxGetPr(mxGetField(MatlabBeam,0,"zFocus"))},
    .u             = {u[0],u[1],u[2]},
    .v             = {v[0],v[1],v[2]} // normal vector to beam center axis
  };
  struct beam const *B = &B_var;

  // Light Collector struct definition
  mxArray *MatlabLC      = mxGetField(MatlabMC,0,"LC");
  bool useLightCollector = mxIsLogicalScalarTrue(mxGetField(MatlabMC,0,"useLightCollector"));
  
  double theta     = *mxGetPr(mxGetField(MatlabLC,0,"theta"));
  double phi       = *mxGetPr(mxGetField(MatlabLC,0,"phi"));
  double f         = *mxGetPr(mxGetField(MatlabLC,0,"f"));
  long   nTimeBins = S? 0: *mxGetPr(mxGetField(MatlabLC,0,"nTimeBins"));

  struct lightCollector const LC_var = (struct lightCollector) {
    .r            = {*mxGetPr(mxGetField(MatlabLC,0,"x")) - (isfinite(f)? f*sin(theta)*cos(phi):0), // xyz coordinates of center of light collector
                     *mxGetPr(mxGetField(MatlabLC,0,"y")) - (isfinite(f)? f*sin(theta)*sin(phi):0),
                     *mxGetPr(mxGetField(MatlabLC,0,"z")) - (isfinite(f)? f*cos(theta)         :0)},
    .theta        =  theta,
    .phi          =  phi,
    .f            =  f,
    .diam         =  *mxGetPr(mxGetField(MatlabLC,0,"diam")),
    .FSorNA       =  *mxGetPr(mxGetField(MatlabLC,0,isfinite(f)?"fieldSize":"NA")),
    .res          = {*mxGetPr(mxGetField(MatlabLC,0,"res")),
                     nTimeBins? nTimeBins+2: 1},
    .tStart       =  S? 0: *mxGetPr(mxGetField(MatlabLC,0,"tStart")),
    .tEnd         =  S? 0: *mxGetPr(mxGetField(MatlabLC,0,"tEnd"))
  };
  struct lightCollector const *LC = &LC_var;

  // Paths definitions (example photon trajectories)
  struct paths Pa_var;
  struct paths *Pa = &Pa_var;
  Pa->nExamplePaths = *mxGetPr(mxGetField(MatlabMC,0,"nExamplePaths")); // How many photon path examples are we supposed to store?
  Pa->nMasterPhotonsLaunched = 0;
  Pa->pathsSize = 10*Pa->nExamplePaths; // Will be dynamically increased later if needed
  Pa->pathsElems = 0;
  Pa->data = malloc(4*Pa->pathsSize*sizeof(double)); // Will be NULL if Pa->pathsSize == 0
  if(Pa->pathsSize && !Pa->data) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");
  
  // Far field angle distribution resolution
  long farFieldRes = *mxGetPr(mxGetField(MatlabMC,0,"farFieldRes"));
  
  // Prepare output MATLAB struct
  plhs[0] = mxDuplicateArray(prhs[0]);
  
  mxArray *MCout = mxGetField(plhs[0],0,S? "FMC": "MC");
  mxArray *LCout = mxGetField(MCout,0,"LC");
  mxDestroyArray(mxGetField(MCout,0,"nPhotons"));
  mxDestroyArray(mxGetField(MCout,0,"nThreads"));
  
  if(useLightCollector) mxDestroyArray(mxGetField(LCout,0,"image"));
  if(calcNFR) mxDestroyArray(mxGetField(MCout,0,"NFR"));
  if(calcNFRdet) mxDestroyArray(mxGetField(MCout,0,"NFRdet"));
  if(farFieldRes) mxDestroyArray(mxGetField(MCout,0,"farField"));
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
  if(useLightCollector) mxSetField(LCout,0,"image", mxCreateNumericArray(3,(mwSize[]){LC->res[0],LC->res[0],LC->res[1]},mxDOUBLE_CLASS,mxREAL));
  if(farFieldRes)       mxSetField(MCout,0,"farField", mxCreateDoubleMatrix(farFieldRes,farFieldRes,mxREAL));
  if(G->boundaryType == 1) {
    mxSetField(MCout,0,"NI_xpos", mxCreateDoubleMatrix(G->n[1],G->n[2],mxREAL));
    mxSetField(MCout,0,"NI_xneg", mxCreateDoubleMatrix(G->n[1],G->n[2],mxREAL));
    mxSetField(MCout,0,"NI_ypos", mxCreateDoubleMatrix(G->n[0],G->n[2],mxREAL));
    mxSetField(MCout,0,"NI_yneg", mxCreateDoubleMatrix(G->n[0],G->n[2],mxREAL));
    mxSetField(MCout,0,"NI_zpos", mxCreateDoubleMatrix(G->n[0],G->n[1],mxREAL));
    mxSetField(MCout,0,"NI_zneg", mxCreateDoubleMatrix(G->n[0],G->n[1],mxREAL));
  } else if(G->boundaryType == 2) mxSetField(MCout,0,"NI_zneg", mxCreateDoubleMatrix(KILLRANGE*G->n[0],KILLRANGE*G->n[1],mxREAL));
  mxSetField(MCout,0,"nPhotons",mxCreateDoubleMatrix(1,1,mxREAL));
  mxSetField(MCout,0,"nThreads",mxCreateDoubleMatrix(1,1,mxREAL));
  double *NFR         = calcNFR? mxGetPr(mxGetField(MCout,0,"NFR")): NULL;
  double *NFRdet      = calcNFRdet? mxGetPr(mxGetField(MCout,0,"NFRdet")): NULL;
  double *image       = useLightCollector? mxGetPr(mxGetField(LCout,0,"image")): NULL;
  double *nPhotonsPtr = mxGetPr(mxGetField(MCout,0,"nPhotons"));
  double *nThreadsPtr = mxGetPr(mxGetField(MCout,0,"nThreads"));
  double *FF          = farFieldRes? mxGetPr(mxGetField(MCout,0,"farField")): NULL;
  double *NI_xpos     = G->boundaryType == 1? mxGetPr(mxGetField(MCout,0,"NI_xpos")): NULL;
  double *NI_xneg     = G->boundaryType == 1? mxGetPr(mxGetField(MCout,0,"NI_xneg")): NULL;
  double *NI_ypos     = G->boundaryType == 1? mxGetPr(mxGetField(MCout,0,"NI_ypos")): NULL;
  double *NI_yneg     = G->boundaryType == 1? mxGetPr(mxGetField(MCout,0,"NI_yneg")): NULL;
  double *NI_zpos     = G->boundaryType == 1? mxGetPr(mxGetField(MCout,0,"NI_zpos")): NULL;
  double *NI_zneg     = G->boundaryType != 0? mxGetPr(mxGetField(MCout,0,"NI_zneg")): NULL;
  
  
  if(!silentMode) {
    // Display progress indicator
    printf(simFluorescence?"-----------Fluorescence Monte Carlo Simulation-----------\n":
                           "-----------------Monte Carlo Simulation------------------\n");
    if(simulationTimed) printf("Simulation duration = %0.2f min\n",simulationTimeRequested);
    else                printf("Requested # of photons = %0.2e\n",(double)nPhotonsRequested);
    printf("Calculating...   0%% done");
    mexEvalString("drawnow;");
  }
  
  // ============================ MAJOR CYCLE ========================
  clock_gettime(CLOCK_MONOTONIC, &simulationTimeStart);
  #ifdef _WIN32
  bool useAllCPUs = mxIsLogicalScalarTrue(mxGetField(MatlabMC,0,"useAllCPUs"));
  *nThreadsPtr = useAllCPUs? omp_get_num_procs(): fmax(omp_get_num_procs()-1,1);
  #pragma omp parallel num_threads((long)*nThreadsPtr)
  #else
  *nThreadsPtr = 1;
  #endif
  {
    struct photon P_var;
    struct photon *P = &P_var;
    
    P->recordSize    = NFRdet? 1000: 0; // If we're supposed to calculate NFRdet, we start the record at a size of 1000 elements - it will be dynamically expanded later if needed
    P->j_record      = malloc(P->recordSize*sizeof(long)); // Will be NULL if P->recordSize == 0
    if(P->recordSize && !P->j_record) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");
    P->weight_record = malloc(P->recordSize*sizeof(double)); // Will be NULL if P->recordSize == 0
    if(P->recordSize && !P->weight_record) mexErrMsgIdAndTxt("MCmatlab:OutOfMemory","Error: Out of memory");
    
    #ifdef _WIN32
    int threadnum = omp_get_thread_num();
    #else
    int threadnum = 0;
    #endif
    dsfmt_init_gen_rand(&P->dsfmt,simulationTimeStart.tv_nsec + threadnum); // Seed the photon's random number generator
    
    while((simulationTimed? (pctProgress < 100) : (*nPhotonsPtr < nPhotonsRequested - threadnum)) && !ctrlc_caught) {
      #ifdef _WIN32
      #pragma omp atomic
      #endif
      *nPhotonsPtr = *nPhotonsPtr + 1;

      launchPhoton(P,B,G,Pa);
      
      while(P->alive) {
        while(P->stepLeft>0) {
          if(!P->sameVoxel) {
            checkEscape(P,G,LC,NFRdet,image,farFieldRes,FF,NI_xpos,NI_xneg,NI_ypos,NI_yneg,NI_zpos,NI_zneg);
            if(!P->alive) break;
            getNewVoxelProperties(P,G); // If photon has just entered a new voxel or has just been launched
          }
          propagatePhoton(P,G,NFR,Pa);
        }
        checkRoulette(P);
        scatterPhoton(P,G);
      }

      #ifdef _WIN32
      #pragma omp master
      #endif
      {
        // Check whether ctrl+c has been pressed
        #ifdef _WIN32
        if(utIsInterruptPending()) {
          ctrlc_caught = true;
          printf("\nCtrl+C detected, stopping.");
        }
        #endif
        clock_gettime(CLOCK_MONOTONIC, &simulationTimeCurrent);
        
        // Print out message about progress.
        if(simulationTimed) {
          newPctProgress = 100*(simulationTimeCurrent.tv_sec  - simulationTimeStart.tv_sec +
                               (simulationTimeCurrent.tv_nsec - simulationTimeStart.tv_nsec)/1e9) / (simulationTimeRequested*60);
        } else {
          newPctProgress = 100*(*nPhotonsPtr)/nPhotonsRequested;
        }
        if((newPctProgress != pctProgress) && !ctrlc_caught) {
          pctProgress = newPctProgress;
          if(!silentMode) {
            printf("\b\b\b\b\b\b\b\b\b%3.i%% done", pctProgress<100? pctProgress: 100);
            mexEvalString("drawnow;");
          }
        }
      }
    }
    
    free(P->j_record); // Will do nothing if P->j_record == NULL
    free(P->weight_record); // Will do nothing if P->weight_record == NULL
  }
  
  clock_gettime(CLOCK_MONOTONIC, &simulationTimeCurrent);
  simulationTimeSpent = simulationTimeCurrent.tv_sec  - simulationTimeStart.tv_sec +
                       (simulationTimeCurrent.tv_nsec - simulationTimeStart.tv_nsec)/1e9;
  if(!silentMode) {
    printf("\nSimulated %0.2e photons at a rate of %0.2e photons per minute\n",*nPhotonsPtr, *nPhotonsPtr*60/simulationTimeSpent);
    mexEvalString("drawnow;");
  }
  normalizeDeposition(B,G,LC,NFR,NFRdet,image,FF,farFieldRes,NI_xpos,NI_xneg,NI_ypos,NI_yneg,NI_zpos,NI_zneg,*nPhotonsPtr); // Convert data to relative fluence rate
  free(G->M);
  
  if(Pa->nExamplePaths) {
    mxDestroyArray(mxGetField(MCout,0,"examplePaths"));
    mxSetField(MCout,0,"examplePaths",mxCreateNumericMatrix(4,Pa->pathsElems,mxDOUBLE_CLASS,mxREAL));
    memcpy(mxGetPr(mxGetField(MCout,0,"examplePaths")),Pa->data,4*sizeof(double)*Pa->pathsElems);
    free(Pa->data);
  }
}
