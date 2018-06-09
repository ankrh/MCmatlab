/********************************************
 *
 *  mcxyz.c,	in the C programming language
 *
 * created 2010, 2012 by
 *	Steven L. JACQUES
 *  Ting LI
 *	Oregon Health & Science University
 *
 * Log:
 *  2010: Written by Ting Li based on Steve's mcsub.c.
 *  2010-12-30: Use Steve Jacques' FindVoxelFace().
 *  2012-05-08: Reorganized by Steve Jacques.
 *  2014-01-09: Edited to included a uniformly distributed light source over the entire surface by Mathias Christensen, DTU Fotonik.
 *  2017-04-20: Overhauled by Anders K. Hansen and Dominik Marti, DTU Fotonik. Fundamental method remained unchanged.
 *      Uses the Mersenne Twister for random number generation.
 *  2017-06-07: Adapted to MATLAB mex file generation by Anders K. Hansen, DTU Fotonik
 *  2018-04-19: Added support for illuminating with an isotropically emitting 3D distribution of sources, such as a collection of fluorescing emitters
 *
 ** COMPILING ON WINDOWS
 * Can be compiled in MATLAB with "mex COPTIMFLAGS='$COPTIMFLAGS -O3 -fopenmp -std=c11 -Wall -pedantic' LDOPTIMFLAGS='$LDOPTIMFLAGS -O3 -fopenmp -std=c11 -Wall -pedantic' -outdir private .\src\mcxyz.c ".\src\libut.lib""
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Go to MATLAB's addon manager and tell it to install the "Support for MinGW-w64 compiler"
 * 2. Type "mex -setup" in the MATLAB command window and ensure that MATLAB has set the C compiler to MinGW64
 * 3. If you're using an older MATLAB version, you may need to copy the files "libgomp.a" and "libgomp.spec" to the folder with a path similar to "C:\ProgramData\MATLAB\SupportPackages\R2017a\MW_MinGW_4_9\lib\gcc\x86_64-w64-mingw32\4.9.2"
 * 4. mex should now be able to compile the code using the above command but in order to run, it needs to have the file "libgomp_64-1.dll" copied to the same folder as the mex file.
 *
 ** COMPILING ON MAC
 * As of June 2017, the macOS compiler doesn't support libut (for ctrl+c 
 * breaking) or openmp (for multithreading).
 * Compile in MATLAB with "mex COPTIMFLAGS='$COPTIMFLAGS -O3 -std=c11 -Wall -pedantic' LDOPTIMFLAGS='$LDOPTIMFLAGS -O3 -std=c11 -Wall -pedantic' -outdir private ./src/mcxyz.c"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Install XCode from the App Store
 * 2. Type "mex -setup" in the MATLAB command window
 ********************************************/

#include "mex.h"
#include <math.h>
#include <time.h>
#include "lambert.c" // For calculating the Lambert W function, originally part of the GNU Scientific Library, created by K. Briggs, G. Jungman and B. Gough and slightly modified by A. Hansen for easier mcxyz integration
#define DSFMT_MEXP 19937 // Mersenne exponent for dSFMT
#include "dSFMT-src-2.2.3/dSFMT.c" // Double precision SIMD oriented Fast Mersenne Twister(dSFMT)
#ifdef _WIN32 // This is defined on both win32 and win64 systems. We use this preprocessor condition to avoid loading openmp or libut on, e.g., Mac
#include <omp.h>
extern bool utIsInterruptPending(); // Allows catching ctrl+c while executing the mex function
#endif

#define PI          3.14159265358979323846
#define THRESHOLD   0.01		// used in roulette
#define CHANCE      0.1  		// used in roulette
#define SIGN(x)     ((x)>=0 ? 1:-1)
#define RandomNum   dsfmt_genrand_open_close(&P->dsfmt) // Calls for a random number in (0,1]
#define KILLRANGE   6.0
    /* KILLRANGE determines the region that photons are allowed to stay 
     * alive in in multiples of the cuboid size (if outside, the probability
     * of returning to the region of interest is judged as too low). When
     * launching an infinite plane wave without boundaries, photons will be
     * launched in this whole extended region. */

struct geometry { // Struct type for the constant geometry definitions
    double         d[3];
    long           n[3];
    int            boundaryFlag;
    double         *muav,*musv,*gv;
    unsigned char  *v;
};

struct beam { // Struct type for the constant beam definitions
    int            beamtypeFlag;
    double         *S;
    double         sourcesum;
    double         waist,divergence;
    double         focus[3];
    double         u0[3];
    double         v0[3];
};

struct photon { // Struct type for parameters describing the thread-specific current state of a photon
    double         i[3],u[3]; // Fractional position index i, ray trajectory unit vector u
    long           j; // Linear index of current voxel (or closest defined voxel if photon outside cuboid)
    double         mua,mus,g; // Absorption, scattering and anisotropy values at current photon position
    double         stepLeft,weight;
    bool           insideVolume,alive,sameVoxel;
    dsfmt_t        dsfmt; // "State" of the Mersenne Twister pseudo-random number generator
    long long      nLaunches; // Number of times this photon has been launched
};

void axisrotate(double const * const r, double const * const u, double const theta, double * const r_out) {
	double st = sin(theta), ct = cos(theta);
	
    r_out[0] = (u[0]*u[0]*(1-ct) +      ct)*r[0] + (u[0]*u[1]*(1-ct) - u[2]*st)*r[1] + (u[0]*u[2]*(1-ct) + u[1]*st)*r[2];
    r_out[1] = (u[1]*u[0]*(1-ct) + u[2]*st)*r[0] + (u[1]*u[1]*(1-ct) +      ct)*r[1] + (u[1]*u[2]*(1-ct) - u[0]*st)*r[2];
    r_out[2] = (u[2]*u[0]*(1-ct) - u[1]*st)*r[0] + (u[2]*u[1]*(1-ct) + u[0]*st)*r[1] + (u[2]*u[2]*(1-ct) +      ct)*r[2];
}

void launchPhoton(struct photon * const P, struct beam const * const B, struct geometry const * const G) {
    double r,phi,rand,costheta,sintheta;
    long   d,j,idx;
    double target[3],w0[3];
    
    P->alive = true;
    P->sameVoxel = false;
    P->weight = 1;
    
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
    } else switch (B->beamtypeFlag) {
        case 0: // Top-hat focus, top-hat far field beam
            phi  	= RandomNum*2*PI;
            axisrotate(B->v0,B->u0,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
            r		= B->waist*sqrt(RandomNum); // for target calculation
            for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
            phi     = RandomNum*2*PI;
            axisrotate(B->v0,B->u0,phi,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.
            phi     = B->divergence*sqrt(RandomNum); // for trajectory calculation. The sqrt is valid within paraxial approximation.
            axisrotate(B->u0,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
            P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
            P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
            P->i[2] = 0;
        break;
        case 1: // Gaussian focus, Gaussian far field beam
            phi     = RandomNum*2*PI;
            axisrotate(B->v0,B->u0,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
            r		= B->waist*sqrt(-0.5*log(RandomNum)); // for target calculation
            for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
            phi     = RandomNum*2*PI;
            axisrotate(B->v0,B->u0,phi,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.
            phi     = B->divergence*sqrt(-0.5*log(RandomNum)); // for trajectory calculation. The sqrt is valid within paraxial approximation.
            axisrotate(B->u0,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
            P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
            P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
            P->i[2] = 0;
        break;
        case 2: // Isotropically emitting point source
            P->i[0] = B->focus[0]/G->d[0] + G->n[0]/2.0;
            P->i[1] = B->focus[1]/G->d[1] + G->n[1]/2.0;
            P->i[2] = B->focus[2]/G->d[2];
            costheta = 1 - 2*RandomNum;
            sintheta = sqrt(1 - costheta*costheta);
            phi = 2*PI*RandomNum;
            P->u[0] = sintheta*cos(phi);
            P->u[1] = sintheta*sin(phi);
            P->u[2] = costheta;
        break;
        case 3: // infinite plane wave
            P->i[0] = ((G->boundaryFlag==1)? 1: KILLRANGE)*G->n[0]*(RandomNum-0.5) + G->n[0]/2.0; // Generates a random ix coordinate within the cuboid
            P->i[1] = ((G->boundaryFlag==1)? 1: KILLRANGE)*G->n[1]*(RandomNum-0.5) + G->n[1]/2.0; // Generates a random iy coordinate within the cuboid
            P->i[2] = 0;
            for(idx=0;idx<3;idx++) P->u[idx] = B->u0[idx];
        break;
        case 4: // pencil beam
            P->i[0] = (B->focus[0] - B->focus[2]*B->u0[0]/B->u0[2])/G->d[0] + G->n[0]/2.0;
            P->i[1] = (B->focus[1] - B->focus[2]*B->u0[1]/B->u0[2])/G->d[1] + G->n[1]/2.0;
            P->i[2] = 0;
            for(idx=0;idx<3;idx++) P->u[idx] = B->u0[idx];
        break;
        case 5: // top-hat focus, Gaussian far field beam
            phi     = RandomNum*2*PI;
            axisrotate(B->v0,B->u0,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
            r		= B->waist*sqrt(RandomNum); // for target calculation
            for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
            phi     = RandomNum*2*PI;
            axisrotate(B->v0,B->u0,phi,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.
            phi     = B->divergence*sqrt(-0.5*log(RandomNum)); // for trajectory calculation. The sqrt is valid within paraxial approximation.
            axisrotate(B->u0,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
            P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
            P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
            P->i[2] = 0;
        break;
        case 6: // Gaussian focus, top-hat far field beam
            phi     = RandomNum*2*PI;
            axisrotate(B->v0,B->u0,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
            r		= B->waist*sqrt(-0.5*log(RandomNum)); // for target calculation
            for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
            phi     = RandomNum*2*PI;
            axisrotate(B->v0,B->u0,phi,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.
            phi     = B->divergence*sqrt(RandomNum); // for trajectory calculation. The sqrt is valid within paraxial approximation.
            axisrotate(B->u0,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
            P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
            P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
            P->i[2] = 0;
        break;
        case 7: // Laguerre-Gaussian LG01 focus, LG01 far field beam
            phi     = RandomNum*2*PI;
            axisrotate(B->v0,B->u0,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
            r		= B->waist*sqrt((gsl_sf_lambert_Wm1(-RandomNum*exp(-1))+1)/(-2))/1.50087; // for target calculation
            for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
            phi     = B->divergence*sqrt((gsl_sf_lambert_Wm1(-RandomNum*exp(-1))+1)/(-2))/1.50087; // for trajectory calculation. The sqrt is valid within paraxial approximation.
            axisrotate(B->u0,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
            P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
            P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
            P->i[2] = 0;
        break;
    }
    
    P->stepLeft	= -log(RandomNum);
    P->nLaunches++;
}

void getNewVoxelProperties(struct photon * const P, struct geometry const * const G) {
    P->insideVolume = P->i[0] < G->n[0] && P->i[0] >= 0 &&
                      P->i[1] < G->n[1] && P->i[1] >= 0 &&
                      P->i[2] < G->n[2] && P->i[2] >= 0;

    switch (G->boundaryFlag) {
        case 0:
            P->alive = (fabs(P->i[0]/G->n[0] - 1.0/2) <  KILLRANGE/2 &&
                        fabs(P->i[1]/G->n[1] - 1.0/2) <  KILLRANGE/2 &&
                        fabs(P->i[2]/G->n[2] - 1.0/2) <  KILLRANGE/2);
            break;
        case 1:
            P->alive = P->insideVolume;
            break;
        case 2:
            P->alive = (fabs(P->i[0]/G->n[0] - 1.0/2) <  KILLRANGE/2 &&
                        fabs(P->i[1]/G->n[1] - 1.0/2) <  KILLRANGE/2 &&
                             P->i[2]/G->n[2] - 1.0/2  <  KILLRANGE/2 &&
                             P->i[2]                  >= 0);
            break;
    }

    /* Get tissue voxel properties of current position.
     * If photon is outside cuboid, properties are those of
     * the closest defined voxel. */
    P->j = ((P->i[2] < 0)? 0: ((P->i[2] >= G->n[2])? G->n[2]-1: floor(P->i[2])))*G->n[0]*G->n[1] +
           ((P->i[1] < 0)? 0: ((P->i[1] >= G->n[1])? G->n[1]-1: floor(P->i[1])))*G->n[0]         +
           ((P->i[0] < 0)? 0: ((P->i[0] >= G->n[0])? G->n[0]-1: floor(P->i[0]))); // Index values are restrained to integers in the interval [0,n-1]
    P->mua = G->muav[G->v[P->j]];
    P->mus = G->musv[G->v[P->j]];
    P->g   = G->gv  [G->v[P->j]];
}

void propagatePhoton(struct photon * const P, struct geometry const * const G, double * const F) {
    long idx;

    long double s = (long double)P->stepLeft/P->mus; // s starts out as the propagation distance until next scattering event (assuming current mus value)...
    for(idx=0;idx<3;idx++) if(P->u[idx]) s = fminl(s,(floor(P->i[idx]) + (P->u[idx]>0) - (long double)P->i[idx])*G->d[idx]/P->u[idx]); // but if any of the distances to the next voxel boundary (yz, xz or xy) are smaller, set s to that smallest value

    for(idx=0;idx<3;idx++) {
        P->i[idx] += s*P->u[idx]/G->d[idx]; // Here it is important that s (a long double) has higher precision than P->i (a double)
        if(!fmod(P->i[idx],1) && P->u[idx]) { // If we have arrived at a voxel boundary without simply traveling parallel to it...
            P->sameVoxel = false; // then we are in a new voxel...
            if(P->u[idx] < 0) P->i[idx] = nextafter(P->i[idx],-INFINITY); // and if we were traveling in the negative direction, we need to nudge the photon the smallest possible amount into the next voxel
        }
    }
    P->stepLeft -= s*P->mus;
    
    double absorb = P->weight*(1 - exp(-P->mua*s));   // photon weight absorbed at this step
    P->weight -= absorb;					   // decrement WEIGHT by amount absorbed
    if(P->insideVolume) {	// only save data if the photon is inside simulation cuboid
        #ifdef _WIN32
        #pragma omp atomic
        #endif
            F[P->j] += absorb;
    }
}

void checkRoulette(struct photon * const P) {
    /**** CHECK ROULETTE 
     If photon weight below THRESHOLD, then terminate photon using Roulette technique.
     Photon has CHANCE probability of having its weight increased by factor of 1/CHANCE,
     and 1-CHANCE probability of terminating. */
    if(P->weight < THRESHOLD) {
        if(RandomNum <= CHANCE) P->weight /= CHANCE;
        else P->alive = false;
    }
}

void scatterPhoton(struct photon * const P, struct geometry const * const G) {
    // Sample for costheta using Henyey-Greenstein scattering
    double costheta = P->g? (1 + P->g*P->g - pow((1 - P->g*P->g)/(1 - P->g + 2*P->g*RandomNum),2))/(2*P->g) : 2*RandomNum - 1;
    double sintheta = sqrt(1 - costheta*costheta);
    double phi = 2*PI*RandomNum;
    double cosphi = cos(phi);
    double sinphi = sin(phi);

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
    
    P->stepLeft	= -log(RandomNum);
}

void normalizeDeposition(struct beam const * const B, struct geometry const * const G, long long nPhotons, double * const F) {
    long j;
    long L = G->n[0]*G->n[1]*G->n[2]; // Total number of voxels in cuboid
    // Normalize deposition to yield either absolute or relative fluence rate (F).
    if(B->S) { // For a 3D source distribution (e.g., fluorescence)
		for(j=0; j<L;j++) F[j] /= nPhotons*G->muav[G->v[j]]/B->sourcesum; // Absolute fluence rate [W/cm^2]
        mxFree(B->S);
    } else if(B->beamtypeFlag == 3 && G->boundaryFlag != 1) { // For infinite plane wave launched into volume without absorbing walls
		for(j=0; j<L;j++) F[j] /= G->d[0]*G->d[1]*G->d[2]*nPhotons*G->muav[G->v[j]]/(KILLRANGE*KILLRANGE); // Normalized fluence rate [W/cm^2/W.incident]
	} else {
		for(j=0; j<L;j++) F[j] /= G->d[0]*G->d[1]*G->d[2]*nPhotons*G->muav[G->v[j]]; // Normalized fluence rate [W/cm^2/W.incident]
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    long      idx;                  // General-purpose non-thread-specific index variable
	int       pctProgress = 0;      // Simulation progress in percent
    int       newPctProgress;    
    bool      ctrlc_caught = false; // Has a ctrl+c been passed from MATLAB?
    long long nPhotons = 0;
    
    // Timekeeping variables
    double          simulationTimeRequested = *mxGetPr(mxGetField(prhs[0],0,"simulationTime"));
	double          simulationTimeSpent;
	struct timespec simulationTimeStart, simulationTimeCurrent;
	
    // Geometry struct definition
    mxArray *tissueList = mxGetField(prhs[0],0,"tissueList");
    int     nT = mxGetN(tissueList);
    double 	muav[nT];            // muav[0:nT-1], absorption coefficient of ith tissue type
	double 	musv[nT];            // scattering coeff. 
	double 	gv[nT];              // anisotropy of scattering
    for(idx=0;idx<nT;idx++) {
        muav[idx] = *mxGetPr(mxGetField(tissueList,idx,"mua"));
        musv[idx] = *mxGetPr(mxGetField(tissueList,idx,"mus"));
        gv[idx]   = *mxGetPr(mxGetField(tissueList,idx,"g"));
    }
    
    mwSize const  *dimPtr = mxGetDimensions(mxGetField(prhs[0],0,"T"));
    struct geometry const G_var = (struct geometry) {
        .d = {*mxGetPr(mxGetField(prhs[0],0,"dx")),
              *mxGetPr(mxGetField(prhs[0],0,"dy")),
              *mxGetPr(mxGetField(prhs[0],0,"dz"))},
        .n = {dimPtr[0], dimPtr[1], dimPtr[2]},
        .boundaryFlag = *mxGetPr(mxGetField(prhs[0],0,"boundaryFlag")),
        .muav = muav,
        .musv = musv,
        .gv = gv,
        .v = mxGetData(mxGetField(prhs[0],0,"T"))
    };
    struct geometry const *G = &G_var;
    
    long L = G->n[0]*G->n[1]*G->n[2]; // Total number of voxels in cuboid
    
	// Beam struct definition
    double *S_PDF       = mxGetData(mxGetField(prhs[0],0,"sourceDistribution")); // Probability density function of 3D source distribution
    double sourcesum    = 0;
    double *S           = NULL; // Cumulative distribution function
	if(S_PDF) {
        S = mxMalloc((L+1)*sizeof(double));
        if(!S) {printf("Out of memory, returning...\n"); return;}
        S[0] = 0;
        for(idx=1;idx<(L+1);idx++) S[idx] = S[idx-1] + S_PDF[idx-1];
        sourcesum = S[L];
        for(idx=1;idx<(L+1);idx++) S[idx] /= sourcesum;
    }
    
    double u0[3] =      {*mxGetPr(mxGetField(prhs[0],0,"ux0")),
                         *mxGetPr(mxGetField(prhs[0],0,"uy0")),
                         *mxGetPr(mxGetField(prhs[0],0,"uz0"))};
    
    struct beam const B_var = (struct beam) {
        .beamtypeFlag =  *mxGetPr(mxGetField(prhs[0],0,"beamtypeFlag")),
        .S            =  S,
        .sourcesum    =  sourcesum,
        .waist        =  *mxGetPr(mxGetField(prhs[0],0,"waist")),
        .divergence   =  *mxGetPr(mxGetField(prhs[0],0,"divergence")),
        .focus        = {*mxGetPr(mxGetField(prhs[0],0,"xFocus")),
                         *mxGetPr(mxGetField(prhs[0],0,"yFocus")),
                         *mxGetPr(mxGetField(prhs[0],0,"zFocus"))},
        .u0           = {u0[0],u0[1],u0[2]},
        .v0           = {(u0[2]==1)? 1: -u0[1]/sqrt(u0[1]*u0[1] + u0[0]*u0[0]),
                         (u0[2]==1)? 0:  u0[0]/sqrt(u0[1]*u0[1] + u0[0]*u0[0]),
                         0} // normal vector to photon trajectory
    };
    struct beam const *B = &B_var;
    
    // Prepare output matrix
    plhs[0] = mxCreateNumericArray(3,dimPtr,mxDOUBLE_CLASS,mxREAL);
    double  *F = mxGetPr(plhs[0]);
    for(idx=0;idx<L;idx++) F[idx] = 0;
    
    // Display parameters
    printf("---------------------------------------------------------\n");
    printf("(nx,ny,nz) = (%d,%d,%d), (dx,dy,dz) = (%0.1e,%0.1e,%0.1e) cm\n",G->n[0],G->n[1],G->n[2],G->d[0],G->d[1],G->d[2]);
    
    if(B->S) printf("Using isotropically emitting 3D distribution as light source\n");
    else {
        switch(B->beamtypeFlag) {
            case 0:
                printf("Using top-hat focus, top-hat far field beam\n");
                break;
            case 1:
                printf("Using Gaussian focus, Gaussian far field beam\n");
                break;
            case 2:
                printf("Using isotropic point source\n");
                break;
            case 3:
                printf("Using infinite plane wave\n");
                break;
            case 4:
                printf("Using pencil beam\n");
                break;
            case 5:
                printf("Using top-hat focus, Gaussian far field beam\n");
                break;
            case 6:
                printf("Using Gaussian focus, top-hat far field beam\n");
                break;
            case 7:
                printf("Using Laguerre-Gaussian LG01 focus, LG01 far field beam\n");
                break;
        }

        printf("Focus location = (%0.1e,%0.1e,%0.1e) cm\n",B->focus[0],B->focus[1],B->focus[2]);
        if(B->beamtypeFlag!=2) printf("Beam direction = (%0.3f,%0.3f,%0.3f)\n",B->u0[0],B->u0[1],B->u0[2]);
        if(B->beamtypeFlag!=2 && B->beamtypeFlag!=3 && B->beamtypeFlag!=4) {
            printf("Beam waist = %0.1e cm, divergence = %0.1e rad\n",B->waist,B->divergence);
        }
    }
    
    if(G->boundaryFlag==0) printf("boundaryFlag = 0, so no boundaries\n");
    else if(G->boundaryFlag==1) printf("boundaryFlag = 1, so escape at all boundaries\n");    
	else if(G->boundaryFlag==2) printf("boundaryFlag = 2, so escape at surface only\n");    
	else {
        printf("Error: Invalid boundaryFlag\n");
        return;
    }
    
    printf("Simulation duration = %0.2f min\n\n",simulationTimeRequested);
	printf("Calculating...   0%% done");
    mexEvalString("drawnow;");
	
    // ============================ MAJOR CYCLE ========================
	clock_gettime(CLOCK_MONOTONIC, &simulationTimeStart);
	
    #ifdef _WIN32
    #pragma omp parallel num_threads(omp_get_num_procs())
    #endif
	{
        struct photon P_var;
        struct photon *P = &P_var;
        
        P->nLaunches = 0;
        
        #ifdef _WIN32
		dsfmt_init_gen_rand(&P->dsfmt,simulationTimeStart.tv_nsec + omp_get_thread_num()); // Seed the photon's random number generator
        #else
		dsfmt_init_gen_rand(&P->dsfmt,simulationTimeStart.tv_nsec); // Seed the photon's random number generator
        #endif
        
		while((pctProgress < 100) && !ctrlc_caught) {
            launchPhoton(P,B,G);
            
			while(P->alive) {
				while(P->stepLeft>0 && P->alive) {
					if(!P->sameVoxel) getNewVoxelProperties(P,G); // If photon has just entered a new voxel or has just been launched
                    if(P->alive) propagatePhoton(P,G,F);
				}
                if(P->alive) checkRoulette(P);
				if(P->alive) scatterPhoton(P,G);
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
                newPctProgress = 100*(simulationTimeCurrent.tv_sec  - simulationTimeStart.tv_sec +
                                     (simulationTimeCurrent.tv_nsec - simulationTimeStart.tv_nsec)/1e9) / (simulationTimeRequested*60);
				if((newPctProgress != pctProgress) && !ctrlc_caught) {
					pctProgress = newPctProgress;
                    printf("\b\b\b\b\b\b\b\b\b%3.i%% done", pctProgress);
                    mexEvalString("drawnow;");
				}
			}
		}
        #ifdef _WIN32
        #pragma omp atomic
        #endif
        nPhotons += P->nLaunches; // Add up the launches made in the individual threads to the total photon counter
	}
    
	clock_gettime(CLOCK_MONOTONIC, &simulationTimeCurrent);
	simulationTimeSpent = simulationTimeCurrent.tv_sec  - simulationTimeStart.tv_sec +
                         (simulationTimeCurrent.tv_nsec - simulationTimeStart.tv_nsec)/1e9;
	printf("\nSimulated %0.2e photons at a rate of %0.2e photons per minute\n",(double)nPhotons, nPhotons*60/simulationTimeSpent);
    printf("---------------------------------------------------------\n");
    mexEvalString("drawnow;");
    normalizeDeposition(B,G,nPhotons,F); // Convert data to either relative or absolute fluence rate
}
