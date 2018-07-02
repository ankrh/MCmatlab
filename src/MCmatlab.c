/********************************************
 *
 *  MCmatlab.c,	in the C programming language, written for MATLAB MEX function generation
 *  Compliant with the ISO C11 standard
 *
 *  Heavily inspired by mcxyz.c by Steven Jacques, Ting Li, Scott Prahl at the Oregon Health & Science University
 *
 ** COMPILING ON WINDOWS
 * Can be compiled in MATLAB with "mex COPTIMFLAGS='$COPTIMFLAGS -O3 -fopenmp -std=c11 -Wall -pedantic' LDOPTIMFLAGS='$LDOPTIMFLAGS -O3 -fopenmp -std=c11 -Wall -pedantic' -outdir private .\src\MCmatlab.c ".\src\libut.lib""
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
 * Compile in MATLAB with "mex COPTIMFLAGS='$COPTIMFLAGS -O3 -std=c11 -Wall -pedantic' LDOPTIMFLAGS='$LDOPTIMFLAGS -O3 -std=c11 -Wall -pedantic' -outdir private ./src/MCmatlab.c"
 *
 * To get the MATLAB C compiler to work, try this:
 * 1. Install XCode from the App Store
 * 2. Type "mex -setup" in the MATLAB command window
 ********************************************/
    #include <windows.h>


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
    double         power;
    double         waist,divergence;
    double         focus[3];
    double         u[3];
    double         v[3];
};

struct lightcollector { // Struct type for the constant light collector definitions. It can be either an objective lens (for 0<f<INFINITY) or a fiber tip or simple aperture (for f=INFINITY)
    double         r[3]; // Position of center of objective focal plane (not the objective itself) or position of center of the fiber tip
    double         theta; // Polar angle that the objective or fiber is facing
    double         phi; // Azimuthal angle of objective or fiber orientation
    double         f; // Focal length of the objective. If the light collector is a fiber tip, this will be INFINITY.
    double         diam; // Diameter of the objective aperture or core diameter of the fiber. For an ideal thin lens objective, this is 2*tan(arcsin(lensNA/f)).
    double         FSorNA; // For an objective lens: Field Size of the imaging system (diameter of area in object plane that gets imaged). For a fiber tip: The fiber's NA.
    long           res[2]; // Resolution of light collector planes in pixels along X and Y axes
};

struct photon { // Struct type for parameters describing the thread-specific current state of a photon
    double         i[3],u[3],D[3]; // Fractional position indices i, ray trajectory unit vector u and distances D to next voxel boundary (yz, xz or yz) along current trajectory
    long           j; // Linear index of current voxel (or closest defined voxel if photon outside cuboid)
    double         mua,mus,g; // Absorption, scattering and anisotropy values at current photon position
    double         stepLeft,weight;
    bool           insideVolume,alive,sameVoxel;
    dsfmt_t        dsfmt; // "State" of the Mersenne Twister pseudo-random number generator
    long long      nThreadPhotons; // Number of times this thread has launched a photon
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
        case 0: // pencil beam
            P->i[0] = (B->focus[0] - B->focus[2]*B->u[0]/B->u[2])/G->d[0] + G->n[0]/2.0;
            P->i[1] = (B->focus[1] - B->focus[2]*B->u[1]/B->u[2])/G->d[1] + G->n[1]/2.0;
            P->i[2] = 0;
            for(idx=0;idx<3;idx++) P->u[idx] = B->u[idx];
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
        break;
        case 2: // infinite plane wave
            P->i[0] = ((G->boundaryFlag==1)? 1: KILLRANGE)*G->n[0]*(RandomNum-0.5) + G->n[0]/2.0; // Generates a random ix coordinate within the cuboid
            P->i[1] = ((G->boundaryFlag==1)? 1: KILLRANGE)*G->n[1]*(RandomNum-0.5) + G->n[1]/2.0; // Generates a random iy coordinate within the cuboid
            P->i[2] = 0;
            for(idx=0;idx<3;idx++) P->u[idx] = B->u[idx];
        break;
        case 3: // Gaussian focus, Gaussian far field beam
            phi     = RandomNum*2*PI;
            axisrotate(B->v,B->u,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
            r		= B->waist*sqrt(-0.5*log(RandomNum)); // for target calculation
            for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
            phi     = RandomNum*2*PI;
            axisrotate(B->v,B->u,phi,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.
            phi     = B->divergence*sqrt(-0.5*log(RandomNum)); // for trajectory calculation. The sqrt is valid within paraxial approximation.
            axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
            P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
            P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
            P->i[2] = 0;
        break;
        case 4: // Gaussian focus, top-hat far field beam
            phi     = RandomNum*2*PI;
            axisrotate(B->v,B->u,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
            r		= B->waist*sqrt(-0.5*log(RandomNum)); // for target calculation
            for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
            phi     = RandomNum*2*PI;
            axisrotate(B->v,B->u,phi,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.
            phi     = B->divergence*sqrt(RandomNum); // for trajectory calculation. The sqrt is valid within paraxial approximation.
            axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
            P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
            P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
            P->i[2] = 0;
        break;
        case 5: // top-hat focus, Gaussian far field beam
            phi     = RandomNum*2*PI;
            axisrotate(B->v,B->u,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
            r		= B->waist*sqrt(RandomNum); // for target calculation
            for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
            phi     = RandomNum*2*PI;
            axisrotate(B->v,B->u,phi,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.
            phi     = B->divergence*sqrt(-0.5*log(RandomNum)); // for trajectory calculation. The sqrt is valid within paraxial approximation.
            axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
            P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
            P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
            P->i[2] = 0;
        break;
        case 6: // top-hat focus, top-hat far field beam
            phi  	= RandomNum*2*PI;
            axisrotate(B->v,B->u,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
            r		= B->waist*sqrt(RandomNum); // for target calculation
            for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
            phi     = RandomNum*2*PI;
            axisrotate(B->v,B->u,phi,w0); // w0 unit vector is now normal to both beam center axis and to ray propagation direction. Angle from v0 to w0 is phi.
            phi     = B->divergence*sqrt(RandomNum); // for trajectory calculation. The sqrt is valid within paraxial approximation.
            axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
            P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
            P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
            P->i[2] = 0;
        break;
        case 7: // Laguerre-Gaussian LG01 beam
            phi     = RandomNum*2*PI;
            axisrotate(B->v,B->u,phi,w0); // w0 unit vector now points in the direction from focus center point to ray target point
            r		= B->waist*sqrt((gsl_sf_lambert_Wm1(-RandomNum*exp(-1))+1)/(-2))/1.50087; // for target calculation
            for(idx=0;idx<3;idx++) target[idx] = B->focus[idx] + r*w0[idx];
            phi     = B->divergence*sqrt((gsl_sf_lambert_Wm1(-RandomNum*exp(-1))+1)/(-2))/1.50087; // for trajectory calculation. The sqrt is valid within paraxial approximation.
            axisrotate(B->u,w0,phi,P->u); // ray propagation direction is found by rotating beam center axis an angle phi around w0
            P->i[0] = (target[0] - target[2]*P->u[0]/P->u[2])/G->d[0] + G->n[0]/2.0; // the coordinates for the ray starting point is the intersection of the ray with the z = 0 surface
            P->i[1] = (target[1] - target[2]*P->u[1]/P->u[2])/G->d[1] + G->n[1]/2.0;
            P->i[2] = 0;
        break;
    }
    
    // Calculate distances to next voxel boundary planes
    for(idx=0;idx<3;idx++) P->D[idx] = P->u[idx]? (floor(P->i[idx]) + (P->u[idx]>0) - P->i[idx])*G->d[idx]/P->u[idx] : INFINITY;
    
    P->stepLeft	= -log(RandomNum);
    P->nThreadPhotons++;
}

void checkEscape(struct photon * const P, struct geometry const * const G, struct lightcollector const * const LC, double *LCP, double *ImP) {
    bool escaped = false;
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
            escaped = !P->insideVolume;
            break;
        case 2:
            P->alive = (fabs(P->i[0]/G->n[0] - 1.0/2) <  KILLRANGE/2 &&
                        fabs(P->i[1]/G->n[1] - 1.0/2) <  KILLRANGE/2 &&
                             P->i[2]/G->n[2] - 1.0/2  <  KILLRANGE/2 &&
                             P->i[2]                  >= 0);
            escaped = (P->i[2] < 0);
            break;
    }
    
    if(escaped && LCP) { // If LCP is not NULL then that's because useLightCollector was set to true (non-zero)
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
                        LCP [(long)(LC->res[0]*(RLCP[0]/LC->diam   + 1.0/2)) + LC->res[0]*(long)(LC->res[1]*(RLCP[1]/LC->diam   + 1.0/2))] += P->weight;
                        ImP [(long)(LC->res[0]*(RImP[0]/LC->FSorNA + 1.0/2)) + LC->res[0]*(long)(LC->res[1]*(RImP[1]/LC->FSorNA + 1.0/2))] += P->weight;
//                         double thetaLCFF = atan(distImP/LC->f); // Light collector far field polar angle
//                         double phiLCFF = atan2(RImP[1],RImP[0]); // Light collector far field azimuthal angle
//                         LCFF[(long)(LC->res[0]*(thetaLCFF/(PI/2*(1+DBL_EPSILON)))) + LC->res[0]*(long)(LC->res[1]*((phiLCFF+PI)/(2*PI*(1+DBL_EPSILON))))] += P->weight;
                    }
                } else { // If the light collector is a fiber tip
                    double thetaLCFF = atan(-sqrt(U[0]*U[0] + U[1]*U[1])/U[2]); // Light collector far field polar angle
                    if(thetaLCFF < asin(fmin(1,LC->FSorNA))) { // If the photon has an angle within the fiber's NA acceptance
                        LCP [(long)(LC->res[0]*(RLCP[0]/LC->diam   + 1.0/2)) + LC->res[0]*(long)(LC->res[1]*(RLCP[1]/LC->diam   + 1.0/2))] += P->weight;
//                         double phiLCFF = atan2(U[1],U[0]); // Light collector far field azimuthal angle
//                         LCFF[(long)(LC->res[0]*(thetaLCFF/(PI/2*(1+DBL_EPSILON)))) + LC->res[0]*(long)(LC->res[1]*((phiLCFF+PI)/(2*PI*(1+DBL_EPSILON))))] += P->weight;
                    }
                }
            }
        }
    }
}

void getNewVoxelProperties(struct photon * const P, struct geometry const * const G) {
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

    P->sameVoxel = true;

    double s = fmin(P->stepLeft/P->mus,fmin(P->D[0],fmin(P->D[1],P->D[2])));
    P->stepLeft -= s*P->mus;
    
    for(idx=0;idx<3;idx++) {
        long i_old = floor(P->i[idx]);
        if(s == P->D[idx]) { // If we're supposed to go to the next voxel along this dimension
            P->i[idx] = (P->u[idx] > 0)? i_old + 1: i_old - DBL_EPSILON*(labs(i_old)+1);
            P->sameVoxel = false;
            P->D[idx] = G->d[idx]/fabs(P->u[idx]); // Reset voxel boundary distance
        } else { // If we were supposed to remain in the same voxel along this dimension
            P->i[idx] += s*P->u[idx]/G->d[idx]; // First take the expected step (including various rounding errors)
            if(floor(P->i[idx]) != i_old) P->i[idx] = (P->u[idx] > 0)? i_old + 1 - DBL_EPSILON*(labs(i_old)+1): i_old ; // If photon due to rounding errors actually crossed the border, set it to be barely in the original voxel
            P->D[idx] -= s;
        }
    }
    
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

    P->stepLeft	= -log(RandomNum);
}

void normalizeDeposition(struct beam const * const B, struct geometry const * const G, struct lightcollector const * const LC, double * const F, double * const LCP, double * const ImP, double nPhotons) {
    long j;
    double V = G->d[0]*G->d[1]*G->d[2]; // Voxel volume
    long L = G->n[0]*G->n[1]*G->n[2]; // Total number of voxels in cuboid
    long L_LC = LC->res[0]*LC->res[1]; // Total number of pixels in light collector planes
    // Normalize deposition to yield either absolute or relative fluence rate (F).
    if(B->S) { // For a 3D source distribution (e.g., fluorescence)
		for(j=0; j<L;j++) F[j] /= V*nPhotons*G->muav[G->v[j]]/B->power; // Absolute fluence rate [W/cm^2].
        if(LCP) for(j=0;j<L_LC;j++) {
            LCP[j] /= LC->diam  *LC->diam  /L_LC*nPhotons/B->power; // Absolute fluence rate [W/cm^2]
            ImP[j] /= LC->FSorNA*LC->FSorNA/L_LC*nPhotons/B->power; // Absolute fluence rate [W/cm^2]
        }
        mxFree(B->S);
    } else if(B->beamtypeFlag == 3 && G->boundaryFlag != 1) { // For infinite plane wave launched into volume without absorbing walls
		for(j=0; j<L;j++) F[j] /= V*nPhotons*G->muav[G->v[j]]/(KILLRANGE*KILLRANGE); // Normalized fluence rate [W/cm^2/W.incident]
        if(LCP) for(j=0;j<L_LC;j++) {
            LCP[j] /= LC->diam  *LC->diam  /L_LC*nPhotons/(KILLRANGE*KILLRANGE); // Normalized fluence rate [W/cm^2/W.incident]
            ImP[j] /= LC->FSorNA*LC->FSorNA/L_LC*nPhotons/(KILLRANGE*KILLRANGE); // Normalized fluence rate [W/cm^2/W.incident]
        }
	} else {
		for(j=0; j<L;j++) F[j] /= V*nPhotons*G->muav[G->v[j]]; // Normalized fluence rate [W/cm^2/W.incident]
        if(LCP) for(j=0;j<L_LC;j++) {
            LCP[j] /= LC->diam  *LC->diam  /L_LC*nPhotons; // Normalized fluence rate [W/cm^2/W.incident]
            ImP[j] /= LC->FSorNA*LC->FSorNA/L_LC*nPhotons; // Normalized fluence rate [W/cm^2/W.incident]
        }
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[]) {
    long      idx;                  // General-purpose non-thread-specific index variable
	int       pctProgress = 0;      // Simulation progress in percent
    int       newPctProgress;    
    bool      ctrlc_caught = false; // Has a ctrl+c been passed from MATLAB?
    
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
    double *S_PDF       = mxGetData(mxGetField(prhs[0],0,"sourceDistribution")); // Power emitted by the individual voxels per unit volume. Can be percieved as an unnormalized probability density function of the 3D source distribution
    double power        = 0;
    double *S           = NULL; // Cumulative distribution function
	if(S_PDF) {
        S = mxMalloc((L+1)*sizeof(double));
        if(!S) {printf("Out of memory, returning...\n"); return;}
        S[0] = 0;
        for(idx=1;idx<(L+1);idx++) S[idx] = S[idx-1] + S_PDF[idx-1];
        power = S[L]*G->d[0]*G->d[1]*G->d[2];
        for(idx=1;idx<(L+1);idx++) S[idx] /= S[L];
    }
    
    double tb = *mxGetPr(mxGetField(prhs[0],0,"thetaBeam"));
    double pb = *mxGetPr(mxGetField(prhs[0],0,"phiBeam"));

    double u[3] = {sin(tb)*cos(pb),
                   sin(tb)*sin(pb),
                   cos(tb)}; // Temporary array
    
    double v[3] = {1,0,0};
    if(u[2]!=1) unitcrossprod(u,(double[]){0,0,1},v);
    
    struct beam const B_var = (struct beam) {
        .beamtypeFlag =  *mxGetPr(mxGetField(prhs[0],0,"beamtypeFlag")),
        .S            =  S,
        .power        =  power,
        .waist        =  *mxGetPr(mxGetField(prhs[0],0,"waist")),
        .divergence   =  *mxGetPr(mxGetField(prhs[0],0,"divergence")),
        .focus        = {*mxGetPr(mxGetField(prhs[0],0,"xFocus")),
                         *mxGetPr(mxGetField(prhs[0],0,"yFocus")),
                         *mxGetPr(mxGetField(prhs[0],0,"zFocus"))},
        .u            = {u[0],u[1],u[2]},
        .v            = {v[0],v[1],v[2]} // normal vector to beam center axis
    };
    struct beam const *B = &B_var;
    
    // Detector struct definition
    bool useLightCollector = *mxGetPr(mxGetField(prhs[0],0,"useLightCollector"));
    
    double theta = *mxGetPr(mxGetField(prhs[0],0,"theta_LC"));
    double phi   = *mxGetPr(mxGetField(prhs[0],0,"phi_LC"));
    double f     = *mxGetPr(mxGetField(prhs[0],0,"f_LC"));

    struct lightcollector const LC_var = (struct lightcollector) {
        .r            = {*mxGetPr(mxGetField(prhs[0],0,"xFPC_LC")) - (isfinite(f)? f*sin(theta)*cos(phi):0), // xyz coordinates of center of light collector
                         *mxGetPr(mxGetField(prhs[0],0,"yFPC_LC")) - (isfinite(f)? f*sin(theta)*sin(phi):0),
                         *mxGetPr(mxGetField(prhs[0],0,"zFPC_LC")) - (isfinite(f)? f*cos(theta):0)},
        .theta        =  theta,
        .phi          =  phi,
        .f            =  f,
        .diam         =  *mxGetPr(mxGetField(prhs[0],0,"diam_LC")),
        .FSorNA       =  *mxGetPr(mxGetField(prhs[0],0,"FSorNA_LC")),
        .res          = {*mxGetPr(mxGetField(prhs[0],0,"resX_LC")),
                         *mxGetPr(mxGetField(prhs[0],0,"resY_LC"))}};
    struct lightcollector const *LC = &LC_var;
    
    // Prepare output MATLAB struct
    double *LCP  = NULL;
    double *ImP  = NULL;
//     double *LCFF = NULL;
    
    if(useLightCollector) {
        plhs[0] = mxCreateStructMatrix(1,1,5,(const char *[]){"F","LCP","ImP","nPhotons","nThreads"});
        mxSetField(plhs[0],0,"LCP", mxCreateDoubleMatrix(LC->res[0],LC->res[1],mxREAL));
        mxSetField(plhs[0],0,"ImP", mxCreateDoubleMatrix(LC->res[0],LC->res[1],mxREAL));
//         mxSetField(plhs[0],0,"LCFF",mxCreateDoubleMatrix(LC->res[0],LC->res[1],mxREAL));
        LCP  = mxGetPr(mxGetField(plhs[0],0,"LCP"));
        ImP  = mxGetPr(mxGetField(plhs[0],0,"ImP"));
//         LCFF = mxGetPr(mxGetField(plhs[0],0,"LCFF"));
    } else {
        plhs[0] = mxCreateStructMatrix(1,1,3,(const char *[]){"F","nPhotons","nThreads"});
    }
    mxSetField(plhs[0],0,"F",       mxCreateNumericArray(3,dimPtr,mxDOUBLE_CLASS,mxREAL));
    mxSetField(plhs[0],0,"nPhotons",mxCreateDoubleMatrix(1,1,mxREAL));
    mxSetField(plhs[0],0,"nThreads",mxCreateDoubleMatrix(1,1,mxREAL));
    double *F           = mxGetPr(mxGetField(plhs[0],0,"F"));
    double *nPhotonsPtr = mxGetPr(mxGetField(plhs[0],0,"nPhotons"));
    double *nThreadsPtr = mxGetPr(mxGetField(plhs[0],0,"nThreads"));
    
    // Display parameters
    printf("---------------------------------------------------------\n");
    printf("(nx,ny,nz) = (%d,%d,%d), (dx,dy,dz) = (%0.1e,%0.1e,%0.1e) cm\n",G->n[0],G->n[1],G->n[2],G->d[0],G->d[1],G->d[2]);
    
    if(B->S) printf("Using isotropically emitting 3D distribution as light source\n");
    else {
        switch(B->beamtypeFlag) {
            case 0:
                printf("Using pencil beam\n");
                break;
            case 1:
                printf("Using isotropically emitting point source\n");
                break;
            case 2:
                printf("Using infinite plane wave\n");
                break;
            case 3:
                printf("Using Gaussian focus, Gaussian far field beam\n");
                break;
            case 4:
                printf("Using Gaussian focus, top-hat far field beam\n");
                break;
            case 5:
                printf("Using top-hat focus, Gaussian far field beam\n");
                break;
            case 6:
                printf("Using top-hat focus, top-hat far field beam\n");
                break;
            case 7:
                printf("Using Laguerre-Gaussian LG01 beam\n");
                break;
        }

        if(B->beamtypeFlag!=2) printf("Focus location = (%0.1e,%0.1e,%0.1e) cm\n",B->focus[0],B->focus[1],B->focus[2]);
        if(B->beamtypeFlag!=1) printf("Beam direction = (%0.3f,%0.3f,%0.3f)\n",B->u[0],B->u[1],B->u[2]);
        if(B->beamtypeFlag >2) printf("Beam waist = %0.1e cm, divergence = %0.1e rad\n",B->waist,B->divergence);
    }
    
    if(G->boundaryFlag==0) printf("boundaryFlag = 0, so no boundaries\n");
    else if(G->boundaryFlag==1) printf("boundaryFlag = 1, so escape at all boundaries\n");    
	else if(G->boundaryFlag==2) printf("boundaryFlag = 2, so escape at surface only\n");    
    
    printf("Simulation duration = %0.2f min\n\n",simulationTimeRequested);
	printf("Calculating...   0%% done");
    mexEvalString("drawnow;");
	
    // ============================ MAJOR CYCLE ========================
	clock_gettime(CLOCK_MONOTONIC, &simulationTimeStart);
	
    #ifdef _WIN32
    *nThreadsPtr = omp_get_num_procs();
    #pragma omp parallel num_threads((long)*nThreadsPtr)
    #else
    *nThreadsPtr = 1;
    #endif
	{
        struct photon P_var;
        struct photon *P = &P_var;
        
        P->nThreadPhotons = 0;
        
        #ifdef _WIN32
		dsfmt_init_gen_rand(&P->dsfmt,simulationTimeStart.tv_nsec + omp_get_thread_num()); // Seed the photon's random number generator
        #else
		dsfmt_init_gen_rand(&P->dsfmt,simulationTimeStart.tv_nsec); // Seed the photon's random number generator
        #endif
        
		while((pctProgress < 100) && !ctrlc_caught) {
            launchPhoton(P,B,G);
            
			while(P->alive) {
				while(P->stepLeft>0) {
					if(!P->sameVoxel) {
                        checkEscape(P,G,LC,LCP,ImP);
                        if(!P->alive) break;
                        getNewVoxelProperties(P,G); // If photon has just entered a new voxel or has just been launched
                    }
                    propagatePhoton(P,G,F);
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
        *nPhotonsPtr += P->nThreadPhotons; // Add up the launches made in the individual threads to the total photon counter
	}
    
	clock_gettime(CLOCK_MONOTONIC, &simulationTimeCurrent);
	simulationTimeSpent = simulationTimeCurrent.tv_sec  - simulationTimeStart.tv_sec +
                         (simulationTimeCurrent.tv_nsec - simulationTimeStart.tv_nsec)/1e9;
	printf("\nSimulated %0.2e photons at a rate of %0.2e photons per minute\n",*nPhotonsPtr, *nPhotonsPtr*60/simulationTimeSpent);
    printf("---------------------------------------------------------\n");
    mexEvalString("drawnow;");
    normalizeDeposition(B,G,LC,F,LCP,ImP,*nPhotonsPtr); // Convert data to either relative or absolute fluence rate
}
